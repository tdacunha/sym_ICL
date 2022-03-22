import lib
import orbit_model
import struct
import array
import numpy as np
import galpy.potential as potential
import scipy.stats as stats
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import palette
import numpy.random as random
from palette import pc
from colossus.cosmology import cosmology

params = { 'flat': True, 'H0': 70.0, 'Om0': 0.286,
           'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.95 }
cosmo = cosmology.setCosmology('myCosmo', params)
    
P_FILE_FMT = "%s/particles/part_%03d.%d"
MP = 2.8e5

TAG_METHOD = "energy"
#TAG_METHOD = "infall"

def read_part_file(base_dir, snap, i, vars_to_read=["x", "v", "phi"]):
    """ read_part_file reads the particles from a file for halo i at the given
    snapshot. You can change which particles are read using vars_to_read. By
    default, all three are read and returned as a tuple of (x, v, phi).
    """
    fname = P_FILE_FMT % (base_dir, snap, i)
    f = open(fname, "rb")

    n = struct.unpack("i", f.read(4))[0]
    x, v, phi = array.array("f"), array.array("f"), array.array("f")

    if "x" in vars_to_read:
        x.fromfile(f, n*3)
        x = np.array(x, dtype=float)
        x = x.reshape((n, 3))
    else:
        x = None
        f.seek(12*n, 1)

    if "v" in vars_to_read:
        v.fromfile(f, n*3)
        v = np.array(v, dtype=float)
        v = v.reshape((n, 3))
    else:
        v = None
        f.seek(12*n, 1)

    if "phi" in vars_to_read:
        phi.fromfile(f, n)
        phi = np.array(phi, dtype=float)
    else:
        phi = None

    f.close()
    return x, v, phi

def v_circ(m, r):
    return 655.8 * (m/1e14)**0.5 * (r/1.0)**-0.5

def profile_info(x, mp, ok):
    r = np.sqrt(np.sum(x**2, axis=1))
    order = np.argsort(r)
    r_sort = r[order]
    dm = np.ones(len(r_sort))*mp
    dm[~ok] = 0
    m_enc = np.cumsum(dm)
    
    v_rot = v_circ(m_enc, r_sort)
    i_max = np.argmax(v_rot)
    rmax, vmax = r_sort[i_max], v_rot[i_max]
    
    dW = np.zeros(len(m_enc))

    dr = r_sort[1:] - r_sort[:-1]
    dr[dr == 0] = np.finfo(dr.dtype).eps
    r_scaled = (r_sort[1:]*r_sort[:-1]) / dr
    dW[:-1] = v_circ(m_enc[:-1], r_scaled)**2

    vesc_lim = v_circ(m_enc[-1], r_sort[-1]) * np.sqrt(2)

    W = (np.cumsum(dW[::-1])[::-1] + 2*vesc_lim**2)/vmax**2
    out = np.zeros(len(W))
    out[order] = W
    
    return r[i_max], v_rot[i_max], out, order

def interp_vesc(x, vesc, order, bins):
    r = np.sqrt(np.sum(x**2, axis=1))

    r_sorted, vesc_sorted = r[order], vesc[order]
    r_bins = np.linspace(np.log10(r_sorted[0]), np.log10(r_sorted[-1]), bins)
    vesc_bins = vesc_sorted[np.searchsorted(np.log10(r_sorted), r_bins)]

    log_interp_f = interpolate.interp1d(
        r_bins, vesc_bins,
        bounds_error=False,
        fill_value=(vesc_bins[0], vesc_bins[-1])
    )    

    r_max, vesc_max = r_sorted[-1], vesc_sorted[-1]
    def f(rr):
        if rr > r_max:
            return np.sqrt(r_max/rr)*vesc_max
        elif rr <= 0:
            return vesc_sorted[0]
        return log_interp_f(np.log10(rr))
        
    return f

def select_trajectories(n_trj_plot, trj_per_plot, trj_ranks, idx, ranks):
    idx_out = np.zeros(n_trj_plot*trj_per_plot, dtype=int)
    for j in range(trj_per_plot):
        idx_ok = idx[ranks[idx] == trj_ranks[j]]
        choices = random.choice(idx_ok, n_trj_plot, replace=False)
        for i in range(len(choices)):
            idx_out[i*trj_per_plot + j] = choices[i]
    return idx_out

def core_properties(x, v, phi, n_core):
    f_cut = n_core/len(phi)
    phi_cut = np.percentile(phi, 100*f_cut)
    ok = phi <= phi_cut
    return np.mean(x[ok,:], axis=0), np.mean(v[ok,:], axis=0)

def main():
    palette.configure("False")

    params = {'flat': True, 'H0': 70.0, 'Om0': 0.286,
              'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.95}
    
    #base_dir = "../tmp_data/Halo169/"
    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo169/"
    
    #p_snaps = np.arange(127, 227 + 1)
    p_snaps = np.arange(177, 227 + 1)
    #p_snaps = [177, 202, 227]
    h_idx = 4

    particle_energies_fmt = "../plots/energies/particle_energies_%03d.png"
    ang_mom_trj_fmt = "../plots/energies/ang_mom_trj.%02d.png"
    vt_trj_fmt = "../plots/energies/vt_trj.%02d.png"
    vr_trj_fmt = "../plots/energies/vt_trj.%02d.png"
    radii_trj_fmt = "../plots/energies/radii_trj.%02d.png"
    energies_trj_fmt = "../plots/energies/energies_trj.%02d.png"
    
    _, m = lib.read_mergers(base_dir)

    scale = lib.scale_factors()
    z = 1/scale - 1
    age = cosmo.age(z)

    valid_x = None

    nfw_pot = potential.NFWPotential(
        a=1/orbit_model._x_max_nfw(), normalize=True)

    ranks = None

    trj_colors = [pc("r"), pc("o"), pc("g"), pc("b"), pc("p"), pc("k")]
    trj_ranks = [1, 2, 3, 4, 5, -1]
    trj_per_plot = 6
    n_trj_plots = 10
    n_trj = trj_per_plot * n_trj_plots
    
    trj_t = np.zeros(len(p_snaps))
    trj_r = np.zeros((n_trj, len(p_snaps)))
    trj_vt = np.zeros((n_trj, len(p_snaps)))
    trj_vr = np.zeros((n_trj, len(p_snaps)))
    trj_l = np.zeros((n_trj, len(p_snaps)))
    trj_r_apo = np.zeros((n_trj, len(p_snaps)))
    trj_r_peri = np.zeros((n_trj, len(p_snaps)))
    trj_TE = np.zeros((n_trj, len(p_snaps)))
    trj_KE = np.zeros((n_trj, len(p_snaps)))

    i_trj = None
    
    

    for i_snap in range(len(p_snaps)):        
        snap = p_snaps[i_snap]
        print(snap)
        rvir, mvir = m["rvir"][h_idx, snap], m["mvir"][h_idx, snap]
        rvir *= scale[snap]

        vvir0 = 655.8 * ((m["mvir"][0,-1]/1e14)**0.5 *
                         (m["rvir"][0,-1]/1.0)**-0.5)
        vvir = 655.8 * (mvir/1e14)**0.5 * (rvir/1.0)**-0.5

        vmax = m["vmax"][h_idx, snap]
        rmax = orbit_model.rmax_nfw(vmax, mvir, rvir)
        
        x, v, phi = read_part_file(base_dir, snap, h_idx, ["x", "v", "phi"])
        ok = x[:,0] > 0
        
        x_core, v_core = core_properties(x[ok], v[ok], phi[ok], 100)        

        print(x_core, m["x"][h_idx,snap,:])
        print(v_core*np.sqrt(scale[snap]), m["v"][h_idx,snap,:])

        for dim in range(3):
            #v[:,dim] -= v_core[dim] #m["v"][h_idx,snap,dim]
            v[:,dim] *= np.sqrt(scale[snap])
            v[:,dim] -= m["v"][h_idx,snap,dim]
            #x[:,dim] -= x_core[dim] #m["x"][h_idx,snap,dim]
            x[:,dim] -= m["x"][h_idx,snap,dim]
            x[:,dim] *= scale[snap]

        idx = np.arange(len(x))
        x, v, idx, phi = x[ok,:], v[ok,:], idx[ok], phi[ok]
        
        ke_nfw, pe_nfw = orbit_model.energy(rmax, vmax, x, v)
        rmax_orig, vmax_orig, pe_orig, _ = profile_info(x, MP, ok[ok])
        ke = np.sum(v**2, axis=1)/2
        rmax_curr, vmax_curr, pe_curr = rmax_orig, vmax_orig, pe_orig
        for i in range(5):
            is_bound = ke/vmax_curr**2 < pe_curr
            rmax_curr, vmax_curr, pe_curr, order = profile_info(x, MP, is_bound)
        rmax_final, vmax_final, pe_final = rmax_curr, vmax_curr, pe_curr
        vesc_final = np.sqrt(2*pe_final)
        vesc_func = interp_vesc(x, vesc_final, order, 100)
        
        r = np.sqrt(np.sum(x**2, axis=1))
        E = ke/vmax_final**2 - pe_final
 
        if ranks is None:
            ranks = np.ones(len(ok))*-1
            if TAG_METHOD == "energy":
                E_001 = np.percentile(E, 0.1)
                E_005 = np.percentile(E, 0.5)
                E_02 = np.percentile(E, 2)
                E_1 = np.percentile(E, 10)
                E_5 = np.percentile(E, 50)

                ranks[idx[E < E_5]] = 5
                ranks[idx[E < E_1]] = 4
                ranks[idx[E < E_02]] = 3
                ranks[idx[E < E_005]] = 2
                ranks[idx[E < E_001]] = 1
            elif TAG_METHOD == "infall":
                snap_001 = np.percentile(infall_snap, 0.1)
                snap_005 = np.percentile(infall_snap, 0.5)
                snap_02 = np.percentile(infall_snap, 2)
                snap_1 = np.percentile(infall_snap, 10)
                snap_5 = np.percentile(infall_snap, 50)

                ranks[idx[infall_snap <= snap_5]] = 5
                ranks[idx[infall_snap <= snap_1]] = 4
                ranks[idx[infall_snap <= snap_02]] = 3
                ranks[idx[infall_snap <= snap_005]] = 2
                ranks[idx[infall_snap <= snap_001]] = 1

            i_trj = select_trajectories(
                n_trj_plots, trj_per_plot, trj_ranks, idx, ranks
            )
            vmax0, rmax0 = vmax_final, rmax_final
            
        j_trj = np.searchsorted(idx, i_trj)
        trj_t[i_snap] = age[snap]
        trj_TE[:,i_snap] = E[j_trj]*vmax_final**2/vmax0**2
        trj_KE[:,i_snap] = ke[j_trj]*vmax_final**2/vmax0**2
        trj_r[:,i_snap] = r[j_trj]/rmax0
        for j in range(len(j_trj)):
            r_peri, r_apo, ok = orbit_model.peri_apo(
                #rmax_final, vmax_final, x[j_trj[j]], v[j_trj[j]], nfw_pot.vesc
                rmax_final, vmax_final, x[j_trj[j]], v[j_trj[j]], vesc_func
            )
            trj_r_peri[j, i_snap] = r_peri/rmax0
            trj_r_apo[j, i_snap] = r_apo/rmax0
            rr, vR, vT = orbit_model.decompose_velocity(
                x[j_trj[j]]/rmax0, v[j_trj[j]]/vmax0
            )
            trj_l[j,i_snap] = rr*vT
            trj_vt[j,i_snap] = vT
            trj_vr[j,i_snap] = np.abs(vR)

        levels = np.zeros((len(E), 6), dtype=bool)
        for k in range(levels.shape[1]):
            if k == levels.shape[1] - 1:
                levels[:,k] = ranks[idx] == -1
            else:
                levels[:,k] = ranks[idx] == k+1
            print("%.2f" % (np.median(r[levels[:,k]]) * 1e3))
 
        plt.figure(0)
        plt.clf()
        
        plt.plot(r[levels[:,2] | levels[:,3] | levels[:,4] | levels[:,5]],
                 E[levels[:,2] | levels[:,3] | levels[:,4] | levels[:,5]],
                 ".", alpha=0.1, c=pc("k"))
        plt.plot(r[levels[:,2]], E[levels[:,2]],
                 ".", c=pc("g"), label=r"$2\%{\rm\ most\ bound}$")
        plt.plot(r[levels[:,1]], E[levels[:,1]],
                 ".", c=pc("o"), label=r"$0.5\%{\rm\ most\ bound}$")
        plt.plot(r[levels[:,0]], E[levels[:,0]],
                 ".", c=pc("r"), label=r"$0.1\%{\rm\ most\ bound}$")
    
        ylo, yhi = plt.ylim(-6, 6)
        xlo, xhi = plt.xlim(5e-5, 0.5)
        plt.plot([rvir]*2, [ylo, yhi], "--", lw=2, c="k")
        plt.plot([xlo, xhi], [0, 0], "--", lw=2, c=pc("r"))

        plt.legend(loc="upper left", frameon=True, fontsize=16)
        plt.xlabel(r"$r\ (h^{-1}{\rm pMpc})$")
        plt.ylabel(r"$E\ (m_p\,V_{\rm max}^2)$")
        plt.xscale("log")
        
        plt.savefig(particle_energies_fmt % snap)

    j_start = 0
    for i in range(n_trj_plots):
        plt.figure(1)
        plt.clf()
        t = trj_t - trj_t[0]

        for j in range(j_start, j_start + trj_per_plot):
            c = trj_colors[j-j_start]
            plt.plot(t, trj_r[j], c=c)
            ok = trj_r_peri[j] > 0
            plt.plot(t[ok], trj_r_apo[j][ok], "--", c=c, lw=2)

        plt.yscale("log")    
        plt.xlabel(r"$t - t_0\ ({\rm Gyr})$")
        plt.ylabel(r"$r/R_{\rm max,0}$")

        plt.savefig(radii_trj_fmt % i)

        plt.figure(2)
        plt.clf()

        for j in range(j_start, j_start + trj_per_plot):
            c = trj_colors[j-j_start]
            plt.plot(t, trj_TE[j], c=c)

        plt.xlabel(r"$t - t_0\ ({\rm Gyr})$")
        plt.ylabel(r"$E/V^2_{\rm max,0}$")
        plt.savefig(energies_trj_fmt % i)

        plt.figure(3)
        plt.clf()

        for j in range(j_start, j_start + trj_per_plot):
            c = trj_colors[j-j_start]
            plt.plot(t, trj_l[j], c=c)

        plt.yscale("log")

        plt.xlabel(r"$t - t_0\ ({\rm Gyr})$")
        plt.ylabel(r"$L/(m_p V_{\rm max,0}R_{\rm max,0})$")
        plt.savefig(ang_mom_trj_fmt % i)

        plt.figure(4)
        plt.clf()

        for j in range(j_start, j_start + trj_per_plot):
            c = trj_colors[j-j_start]
            plt.plot(t, trj_vr[j], c=c)

        plt.yscale("log")

        plt.xlabel(r"$t - t_0\ ({\rm Gyr})$")
        plt.ylabel(r"$v_R/V_{\rm max,0}$")
        plt.savefig(vr_trj_fmt % i)


        j_start += trj_per_plot
        
        
if __name__ == "__main__": main()
