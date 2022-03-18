import lib
import orbit_model
import struct
import array
import numpy as np
import galpy.potential as potential
import scipy.stats as stats
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
    
    return r[i_max], v_rot[i_max], out

def select_trajectories(n_trj_plot, trj_per_plot, trj_ranks, idx, ranks):
    idx_out = np.zeros(n_trj_plot*trj_per_plot, dtype=int)
    for j in range(trj_per_plot):
        idx_ok = idx[ranks[idx] == trj_ranks[j]]
        choices = random.choice(idx_ok, n_trj_plot, replace=False)
        for i in range(len(choices)):
            idx_out[i*trj_per_plot + j] = choices[i]
    return idx_out

def main():
    palette.configure("False")

    params = {'flat': True, 'H0': 70.0, 'Om0': 0.286,
              'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.95}

    base_dir = "../tmp_data/Halo169/"
    p_snaps = [127, 177, 227]
    h_idx = 4
    p_files = [
        "../tmp_data/Halo169/particles/part_127.4",
        "../tmp_data/Halo169/particles/part_177.4",
        "../tmp_data/Halo169/particles/part_227.4"
    ]

    particle_energies_fmt = "../plots/energies/particle_energies_%03d.png"
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
        
        for dim in range(3):
            v[:,dim] *= scale[snap] 
            v[:,dim] -= m["v"][h_idx,snap,dim]
            x[:,dim] -= m["x"][h_idx,snap,dim]
            x[:,dim] *= scale[snap]

        idx = np.arange(len(x))
        x, v, idx = x[ok,:], v[ok,:], idx[ok]
        
        ke_nfw, pe_nfw = orbit_model.energy(rmax, vmax, x, v)
        rmax_orig, vmax_orig, pe_orig = profile_info(x, MP, ok[ok])
        ke = np.sum(v**2, axis=1)/2
        rmax_curr, vmax_curr, pe_curr = rmax_orig, vmax_orig, pe_orig
        for i in range(5):
            is_bound = ke/vmax_curr**2 < pe_curr
            rmax_curr, vmax_curr, pe_curr = profile_info(x, MP, is_bound)
        rmax_final, vmax_final, pe_final = rmax_curr, vmax_curr, pe_curr
        
        r = np.sqrt(np.sum(x**2, axis=1))
        E = ke/vmax_final**2 - pe_final

        if ranks is None:
            ranks = np.ones(len(ok))*-1
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
                rmax_final, vmax_final, x[j_trj[j]], v[j_trj[j]], nfw_pot.vesc
            )
            trj_r_peri[j, i_snap] = r_peri/rmax0
            trj_r_apo[j, i_snap] = r_apo/rmax0
        
        levels = np.zeros((len(E), 6), dtype=bool)
        for k in range(levels.shape[1]):
            if k == levels.shape[1] - 1:
                levels[:,k] = ranks[idx] == -1
            else:
                levels[:,k] = ranks[idx] == k+1
            print("%.2f" % (np.median(r[levels[:,k]]) * 1e3))
        print()
            
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
        t = trj_t - trj_t[0]

        for j in range(j_start, j_start + trj_per_plot):
            c = trj_colors[j-j_start]
            plt.plot(t, trj_TE[j], c=c)

        plt.xlabel(r"$t - t_0\ ({\rm Gyr})$")
        plt.ylabel(r"$E/V^2_{\rm max,0}$")
        plt.savefig(energies_trj_fmt % i)

        j_start += trj_per_plot
        
        
if __name__ == "__main__": main()
