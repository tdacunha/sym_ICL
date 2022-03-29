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
import lib
from palette import pc
from colossus.cosmology import cosmology
from colossus.halo import mass_so
import numpy.linalg as linalg

params = { 'flat': True, 'H0': 70.0, 'Om0': 0.286,
           'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.95 }
cosmo = cosmology.setCosmology('myCosmo', params)

QUANTILES = 10**np.linspace(-4, 0, 12)

P_FILE_FMT = "%s/particles/part_%03d.%d"
MP = 2.8e5

MIN_TRACERS = 10

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

def clean_particles(x, v, phi, h, scale, ok=None):
    """ clean_particles centers particles and covnerts them to phsycial
    units and removes unused particles. x, v, and phi come from the particle
    file. v and phi may be None. h is a halo (a single element of the mergers
    array) and scale is the scale factor. Returns cleaned x, v, and phi (None if
    those arrays were intially None), and the indices of these particles into
    the original array.
    """
    if ok is None:
        ok = x[:,0] > 0
    else:
        ok = ok & (x[:,0] > 0)
    
    for dim in range(3):
        if v is not None:
            v[:,dim] *= np.sqrt(scale)
            v[:,dim] -= h["v"][dim]
        x[:,dim] -= h["x"][dim]
        x[:,dim] *= scale

    idx = np.arange(len(x))
    x, idx = x[ok], idx[ok]
    if v is not None:
        v = v[ok]
    if phi is not None:
        phi = phi[ok]

    return x, v, phi, idx

def rank_by_quantile(x, quantiles, idx, n_max):
    cutoffs = np.zeros(len(quantiles)+1)
    cutoffs[0] = np.min(x)
    cutoffs[1:] = np.quantile(x, quantiles)
    ranks = np.ones(n_max, dtype=int)*-1
    for i in range(len(quantiles)):
        ok = (ranks[idx] == -1) & (x < cutoffs[i+1])
        ranks[idx[ok]] = i
        
    return ranks

def rank_by_radius(base_dir, snap, h_idx, ok=None):
    _, m = lib.read_mergers(base_dir)
    scale = lib.scale_factors()
    
    x, _, _ = lib.read_part_file(base_dir, snap, h_idx, ["x"])
    n_max, n_ranks = len(x), len(QUANTILES)
    x, _, _, idx = clean_particles(x, None, None,
                                   m[h_idx,snap], scale[snap], ok=ok)
    r = np.sqrt(np.sum(x**2, axis=1))
    
    ranks = rank_by_quantile(r, QUANTILES, idx, n_max)
    return ranks

def rank_by_nfw_energy(base_dir, snap, h_idx):
    _, m = lib.read_mergers(base_dir)
    scale = lib.scale_factors()

    x, v, _ = lib.read_part_file(base_dir, snap, h_idx, ["x", "v"])
    n_max, n_ranks = len(x), len(QUANTILES)
    x, v, _, idx = clean_particles(x, v, None, m[h_idx,snap], scale[snap])

    mvir, rvir = m[h_idx,snap]["mvir"], scale[snap]*m[h_idx,snap]["rvir"]
    vmax = m[h_idx,snap]["vmax"]
    rmax = orbit_model.rmax_nfw(vmax, mvir, rvir)
    ke, pe = orbit_model.energy(rmax, vmax, x, v)
    E = ke - pe

    ranks = rank_by_quantile(E, QUANTILES, idx, n_max)
    return ranks


def rank_by_radial_energy(base_dir, snap, h_idx):
    _, m = lib.read_mergers(base_dir)
    scale = lib.scale_factors()
    
    x, v, _ = lib.read_part_file(base_dir, snap, h_idx, ["x", "v"])
    n_max, n_ranks = len(x), len(QUANTILES)
    x, v, _, idx = clean_particles(x, v, None, m[h_idx,snap], scale[snap])

    ke = np.sum(v**2, axis=1)/2

    is_bound = np.ones(len(x), dtype=bool)
    rmax, vmax, pe, _ = profile_info(x, MP, is_bound)
    for i in range(5):
        is_bound = ke/vmax**2 < pe
        rmax, vmax, pe, order = profile_info(x, MP, is_bound)
    E = ke/vmax**2 - pe

    ranks = rank_by_quantile(E, QUANTILES, idx, n_max)
    return ranks

    
def rank_by_tree_energy(base_dir, snap, h_idx):
    _, m = lib.read_mergers(base_dir)
    scale = lib.scale_factors()
    
    x, v, phi = lib.read_part_file(base_dir, snap, h_idx, ["x", "v", "phi"])
    n_max, n_ranks = len(x), len(QUANTILES)
    x, v, phi, idx = clean_particles(x, v, phi, m[h_idx,snap], scale[snap])

    ke = np.sum(v**2, axis=1)/2
    is_bound = np.ones(len(x), dtype=bool)
    rmax, vmax, pe_rad, _ = profile_info(x, MP, is_bound)
    ok = phi != -1
    pe_rad, phi, x, idx, ke = pe_rad[ok], phi[ok], x[ok], idx[ok], ke[ok]

    phi *= np.mean(pe_rad/phi)
    E = ke/vmax**2 - phi

    ranks = rank_by_quantile(E, QUANTILES, idx, n_max)
    return ranks

def rank_by_infall(base_dir, snap, h_idx):    
    _, infall = lib.read_tags(base_dir, h_idx)
    n_max, n_ranks = len(infall), len(QUANTILES)
    idx = np.arange(len(infall), dtype=int)
    infall, idx = infall[infall <= snap], idx[infall<= snap]

    ranks = rank_by_quantile(infall, QUANTILES, idx, n_max)
    return ranks
    

def rank_by_average_radius(base_dir, snaps, h_idx):
    _, m = lib.read_mergers(base_dir)
    scale = lib.scale_factors()
    
    for i_snap, snap in enumerate(snaps):
        x, _, _ = lib.read_part_file(base_dir, snap, h_idx, ["x"])
        if i_snap == 0:
            n_max, n_ranks = len(x), len(QUANTILES)

            r_sum = np.zeros(len(x))
            n_snap = np.zeros(len(x), dtype=int)

        x, _, _, idx = clean_particles(
            x, None, None, m[h_idx,snap], scale[snap])

        r = np.sqrt(np.sum(x**2, axis=1))
        r_sum[idx] += r
        n_snap[idx] += 1
    
    r_avg = np.ones(len(r_sum))*np.inf
    r_avg[n_snap > 0] = r_sum[n_snap > 0] / n_snap[n_snap > 0]

    ranks = rank_by_quantile(r_avg[idx], QUANTILES, idx, n_max)
    return ranks

def rank_by_average_nfw_energy(base_dir, snaps, h_idx):
    _, m = lib.read_mergers(base_dir)
    scale = lib.scale_factors()

    for i_snap, snap in enumerate(snaps):
        x, v, _ = lib.read_part_file(base_dir, snap, h_idx, ["x", "v"])

        if i_snap == 0:
            n_max, n_ranks = len(x), len(QUANTILES)
            E_sum = np.zeros(len(x))
            n_snap = np.zeros(len(x), dtype=int)
            
            mvir = m[h_idx,snap]["mvir"]
            rvir = scale[snap]*m[h_idx,snap]["rvir"]
            vmax = m[h_idx,snap]["vmax"]
            rmax = orbit_model.rmax_nfw(vmax, mvir, rvir)

        x, v, _, idx = clean_particles(
            x, v, None, m[h_idx,snap], scale[snap])

        ke, pe = orbit_model.energy(rmax, vmax, x, v)
        E = ke - pe

        E_sum[idx] += E
        n_snap[idx] += 1
    
    E_avg = np.ones(len(E_sum))*np.inf
    E_avg[n_snap > 0] = E_sum[n_snap > 0] / n_snap[n_snap > 0]
    
    ranks = rank_by_quantile(E_avg[idx], QUANTILES, idx, n_max)
    return ranks


def rank_radii(m, base_dir, snap, h_idx, ranks, q):
    scale = lib.scale_factors()
    x, _, _ = lib.read_part_file(base_dir, snap, h_idx, ["x"])
    x, _, _, idx = clean_particles(x, None, None, m[h_idx,snap], scale[snap])

    r = np.sqrt(np.sum(x**2, axis=1))
    rvir = scale[snap]*m[h_idx,snap]["rvir"]

    rr = np.zeros((len(ranks), np.max(ranks[0])+1, len(q)))

    for i in range(rr.shape[0]):
        for j in range(rr.shape[1]):
            in_rank = np.where(ranks[i][idx] == j)[0]

            if len(in_rank) < MIN_TRACERS:
                rr[i,j,:] = -1
            else:
                for k in range(rr.shape[2]):
                    rr[i,j,k] = np.quantile(r[in_rank], q[k])

    return rr, x, idx

def decompose_velocity(x, v):
    r = np.sqrt(np.sum(x**2, axis=1))
    x_hat = np.zeros(x.shape)
    for dim in range(3): x_hat[:,dim] = x[:,dim]/r

    vr = np.sum(x_hat*v, axis=1)
    vTot = np.sqrt(np.sum(v**2, axis=1))
    vT = np.sqrt(vTot**2 - vr**2)

    return r, vr, vT

def mean_properties_by_rank(x, v, idx, rank, m, snap, h_idx):
    scale = lib.scale_factors()

    r, vr, vT = decompose_velocity(x, v)
    c = vT/np.sqrt(vr**2 + vT**2)

    mvir, rvir = m[h_idx,snap]["mvir"], scale[snap]*m[h_idx,snap]["rvir"]
    vmax = m[h_idx,snap]["vmax"]
    rmax = orbit_model.rmax_nfw(vmax, mvir, rvir)
    ke, pe = orbit_model.energy(rmax, vmax, x, v)

    E = ke - pe

    n_max = np.max(rank)

    mean_c = np.zeros(n_max+1)
    mean_E = np.zeros(n_max+1)
    mean_r = np.zeros(n_max+1)

    for i in range(len(mean_c)):
        ok = np.where(rank[idx] == i)[0]
        mean_c[i] = np.median(c[ok])
        mean_E[i] = np.median(E[ok])
        mean_r[i] = np.median(r[ok])

    return mean_c, mean_E, mean_r

def get_circularity(base_dir, snap, h_idx, rank):
    _, m = lib.read_mergers(base_dir)
    scale = lib.scale_factors()
    
    x, v, _ = lib.read_part_file(base_dir, snap, h_idx, ["x", "v"])
    n_max, n_ranks = len(x), len(QUANTILES)
    x, v, _, idx = clean_particles(x, v, None, m[h_idx,snap], scale[snap])

    return mean_properties_by_rank(x, v, idx, rank, m, snap, h_idx)

def ellipse(points, w=None):
    x, y, z = points[:,0], points[:,1], points[:,2]
    if w is None: w = np.ones(len(x))

    M11, M12, M13 = np.sum(x*x*w), np.sum(x*y*w), np.sum(x*z*w)
    M21, M22, M23 = M12,           np.sum(y*y*w), np.sum(y*z*w)
    M31, M32, M33 = M13,           M23,           np.sum(z*z*w)
    M = np.array([
        [M11,M12,M13],
        [M21,M22,M23],
        [M31,M32,M33]
    ])

    lamb, vec = linalg.eig(M)
    axes = np.sqrt(5*lamb/np.sum(w))
    return axes, vec.T

def axis_ratio_by_rank(x, idx, rank, weighting=None):
    if weighting == "1/r^2":
        w = 1/(x[:,0]**2 + x[:,1]**2 + x[:,2]**2)
    elif weighting == "1/r":
        w = 1/np.sqrt(x[:,0]**2 + x[:,1]**2 + x[:,2]**2)
    else:
        w = np.ones(len(x))
    
    max_rank = np.max(rank)
    c_to_a, b_to_a = np.zeros(max_rank+1), np.zeros(max_rank+1)
    
    for i in range(len(c_to_a)):
        ok = np.where(rank[idx] == i)[0]
        if len(ok) < MIN_TRACERS:
            c_to_a[i], b_to_a[i] = -1, -1
            continue

        axes, _ = ellipse(x[ok], w[ok])
        c, a = np.min(axes), np.max(axes)
        b = np.sum(axes) - c - a

        c_to_a[i], b_to_a[i] = c/a, b/a

    return c_to_a, b_to_a

def main():
    palette.configure(False) 
    host_id = 169
    base_dir = ("/oak/stanford/orgs/kipac/users/phil1/" + 
                "simulations/MWest/Halo%03d/" % host_id)

    h_idx, scale_start, scale_infall = {
        169: (4, 0.25, 0.9),
        222: (6, 0.3, 0.95),
        523: (10, 0.35, 0.9),
        719: (1, 0.4, 0.85)
    }[host_id]
    scale = lib.scale_factors()
    snap_start, snap_infall = np.searchsorted(
        scale, [scale_start, scale_infall])

    pre_snaps = np.arange(snap_start - 10, snap_start)

    _, m = lib.read_mergers(base_dir)

    scale = lib.scale_factors()
    z = 1/scale - 1
    age = cosmo.age(z)
    t_orbit = mass_so.dynamicalTime(z, "vir", "orbit")

    dt = (age - age[snap_start])
    dt_orbit = dt/t_orbit[snap_start]

    comp_snaps = np.searchsorted(dt_orbit, [0, 0.5, 1, 2, 4])
    comp_snaps = comp_snaps[comp_snaps < snap_infall]
    print("comp_snaps:", comp_snaps)
    comp_snaps = np.arange(snap_start, snap_infall+1)

    snap0 = comp_snaps[0]
    ranks = [
        rank_by_radius(base_dir, snap0, h_idx),
        rank_by_radial_energy(base_dir, snap0, h_idx),
        rank_by_nfw_energy(base_dir, snap0, h_idx),
        rank_by_infall(base_dir, snap0, h_idx),
        rank_by_average_radius(base_dir, comp_snaps, h_idx),
    ]

    colors = [pc("r"), pc("o"), pc("g"),  pc("b"), pc("p"), pc("n")]
    names = [r"$r$",
             r"${\rm KE} + {\rm PE}(r)$",
             r"${\rm KE} + {\rm PE}_{\rm NFW}$",
             r"$t_{\rm infall}$",
             r"${\rm optimal}$"]

    rr0, x0, idx0 = rank_radii(
        m, base_dir, snap_start, h_idx, ranks, [0.1, 0.5, 0.9])
    rvir0 = m[h_idx,snap_start]["rvir"]*scale[snap_start]
    
    c_to_a_weighting = "1/r^2"
    c_to_a_0 = np.zeros(rr0.shape[:2])
    for i in range(len(ranks)):
        c_to_a_0[i,:], _ = axis_ratio_by_rank(
            x0, idx0, ranks[i], c_to_a_weighting)

    for k in range(len(comp_snaps)):
        plt.figure(0)
        plt.clf()
        
        snap = comp_snaps[k]
        dt_string = ((r"$\Delta t = %.1f\ {\rm Gyr} = " +
                      r"%.1f\cdot t_{\rm orbit,vir}$") %
                     (dt[snap], dt_orbit[snap]))

        rr, x, idx = rank_radii(
            m, base_dir, snap, h_idx, ranks, [0.1, 0.5, 0.9])

        for i in range(rr.shape[0]):
            ok = rr[i,:,1] >= 0
            plt.plot(QUANTILES[ok], rr[i,ok,1]/rvir0, label=names[i],
                     c=colors[i], lw=3)
            plt.plot(QUANTILES[ok], rr0[i,ok,1]/rvir0, "--",
                     c=colors[i], lw=2) 


        rank_ref = rank_by_radius(base_dir, snap, h_idx, ok=ranks[0]>=0)
        rr_ref, _, _ = rank_radii(m, base_dir, snap, h_idx,
                            [rank_ref], [0.1, 0.5, 0.9])
        ok = rr_ref[0,:,1] >= 0
        plt.plot(QUANTILES[ok], rr_ref[0,ok,1]/rvir0,
                 "--", c="k", lw=3, label=r"$r_{\rm current}$") 

        plt.xscale("log")
        plt.yscale("log")
        plt.title(dt_string)

        lo, hi = plt.xlim()
        plt.xlim(lo, hi)
        ylo, yhi = plt.ylim(2e-3, 2)
        plt.ylim(ylo, yhi)
    
        e_rvir = 4*170e-6/rvir0 * scale[snap]

        plt.fill_between([lo, hi], [ylo]*2, [e_rvir]*2, alpha=0.2, color="k")

        plt.legend(loc="upper left", fontsize=16)
        plt.xlabel(r"$p$")
        plt.ylabel(r"$R_{\rm 1/2}/R_{\rm vir,0}$")
        plt.title(dt_string)

        plt.savefig(
            ("../plots/m_star_profs/profile_evolution/fine/" + 
             "quantile_evolution.%d.%d.%d.png") %
            (host_id, h_idx, snap)
        )

        plt.figure(1)
        plt.clf()

        for i in range(3):
            c_to_a, _ = axis_ratio_by_rank(x, idx, ranks[i], c_to_a_weighting)
            ok = c_to_a >= 0

            plt.plot(QUANTILES[ok], c_to_a[ok], label=names[i],
                     c=colors[i], lw=3)
            plt.plot(QUANTILES[ok], c_to_a_0[i,ok], "--",
                     c=colors[i], lw=2) 

        plt.legend(loc="upper left", fontsize=16)
        plt.xlabel(r"$p$")
        plt.ylabel(r"$c/a$")
        plt.xscale("log")
        plt.title(dt_string)
        plt.ylim(0, 1.2)

        plt.savefig(
            ("../plots/m_star_profs/profile_evolution/fine/" + 
             "c_to_a_evolution.%d.%d.%d.png") %
            (host_id, h_idx, snap)
        )

if __name__ == "__main__": main()
