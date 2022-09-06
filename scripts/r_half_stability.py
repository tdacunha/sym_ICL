import star_tagging
import lib
import numpy as np
import matplotlib.pyplot as plt
import os.path as path
import sys
import struct
import array
import scipy.stats as stats

from colossus.cosmology import cosmology
from colossus.halo import mass_so

import palette
from palette import pc
palette.configure(False)

P_FILE_FMT = "%s/particles/part_%03d.%d"
BASE_DIR_FMT ="/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo%03d/"
N_TRIALS = 50
DM_MP = 2.8e5
h100 = 0.7

cosmo = cosmology.setCosmology("chinchilla",
                               {"flat": True, "H0": 70, "Om0": 0.286,
                                'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.96})


gal_halo = star_tagging.GalaxyHaloModel(
    star_tagging.UniverseMachineMStar(),
    star_tagging.Nadler2020RHalf(),
    star_tagging.PlummerProfile(),
)

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
        f.seek(12*n, 1)
    else:
        v = None

    if "phi" in vars_to_read:
        phi.fromfile(f, n)
        phi = np.array(phi, dtype=float)
    else:
        phi = None

    f.close()
    return x, v, phi

def mass_profile(r_bins, x0, xp, mp):
    """ mass_profile computes the amount of mass in a set of bins defined by
    the bin edges r_bins. The center of the profile is x0, the position of
    each particle is xp, and the particle masses are mp. 
    """
    dxp = np.zeros(xp.shape)
    for dim in range(3):
        dxp[:,dim] = xp[:,dim] - x0[dim]
    rp = np.sqrt(np.sum(dxp**2, axis=1))

    m_sum, _, _ = stats.binned_statistic(rp, mp, "sum", bins=r_bins)

    return m_sum

def mass_to_density(r_bins, mass):
    """ mass_to_density converts 
    """
    vol = 4*np.pi/3 * (r_bins[1:]**3 - r_bins[:-1]**3)
    return mass / vol

def halfmass_radius(xp, x0, mp):
    dxp = np.zeros(xp.shape)
    for dim in range(3):
        dxp[:,dim] = xp[:,dim] - x0[dim]
    rp = np.sqrt(np.sum(dxp**2, axis=1))

    order = np.argsort(rp)
    rp = rp[order]
    m_enc = np.cumsum(mp[order])
    return rp[np.searchsorted(m_enc, m_enc[-1]/2)]


def main():
    # Figure out what halo we're dealing with
    halo_id, sub_idx = int(sys.argv[1]), int(sys.argv[2])
    a_start = float(sys.argv[3])

    base_dir = BASE_DIR_FMT % halo_id

    # Read tree header information
    b_idx, m = lib.read_mergers(base_dir)
    b = lib.read_branches(base_dir)

    # Read the data we need from the tree.
    x, mvir, snap = lib.read_tree(base_dir, ["X", "Mvir", "Snap"])

    # Calculate scale-dependent properties
    scale = lib.scale_factors()
    z = 1/scale - 1
    t_orbit = mass_so.dynamicalTime(z, "vir", "orbit")
    age = cosmo.age(z)

    # calculate merger stats so we can figure out the infall time.
    _, merger_snap, _ = lib.merger_stats(b, m, x, mvir, snap)
    p_sub_idx = lib.pristine_merger_indices(b)

    # Set up buffers for keeping track of mass profiles.
    n_bins = 30
    log_r_bins = np.linspace(np.log10(1e-3), np.log10(100e-3), n_bins+1)
    r_bins = 10**log_r_bins

    star_masses = np.zeros((N_TRIALS, n_bins))
    lmc_masses = np.zeros((N_TRIALS, n_bins))
    gse_masses = np.zeros((N_TRIALS, n_bins))
    dm_masses = np.zeros(n_bins)    

    x_central = m[0]["x"][235]

    # Find halo properties at infall snapshot.
    j = b_idx[sub_idx]
    snap_start = np.searchsorted(scale, a_start)
    snap_end = merger_snap[np.searchsorted(p_sub_idx, j)]


    print(snap_start, snap_end)
    print("%.4f %.4f" % (scale[snap_start], scale[snap_end]))
    print("delta t/t_orbit = %.3f" %
          ((age[snap_end] - age[snap_start])/t_orbit[snap_start]))

    x0 = m["x"][sub_idx][snap_start]
    m0 = m["mvir"][sub_idx][snap_start]
    r0 = m["rvir"][sub_idx][snap_start]
    z0 = 1/scale[snap_start] - 1
    
    n_snap = snap_end - snap_start
    radius_ratios = np.array([0.003, 0.01, 0.03, 0.1])
    n_ratio = len(radius_ratios)

    # read the particle positions at the time of the merger and at the
    # time we want to analyze them
    xp0, _, _ = read_part_file(base_dir, snap_start, sub_idx, ["x"])
    ok = np.zeros((n_ratio, len(xp0)), dtype=bool)
    star_mp = np.zeros((n_ratio, len(xp0)))
    r_half = np.zeros((n_ratio, n_snap))
    
    for j in range(n_ratio):
        ok[j,:] = xp0[:,0] > 0
        star_mp_j = gal_halo.set_m_star(x0, m0, r0, z0, xp0[ok[j,:]],
                                        r_half=r0*radius_ratios[j])
        star_mp[j,ok[j]] = star_mp_j
        r_half[j,0] = halfmass_radius(xp0[ok[j,:]], x0, star_mp_j)


    for snap in range(snap_start+1, snap_end):
        k = snap - snap_start
        if k % 10 == 0: print(snap)
        xk = m["x"][sub_idx][snap]
        xpk, _, _ = read_part_file(base_dir, snap, sub_idx, ["x"])
        for j in range(n_ratio):
            r_half[j,k] = halfmass_radius(xpk[ok[j,:]], xk, star_mp[j,ok[j,:]])


    snap = np.arange(snap_start, snap_end, dtype=int)
    dt = (age[snap] - age[snap_start]) / t_orbit[snap_start]
    
    colors = [pc("r"), pc("o"), pc("b"), pc("p")]
    for j in range(n_ratio):
        scaled_r_half = ((r_half[j,:]*scale[snap]) /
                         (r_half[j,0]*scale[snap_start]))
        plt.plot(dt, scaled_r_half, color=colors[j],
                 label=(r"$R_{\rm half,0}=%g\cdot R_{\rm vir,0}$" %
                        radius_ratios[j]))

    plt.legend(loc="upper left", fontsize=15, frameon=True)
    plt.ylabel(r"$R_{\rm half}/R_{\rm half,0}$")
    plt.xlabel(r"$(t - t_0)/t_{\rm orbit,0}$")
    plt.yscale("log")
    
    ylo, yhi = plt.ylim()
    plt.ylim(ylo, yhi*1.2)
    xlo, xhi =plt.xlim()
    plt.xlim(xlo, xhi)
    plt.plot([xlo, xhi], [1, 1], "--", lw=1.5, c="k")

    plt.savefig("../plots/m_star_profs/m_star_stability_%03d_%d.png" %
                (halo_id, sub_idx))
    
if __name__ == "__main__": main()
