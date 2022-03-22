import star_tagging
import lib
import numpy as np
import matplotlib.pyplot as plt
import os.path as path
import sys
import struct
import array
import scipy.stats as stats

import palette
from palette import pc
palette.configure(False)

P_FILE_FMT = "%s/particles/part_%03d.%d"
BASE_DIR_FMT ="/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo%03d/"
N_TRIALS = 500
DM_MP = 2.8e5
h100 = 0.7
planck_fb = 0.02212 / 0.1206

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

def main():
    # Figure out what halo we're dealing with
    halo_id = int(sys.argv[1])
    base_dir = BASE_DIR_FMT % halo_id

    # Read tree header information
    b_idx, m = lib.read_mergers(base_dir)
    b = lib.read_branches(base_dir)

    # Read the data we need from the tree.
    x, mvir, snap = lib.read_tree(base_dir, ["X", "Mvir", "Snap"])
    scale = lib.scale_factors()

    # calculate merger stats so we can figure out the infall time.
    _, merger_snap, _ = lib.merger_stats(b, m, x, mvir, snap)
    p_sub_idx = lib.pristine_merger_indices(b)

    # Set up buffers for keeping track of mass profiles.
    n_bins = 30
    log_r_bins = np.linspace(np.log10(2e-3), np.log10(300e-3), n_bins+1)
    r_bins = 10**log_r_bins
    star_masses = np.zeros((N_TRIALS, n_bins))
    lmc_masses = np.zeros((N_TRIALS, n_bins))
    gse_masses = np.zeros((N_TRIALS, n_bins))
    dm_masses = np.zeros(n_bins)    

    x_central = m[0]["x"][235]

    # loop over haloes we want to analyze
    for i in range(0, len(m)):
        print(i)
        if i > 0:
            # Find halo properties at infall snapshot.
            j = b_idx[i]
            merger_snap_i = merger_snap[np.searchsorted(p_sub_idx, j)]
            x0 = m["x"][i][merger_snap_i]
            m0 = m["mvir"][i][merger_snap_i]
            r0 = m["rvir"][i][merger_snap_i]
            z0 = 1/scale[merger_snap_i] - 1
            
            # read the particle positions at the time of the merger and at the
            # time we want to analyze them
            xp_merger, _, _ = read_part_file(base_dir, merger_snap_i, i, ["x"])
            xp, _, _ = read_part_file(base_dir, 235, i, ["x"])

            # for each trial, generate a new set of star particles according to
            # a randomly drawn r_half and m_star.
            for k in range(N_TRIALS):
                star_mp = gal_halo.set_m_star(x0, m0, r0, z0, xp_merger)
                star_mass = mass_profile(r_bins, x_central, xp, star_mp)
                star_masses[k,:] += star_mass
                if i == 1:
                    lmc_masses[k,:] += star_mass
                elif i == 2:
                    gse_masses[k,:] += star_mass
        else:
            # We don't track star particles for the host, so just read its 
            # dark matter particles
            xp, _, _ = read_part_file(base_dir, 235, i, ["x"])
            
        dm_mp = np.ones(len(xp))*DM_MP
        dm_masses += mass_profile(r_bins, x_central, xp, dm_mp)

    # Nothing interesting happens form here on: just a bunch of conversion 
    # from mass profiles to density profiles, unit conversions, percentile
    # calculations, plotting functions etc.
    dm_prof = mass_to_density(r_bins, dm_masses)
    star_prof = np.zeros(star_masses.shape)
    lmc_prof = np.zeros(star_masses.shape)
    gse_prof = np.zeros(star_masses.shape)
    for k in range(N_TRIALS):
        star_prof[k,:] = mass_to_density(r_bins, star_masses[k,:])
        lmc_prof[k,:] = mass_to_density(r_bins, lmc_masses[k,:])
        gse_prof[k,:] = mass_to_density(r_bins, gse_masses[k,:])
   
    r = 10**((log_r_bins[1:] + log_r_bins[:-1]) / 2)
    r *= 1e3/h100
    star_prof *= 1e-6 * h100**2
    lmc_prof *= 1e-6 * h100**2
    gse_prof *= 1e-6 * h100**2
    dm_prof *= 1e-6 * h100**2

    f_star = np.median(np.sum(star_masses, axis=1)) / np.sum(dm_masses)

    lo_star_prof = np.percentile(star_prof, 50-68/2, axis=0)
    med_star_prof = np.percentile(star_prof, 50, axis=0)
    hi_star_prof = np.percentile(star_prof, 50+68/2, axis=0)

    lo_gse_prof = np.percentile(gse_prof, 50-68/2, axis=0)
    med_gse_prof = np.percentile(gse_prof, 50, axis=0)
    hi_gse_prof = np.percentile(gse_prof, 50+68/2, axis=0)

    lo_lmc_prof = np.percentile(lmc_prof, 50-68/2, axis=0)
    med_lmc_prof = np.percentile(lmc_prof, 50, axis=0)
    hi_lmc_prof = np.percentile(lmc_prof, 50+68/2, axis=0)

    plt.yscale("log")
    plt.xscale("log")

    plt.plot(r, dm_prof*r*r*f_star, "--", lw=2, c=pc("k"),
             label=r"$\rho_{\rm DM}\cdot(M_\star/M_{\rm DM})$")
    plt.plot(r, med_star_prof*r*r, c=pc("r"), label=r"$\rho_\star$")
    plt.fill_between(r, lo_star_prof*r*r, hi_star_prof*r*r,
                     alpha=0.2, color=pc("r"))

    plt.plot(r, med_lmc_prof*r*r, c=pc("b"), label=r"$\rho_{\star,1}$")
    plt.fill_between(r, lo_lmc_prof*r*r, hi_lmc_prof*r*r,
                     alpha=0.2, color=pc("b"))

    plt.plot(r, med_gse_prof*r*r, c=pc("o"), label=r"$\rho_{\star,2}$")
    plt.fill_between(r, lo_gse_prof*r*r, hi_gse_prof*r*r,
                     alpha=0.2, color=pc("o"))

    plt.xlabel(r"$r\ ({\rm kpc})$")
    plt.ylabel(r"$r^2\cdot\rho\ (M_\odot/{\rm kpc})$")
    plt.legend(loc="lower right", frameon=True)
    plt.savefig("../plots/m_star_profs/m_star_prof_%03d.png" % halo_id)

if __name__ == "__main__": main()
