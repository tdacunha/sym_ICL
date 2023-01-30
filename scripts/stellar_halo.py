import numpy as np
import matplotlib.pyplot as plt
import symlib

FE_H_LOW = -1.5
FE_H_HIGH = -0.5

try:
    import palette
    palette.configure(False)
except:
    pass

def final_radii(sim_dir, part_info, target_subs=None):
    # Standard boiler-plate code that reads in various parameters and arrays
    param = symlib.simulation_parameters(sim_dir)
    scale = symlib.scale_factors(sim_dir)
    # Subhalos
    h, _ = symlib.read_subhalos(sim_dir)
    # Subhalos in comoving units
    h_cmov, _ = symlib.read_subhalos(sim_dir, comoving=True)
    # Tracked cores
    c = symlib.read_cores(sim_dir)

    if target_subs is None:
        target_subs = np.arange(1, len(h))

    last_snap = param["n_snap"] - 1

    # The cut that determines whether a tracked core is still intact
    core_r = np.sqrt(np.sum(c["x"]**2, axis=2))
    mp = param["mp"]/param["h100"]
    # core tracking must have started, core needs 32 particles, and it can't
    # overlapping qith the host's center.
    core_intact = c["ok"] & (c["m_bound"] > 32*mp) & (core_r > c["r50_bound"])

    rf, in_halo = [None]*len(h), [None]*len(h)

    for i in target_subs:
        # Only work with particles where this flag is true.
        ok = symlib.read_particles(part_info, sim_dir, last_snap,
                                   "valid", owner=i)
        x = symlib.read_particles(part_info, sim_dir, last_snap, "x", owner=i)

        # Correct the units. There's an analogous function for velocities.
        x_i = symlib.set_units_x(x, h_cmov[0,-1], scale[-1], param)

        r_host = np.sqrt(np.sum(x_i**2, axis=1))
        rf[i] = np.ones(len(ok))*-1
        rf[i][ok] = r_host[ok]

        r_sub = np.sqrt(np.sum((x_i - c[i,-1]["x"])**2, axis=1))
        in_halo[i] = np.zeros(len(ok), dtype=bool)

        if core_intact[i,-1]:
            in_halo[i][ok] = r_sub[ok] > c[i,-1]["r_tidal"]
        else:
            in_halo[i][ok] = True

    return rf, in_halo

def density_profile(sim_dir, r_bins, r, mp_star, in_halo,
                    include_sub=None, include_part=None,
                    target_subs=None):
    """ This function just generates a density profile. Nothing too
    interesting.
    """
    h, _ = symlib.read_subhalos(sim_dir)

    if target_subs is None:
        target_subs = np.arange(1, len(h))

    r_host = h["rvir"][0,-1]
    r_bins = r_bins*r_host # Locally, we'll convert out fo normalized units
    V_bins = (r_bins[1:]**3 - r_bins[:-1]**3)*4*np.pi/3

    n_bins = len(r_bins) - 1
    mass_hist = np.zeros(n_bins)

    for i in target_subs:
        if include_sub is not None and not include_sub[i]: continue

        ok = in_halo[i] & (mp_star[i] > 0)
        if include_part is not None:
            ok = ok & include_part[i]

        m, _ = np.histogram(r[i][ok], bins=r_bins, weights=mp_star[i][ok])
        mass_hist += m

    return mass_hist / V_bins

def weighted_mean_profile(sim_dir, r_bins, r, mp_star, x,
                          in_halo, include_sub,
                          target_subs=None):
    """ This function just generates a density profile. Nothing too
    interesting.
    """
    h, _ = symlib.read_subhalos(sim_dir)

    if target_subs is None:
        target_subs = np.arange(1, len(h))

    r_host = h["rvir"][0,-1]
    r_bins = r_bins*r_host # Locally, we'll convert out fo normalized units

    n_bins = len(r_bins) - 1
    mass_hist = np.zeros(n_bins)
    x_hist = np.zeros(n_bins)

    for i in target_subs:
        if not include_sub[i]: continue

        ok = in_halo[i] & (mp_star[i] > 0)
        idx = np.arange(np.sum(ok))
        m, _ = np.histogram(r[i][ok], bins=r_bins, weights=mp_star[i][ok])
        x_i, _ = np.histogram(r[i][ok], bins=r_bins, weights=mp_star[i][ok]*x[i][ok])
        mass_hist += m
        x_hist += x_i

    ok = mass_hist > 0
    ratio = np.ones(n_bins)*-np.inf
    ratio[ok] = x_hist[ok]/mass_hist[ok]

    return ratio


def main():
    base_dir ="/sdf/scratch/phil1/data/" #"/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/\
    suite = "SymphonyMilkyWay"
    plot_dir = "/sdf/home/t/tdacunha/plots/stellar_halo/"

    param = symlib.simulation_parameters(suite)
    n_hosts = 1#symlib.n_hosts(suite)

    # Set up density profiles
    n_bins = 50
    r_bins = np.logspace(-2, 0, n_bins + 1) # In units of Rvir(z=0)
    rho_high_mass = []
    rho_low_mass = []
    rho_high_Fe_H = []
    rho_mid_Fe_H = []
    rho_low_Fe_H = []
    mean_Fe_H_low_mass = []
    mean_Fe_H_high_mass = []

    # Our galaxy-halo model
    gal_halo = symlib.GalaxyHaloModel(
        symlib.UniverseMachineMStarFit(),
        symlib.Jiang2019RHalf(),
        symlib.PlummerProfile(),
        symlib.Kirby2013Metallicity(),
        no_scatter=False
    )

    # Arrays needed to store plot data
    m_star_halo = []
    m_star_gal = []
    mpeak_sub = []
    n_hosts_used = 0
    for i in [n_hosts]:#range(n_hosts):
        print("Host %2d/%d" % (i, n_hosts))

        # This function lets you loop over all the subhalo directories without
        # needing to know them by name.
        sim_dir = symlib.get_host_directory(base_dir, suite, "Halo023")

        # Shared information on ownership and offsets needed to
        # decompress files.
        part_info = symlib.ParticleInfo(sim_dir)

        h, hist = symlib.read_subhalos(sim_dir)

        # Set to whatever halos you want to analyze. I'd suggest starting with
        # a couple small ones for debugging and then changing it back to
        # everything later. Don't include halo 0.
        target_subs = np.arange(1, len(h))

        target_subs_ok = np.zeros(len(target_subs), dtype=bool)
        for i in range(len(target_subs)):
            target_subs_ok[i] = symlib.is_real_confirmed(
                part_info, h, target_subs[i])
        target_subs = target_subs[target_subs_ok]
        print(target_subs)
        # mp_star is a list of arrays giving the stellar mass of each particle.
        mp_star, _, m_star, r_half, Fe_H = symlib.tag_stars(
            sim_dir, gal_halo, target_subs=target_subs)

        # Example analysis function. rf and in_halo are both lists of arrays
        # rf gives the final radius of particles at z=0 relative to the host
        # center and in_halo is true if the particle is in the host's halo
        # and false otherwise. rf will be -1 for invalid particles.
        rf, in_halo = final_radii(sim_dir, part_info,
                                  target_subs=target_subs)

        is_high_mass = hist["merger_ratio"] > 0.15

        # Keep track of contribution to total halo/satellite mass.
        is_low_Fe_H = [None]*len(h)
        is_mid_Fe_H = [None]*len(h)
        is_high_Fe_H = [None]*len(h)
        for j in target_subs:
            m_star_halo.append(np.sum(mp_star[j][in_halo[j]]))
            m_star_gal.append(np.sum(mp_star[j][~in_halo[j]]))
            mpeak_sub.append(hist["mpeak"][j])
            is_high_Fe_H[j] = Fe_H[j] > FE_H_HIGH
            is_low_Fe_H[j] = Fe_H[j] < FE_H_LOW
            is_mid_Fe_H[j] = (~is_high_Fe_H[j]) & (~is_low_Fe_H[j])

        # Compute density profiles, add to running averages.
        rho_high_mass.append(density_profile(
            sim_dir, r_bins, rf, mp_star, in_halo,
            include_sub=is_high_mass,
            target_subs=target_subs))
        rho_low_mass.append(density_profile(
            sim_dir, r_bins, rf, mp_star, in_halo,
            include_sub=~is_high_mass,
            target_subs=target_subs))

        rho_low_Fe_H.append(density_profile(
            sim_dir, r_bins, rf, mp_star, in_halo,
            include_part=is_low_Fe_H,
            target_subs=target_subs))
        rho_mid_Fe_H.append(density_profile(
            sim_dir, r_bins, rf, mp_star, in_halo,
            include_part=is_mid_Fe_H,
            target_subs=target_subs))
        rho_high_Fe_H.append(density_profile(
            sim_dir, r_bins, rf, mp_star, in_halo,
            include_part=is_high_Fe_H,
            target_subs=target_subs))

        mean_Fe_H_high_mass.append(weighted_mean_profile(
            sim_dir, r_bins, rf, mp_star, Fe_H, in_halo, is_high_mass,
            target_subs=target_subs))
        mean_Fe_H_low_mass.append(weighted_mean_profile(
            sim_dir, r_bins, rf, mp_star, Fe_H, in_halo, ~is_high_mass,
            target_subs=target_subs))

        n_hosts_used += 1

    m_star_halo = np.array(m_star_halo)
    m_star_gal = np.array(m_star_gal)
    mpeak_sub = np.array(mpeak_sub)

    # Plotting nonsense from here on
    fig, ax = plt.subplots()
    r_mid = np.sqrt(r_bins[1:]*r_bins[:-1])

    med_rho_high_mass = np.median(rho_high_mass, axis=0)
    med_rho_low_mass = np.median(rho_low_mass, axis=0)

    for i in range(len(rho_high_mass)):
        ax.plot(r_mid, rho_high_mass[i], lw=1, alpha=0.5, c="tab:red")
        ax.plot(r_mid, rho_low_mass[i], lw=1, alpha=0.5, c="tab:blue")

    ax.plot(r_mid, med_rho_high_mass, c="tab:red",
             label=r"$M_{\rm sub,infall}/M_{\rm host,infall} > 0.15$")
    ax.plot(r_mid, med_rho_low_mass, c="tab:blue",
             label=r"$M_{\rm sub,infall}/M_{\rm host,infall} \leq 0.15$")
    ax.plot(r_mid, med_rho_high_mass + med_rho_low_mass, c="k")
    ax.set_xlabel(r"$r/R_{\rm vir}$")
    ax.set_ylabel(r"$\rho\ (M_\odot\,{\rm kpc}^{-3})$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, None)
    ax.legend(loc="upper right", fontsize=17)

    fig.savefig(plot_dir+"average_density_%s.png" % suite)

    fig, ax = plt.subplots()
    r_mid = np.sqrt(r_bins[1:]*r_bins[:-1])

    med_rho_high_Fe_H = np.median(rho_high_Fe_H, axis=0)
    med_rho_mid_Fe_H = np.median(rho_mid_Fe_H, axis=0)
    med_rho_low_Fe_H = np.median(rho_low_Fe_H, axis=0)

    for i in range(len(rho_high_mass)):
        ax.plot(r_mid, rho_high_Fe_H[i], lw=1, alpha=0.5, c="tab:red")
        ax.plot(r_mid, rho_mid_Fe_H[i], lw=1, alpha=0.5, c="tab:orange")
        ax.plot(r_mid, rho_low_Fe_H[i], lw=1, alpha=0.5, c="tab:blue")

    ax.plot(r_mid, med_rho_high_Fe_H, c="tab:red",
             label=r"$%.1f < [{\rm Fe/H}]$" % FE_H_HIGH)
    ax.plot(r_mid, med_rho_mid_Fe_H, c="tab:orange",
             label=r"$%.1f < [{\rm Fe/H}] < %.1f$" % (FE_H_LOW, FE_H_HIGH))
    ax.plot(r_mid, med_rho_low_Fe_H, c="tab:blue",
             label=r"$[{\rm Fe/H} < %.1f]$" % FE_H_LOW)

    ax.plot(r_mid, med_rho_high_mass + med_rho_low_mass, "--", c="k")
    ax.set_xlabel(r"$r/R_{\rm vir}$")
    ax.set_ylabel(r"$\rho\ (M_\odot\,{\rm kpc}^{-3})$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, None)
    ax.legend(loc="upper right", fontsize=17)

    fig.savefig(plot_dir+"Fe_H_density_%s.png" % suite)

    fig, ax = plt.subplots()
    r_mid = np.sqrt(r_bins[1:]*r_bins[:-1])

    med_Fe_H_high_mass = np.median(mean_Fe_H_high_mass, axis=0)
    med_Fe_H_low_mass = np.median(mean_Fe_H_low_mass, axis=0)

    rho_cut = 0

    for i in range(len(mean_Fe_H_high_mass)):
        ok = (rho_high_mass[i] > rho_cut) & (~np.isinf(mean_Fe_H_high_mass[i]))
        ax.plot(r_mid[ok], mean_Fe_H_high_mass[i][ok],
                lw=1, alpha=0.5, c="tab:red")
        ok = (rho_low_mass[i] > rho_cut) & (~np.isinf(mean_Fe_H_low_mass[i]))
        ax.plot(r_mid[ok], mean_Fe_H_low_mass[i][ok],
                lw=1, alpha=0.5, c="tab:blue")

    ok = (med_rho_high_mass > rho_cut) & (~np.isinf(med_Fe_H_high_mass))
    ax.plot(r_mid[ok], med_Fe_H_high_mass[ok], c="tab:red",
             label=r"$M_{\rm sub,infall}/M_{\rm host,infall} > 0.15$")
    ok = (med_rho_low_mass > rho_cut) & (~np.isinf(med_Fe_H_low_mass))
    ax.plot(r_mid[ok], med_Fe_H_low_mass[ok], c="tab:blue",
             label=r"$M_{\rm sub,infall}/M_{\rm host,infall} \leq 0.15$")

    total_rho = med_rho_low_mass + med_rho_high_mass
    ok = (total_rho > rho_cut) & (~np.isinf(med_Fe_H_low_mass)) & (~np.isinf(med_Fe_H_high_mass))
    mean_Fe_H = (med_rho_high_mass[ok]*med_Fe_H_high_mass[ok] +
                 med_rho_low_mass[ok]*med_Fe_H_low_mass[ok]) / total_rho[ok]
    ax.plot(r_mid[ok], mean_Fe_H, c="k")
    ax.set_xlabel(r"$r/R_{\rm vir}$")
    ax.set_ylabel(r"$\langle[{\rm Fe/H}]\rangle$")
    ax.set_xscale("log")
    ax.legend(loc="upper right", fontsize=17)
    ylo, yhi = plt.ylim()
    plt.ylim(ylo, yhi + 0.15*(yhi - ylo))
    fig.savefig(plot_dir+"average_Fe_H_%s.png" % suite)

    fig, ax = plt.subplots()

    order = np.argsort(mpeak_sub)

    m_star_gal_sum = np.cumsum(m_star_gal[order])
    m_star_halo_sum = np.cumsum(m_star_halo[order])

    ax.plot(mpeak_sub[order], m_star_halo_sum/n_hosts_used,
            c="tab:red", label=r"${\rm halo\ stars}$")
    ax.plot(mpeak_sub[order], m_star_gal_sum/n_hosts_used,
            c="tab:blue", label=r"${\rm satellite\ stars}$")
    ax.set_xlabel(r"$M_{\rm sub,peak}$")
    ax.set_ylabel(r"$M_\star(<M_{\rm sub,peak})/M_{*,{\rm tot}}$")
    ax.legend(loc="lower right", fontsize=17)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e4, None)

    fig.savefig(plot_dir+"star_contribution_cdf_%s.png" % suite)

if __name__ == "__main__": main()
