import numpy as np
import matplotlib.pyplot as plt
import symlib

BUGGY_HOSTS = [6, 9, 10, 16, 17, 31, 36, 37, 40, 42, 43]

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
    
    # Each of these returns a list of arrays, one for each subhalo.

    # owner is 0 if the particle is owned by the subhalo and a positive
    # integer otherwise.
    owner = symlib.read_particles(part_info, sim_dir, last_snap, "ownership")
    # True if the particle has been accreted onto the subhalo in the current
    # snapshot and False otherwise.
    valid = symlib.read_particles(part_info, sim_dir, last_snap, "valid")
    x = symlib.read_particles(part_info, sim_dir, last_snap, "x")

    # The cut that determines whether a tracked core is still intact
    core_r = np.sqrt(np.sum(c["x"]**2, axis=2))
    mp = param["mp"]/param["h100"]
    # core tracking must have started, core needs 32 particles, and it can't 
    # overlapping qith the host's center.
    core_intact = c["ok"] & (c["m_bound"] > 32*mp) & (core_r > c["r50_bound"])

    rf, in_halo = [None]*len(h), [None]*len(h)

    for i in target_subs:
        # Only work with particles where this flag is true.
        ok = (owner[i] == 0) & valid[i]
        # Correct the units. There's an analogous function for velocities.
        x_i = symlib.set_units_x(x[i], h_cmov[0,-1], scale[-1], param)

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

def density_profile(sim_dir, r_bins, r, mp_star, in_halo, include_sub,
                    target_subs=None):
    """ This function just generates a density profile. Nothing too
    interesting.
    """
    h, _ = symlib.read_subhalos(sim_dir)

    if target_subs is None:
        target_subs = np.arange(1, len(h))

    r_host = h["rvir"][0,-1]
    r_bins = r_bins*r_host # Locally', we'll convert out fo normalized units
    V_bins = (r_bins[1:]**3 - r_bins[:-1]**3)*4*np.pi/3

    n_bins = len(r_bins) - 1
    mass_hist = np.zeros(n_bins)

    for i in target_subs:
        if not include_sub[i]: continue

        ok = in_halo[i] & (mp_star[i] > 0)
        m, _ = np.histogram(r[i][ok], bins=r_bins, weights=mp_star[i][ok])
        mass_hist += m

    return mass_hist / V_bins

def main():
    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    suite = "SymphonyMilkyWay"

    param = symlib.simulation_parameters(suite)
    n_hosts = symlib.n_hosts(suite)

    # Set up density profiles
    n_bins = 30
    r_bins = np.logspace(-2, 0, n_bins + 1) # In units of Rvir(z=0)
    rho_high_mass = np.zeros(n_bins)
    rho_low_mass = np.zeros(n_bins)

    # Our galaxy-halo model
    gal_halo = symlib.GalaxyHaloModel(
        symlib.UniverseMachineMStar(),
        symlib.Jiang2019RHalf(),
        symlib.PlummerProfile()
    )

    # Arrays needed to store plot data
    m_star_halo = []
    m_star_gal = []
    mpeak_sub = []
    n_hosts_used = 0
    for i in range(n_hosts):
        if i in BUGGY_HOSTS: continue
        print("Host %2d/%d" % (i, n_hosts))

        # This function lets you loop over all the subhalo directories without
        # needing to know them by name.
        sim_dir = symlib.get_host_directory(base_dir, suite, i)

        # Shared information on ownership and offsets needed to
        # decompress files.
        part_info = symlib.ParticleInfo(sim_dir)

        h, hist = symlib.read_subhalos(sim_dir)

        # Set to whatever halos you want to analyze. I'd suggest starting with
        # a couple small ones for debugging and then changing it back to
        # everything later. Don't include halo 0.
        target_subs = np.arange(1, len(h))

        # mp_star is a list of arrays giving the stellar mass of each particle.
        mp_star, _ = symlib.tag_stars(sim_dir, gal_halo, 
                                      target_subs=target_subs)

        # Example analysis function. rf and in_halo are both lists of arrays
        # rf gives the final radius of particles at z=0 relative to the host
        # center and in_halo is true if the particle is in the host's halo
        # and false otherwise. rf will be -1 for invalid particles.
        rf, in_halo = final_radii(sim_dir, part_info,
                                  target_subs=target_subs)

        is_high_mass = hist["merger_ratio"] > 0.05

        # Keep track of contribution to total halo/satellite mass.
        for j in target_subs:
            m_star_halo.append(np.sum(mp_star[j][in_halo[j]]))
            m_star_gal.append(np.sum(mp_star[j][~in_halo[j]]))
            mpeak_sub.append(hist["mpeak"][j])

        # Compute density profiles, add to running averages.
        rho_high_mass += density_profile(
            sim_dir, r_bins, rf, mp_star, in_halo, is_high_mass,
            target_subs=target_subs)
        rho_low_mass += density_profile(
            sim_dir, r_bins, rf, mp_star, in_halo, ~is_high_mass,
            target_subs=target_subs)
    
        n_hosts_used += 1
        # Uncomment this if you only want to look at one halo.
        break

    m_star_halo = np.array(m_star_halo)
    m_star_gal = np.array(m_star_gal)
    mpeak_sub = np.array(mpeak_sub)

    # Plotting nonsense from here on
    fig, ax = plt.subplots()
    r_mid = np.sqrt(r_bins[1:]*r_bins[:-1])
    
    ax.plot(r_mid, rho_high_mass/n_hosts_used, c="tab:red",
             label=r"$M_{\rm sub,infall}/M_{\rm host,infall} > 0.05$")
    ax.plot(r_mid, rho_low_mass/n_hosts_used, c="tab:blue",
             label=r"$M_{\rm sub,infall}/M_{\rm host,infall} \leq 0.05$")
    ax.plot(r_mid, (rho_high_mass + rho_low_mass)/n_hosts_used, c="k")
    ax.set_xlabel(r"$r\ ({\rm kpc})$")
    ax.set_ylabel(r"$\rho\ (M_\odot\,{\rm kpc}^{-3})$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, None)
    ax.legend(loc="upper right", fontsize=17)

    fig.savefig("../plots/stellar_halo/average_density.png")

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

    fig.savefig("../plots/stellar_halo/star_contribution_cdf.png")

if __name__ == "__main__": main()
