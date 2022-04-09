import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import lib
from colossus.cosmology import cosmology
from colossus.halo import mass_so
import star_tagging
import matplotlib.colors as mpl_colors
import relax
import numpy.random as random

BASE_DIR_FMT = "../tmp_data/Halo%03d"
MP = 2.8e5
EPS = 170e-3

params = { 'flat': True, 'H0': 70.0, 'Om0': 0.286,
           'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.95 }
cosmo = cosmology.setCosmology('myCosmo', params)
h100 = params["H0"]/100.0

gal_halo = star_tagging.GalaxyHaloModel(
    star_tagging.UniverseMachineMStar(),
    star_tagging.Jiang2019RHalf(),
    star_tagging.PlummerProfile()
)

BASE_HALO = 169

####################
# Helper functions #
####################

def calc_dt_orbit(snap_start):
    scale = lib.scale_factors()
    z = 1/scale - 1
    age = cosmo.age(z)
    t_orbit = mass_so.dynamicalTime(z, "vir", "orbit")
    dt = (age - age[snap_start])
    dt_orbit = dt/t_orbit[snap_start]

    return dt, dt_orbit

def subhalo_indices(host_id):
    host_id = 169
    base_dir = BASE_DIR_FMT % host_id

    h_idx, scale_start, scale_infall = {
        169: (4, 0.25, 0.9),
        222: (6, 0.3, 0.95),
        523: (10, 0.35, 0.9),
        719: (1, 0.4, 0.85)
    }[host_id]

    scale = lib.scale_factors()
    snap_start = np.searchsorted(scale, scale_start)
    snap_infall = np.searchsorted(scale, scale_infall)

    dt, dt_orbit = calc_dt_orbit(snap_start)
    comp_snaps = np.searchsorted(dt_orbit, [0.1, 0.25, 0.5, 1.0, 2.0, 4.0])
    comp_snaps = comp_snaps[comp_snaps <= snap_infall]
    
    return base_dir, h_idx, snap_start, comp_snaps

######################
# Plotting functions #
######################

def plot_energy_levels(base_dir, h_idx, snap_start, comp_snaps):
    _, m = lib.read_mergers(base_dir)
    m["rvir"] *= 1e3
    
    dt, dt_orbit = calc_dt_orbit(snap_start)
    scale = lib.scale_factors()
    all_snaps = np.array([snap_start] + list(comp_snaps))

    x, v, _ = lib.read_part_file(base_dir, snap_start, h_idx, ["x", "v"])
    n_max = len(x)

    a0, h0 = scale[snap_start], m[h_idx,snap_start]
    x, v, idx = star_tagging.clean_particles(x, v, h0, h100, a0)
    
    R0 = star_tagging.RadialEnergyRanking(
        MP/h100, x, v, idx, n_max,
    )

    ones = np.ones(1000)
    rhalf = gal_halo.r_half_model.r_half(
        h0["rvir"]*a0/h100*ones, h0["cvir"]*ones, 1/a0 - 1)

    rhalf_rvir = np.mean(rhalf/h0["rvir"]/a0)
    
    colors = [pc("r"), pc("o"), pc("b")]

    orbit_names = ["00", "01", "02", "05", "10", "20", "40"]
    
    for i_snap, snap in enumerate(all_snaps):
        plt.figure()
        
        x, v, _ = lib.read_part_file(base_dir, snap, h_idx, ["x", "v"])
        x, v, idx = star_tagging.clean_particles(
            x, v, m[h_idx,snap], h100, scale[snap])
        
        r = np.sqrt(np.sum(x**2, axis=1))/R0.rmax
        vr = np.sum(x*v, axis=1)/r/R0.vmax/R0.rmax

        rank_max = len(colors) - 1
        bulk = ((R0.ranks[idx] == star_tagging.NIL_RANK) |
                (R0.ranks[idx] > 2*rank_max))

        extent=[-2, np.log10(h0["rvir"]*a0/R0.rmax), -2, 2]
        vmin, vmax = 1, 50
        binsize=50
        
        plt.hexbin(np.log10(r[bulk]), vr[bulk], gridsize=100, cmap="Greys",
                   norm=mpl_colors.LogNorm(vmin=vmin, vmax=vmax),
                   extent=extent)
        
        for i in range(rank_max, -1, -1):
            ok = (R0.ranks[idx] == 2*i) | (R0.ranks[idx] == 2*i + 1)
            plt.plot(np.log10(r[ok]), vr[ok], "o", c=colors[i])
            
        plt.title((r"$\Delta t = %.1f\ {\rm Gyr} = " +
                   r"%.2f\cdot t_{\rm orbit,vir}$") %
                  (dt[snap], dt_orbit[snap]))

        plt.plot([np.log10(rhalf_rvir)]*2, [-2, 2], "--", c="k")
        
        plt.xlim(extent[0], extent[1])
        plt.ylim(extent[2], extent[3])
        plt.ylabel(r"$v/V_{\rm max,0}$")
        plt.xlabel(r"$r/R_{\rm max,0}$")

        plt.savefig("../plots/star_tagging_plots/energy_levels.%03d.%d.%s.png" %
                    (BASE_HALO, h_idx, orbit_names[i_snap]))
        
def plot_relaxation(base_dir, h_idx, snap_start, comp_snaps):
    all_snaps = np.array([snap_start] + list(comp_snaps))
    i_tag = 1
    
    _, m = lib.read_mergers(base_dir)
    m["rvir"] *= 1e3
    
    dt, dt_orbit = calc_dt_orbit(snap_start)
    scale = lib.scale_factors()

    x, v, _ = lib.read_part_file(base_dir, snap_start, h_idx, ["x", "v"])    
    n_max = len(x)

    a0, h0 = scale[snap_start], m[h_idx,snap_start]
    x, v, idx = star_tagging.clean_particles(x, v, h0, h100, a0)
    R0 = star_tagging.RadialEnergyRanking(MP/h100, x, v, idx, n_max)

    snap_tag = all_snaps[i_tag]
    a_tag, h_tag = scale[snap_tag], m[h_idx,snap_tag]
    x_tag, v_tag, _ = lib.read_part_file(base_dir, snap_tag, h_idx, ["x", "v"])
    x_tag, v_tag, idx_tag = star_tagging.clean_particles(
        x_tag, v_tag, h_tag, h100, a_tag)

    R0.load_particles(x_tag, v_tag, idx_tag)    
    rhalf_tag = R0.ranked_halfmass_radius()
    t_relax_tag = R0.ranked_relaxation_time(MP/h100, EPS/h100)

    colors = [pc("r"), pc("o"), pc("g"), pc("b"), pc("p"), pc("k")]
    
    for i in range(i_tag+1, len(all_snaps)):
        snap = all_snaps[i]
        a, h = scale[snap], m[h_idx,snap]
        x, v, _ = lib.read_part_file(base_dir, snap, h_idx, ["x", "v"])
        x, v, idx = star_tagging.clean_particles(x, v, h, h100, a)

        R0.load_particles(x, v, idx)
        
        rhalf = R0.ranked_halfmass_radius()
        dt_t_relax = (dt[all_snaps[i]] - dt[snap_start]) / t_relax_tag

        plt.plot(dt_t_relax, rhalf/rhalf_tag, c=colors[i-1-i_tag],
                 label=r"$\Delta t = %.2f\,t_{\rm orbit}$" % dt_orbit[snap])
        plt.plot(dt_t_relax, rhalf/rhalf_tag, "o", c=colors[i-1-i_tag])

    plt.xscale("log")
        
    ylo, yhi = plt.ylim()
    xlo, xhi = plt.xlim()
    plt.ylim(ylo, yhi)
    plt.xlim(xlo, xhi)
    plt.legend(loc="upper left", fontsize=17)

    plt.fill_between([0.1, xhi], [ylo, ylo], [yhi,yhi], alpha=0.25, color="k")
    plt.plot([xlo, xhi], [1, 1], "--", lw=2, c="k")
    
    plt.xlabel(r"$\Delta t/t_{\rm relax}(E_0)$")
    plt.ylabel(r"$R_{1/2}(E_0)/R_{1/2,0}(E_0)$")

    plt.savefig("../plots/star_tagging_plots/relaxation.%03d.%d.png" %
                (BASE_HALO, h_idx))
        
def plot_tagging_steps(base_dir, h_idx, snap_start, comp_snaps):
    # Boiler plate to set up ranks
    i_tag = 1
    i_final_1 = -2
    i_final_2 = -1
    
    scale = lib.scale_factors()
    z = 1/scale - 1
    dt, dt_orbit = calc_dt_orbit(snap_start)
    
    all_snaps = np.array([snap_start] + list(comp_snaps))
    _, m = lib.read_mergers(base_dir)
    m["rvir"] *= 1e3
    a0, h0 = scale[snap_start], m[h_idx,snap_start]

    x, v, _ = lib.read_part_file(base_dir, snap_start, h_idx, ["x", "v"])
    n_max = len(x)
    x, v, idx = star_tagging.clean_particles(x, v, h0, h100, a0)
    R0 = star_tagging.RadialEnergyRanking(MP/h100, x, v, idx, n_max)
    
    # Load in particles at the time of evaluation
    snap = all_snaps[i_tag]
    a, h = scale[snap], m[h_idx,snap]
    x, v, _ = lib.read_part_file(base_dir, snap, h_idx, ["x", "v"])
    x, v, idx = star_tagging.clean_particles(x, v, h, h100, a)

    snap_final_1 = all_snaps[i_final_1]
    a_final_1, h_final_1 = scale[snap_final_1], m[h_idx,snap_final_1]
    x_final_1, v_final_1, _ = lib.read_part_file(
        base_dir, snap_final_1, h_idx, ["x", "v"])
    x_final_1, v_final_1, idx_final_1 = star_tagging.clean_particles(
        x_final_1, v_final_1, h_final_1, h100, a_final_1)
    
    R0.load_particles(x_final_1, v_final_1, idx_final_1)
    r_final_1 = np.sqrt(np.sum(R0.x**2, axis=1))


    snap_final_2 = all_snaps[i_final_2]
    a_final_2, h_final_2 = scale[snap_final_2], m[h_idx,snap_final_2]
    x_final_2, v_final_2, _ = lib.read_part_file(
        base_dir, snap_final_2, h_idx, ["x", "v"])
    x_final_2, v_final_2, idx_final_2 = star_tagging.clean_particles(
        x_final_2, v_final_2, h_final_2, h100, a_final_2)
    
    R0.load_particles(x_final_2, v_final_2, idx_final_2)
    r_final_2 = np.sqrt(np.sum(R0.x**2, axis=1))
    
    R0.load_particles(x, v, idx)
    r = np.sqrt(np.sum(R0.x**2, axis=1))

    c_min = 0.5
    c_max = 0.9

    bins = 10**np.linspace(np.log10(np.min(r)), np.log10(np.max(r)), 200)

    plt.figure(1)
    
    for i in range(R0.n_max, -1, -1):
        ok = np.where(R0.ranks[idx] == i)[0]

        c_val = c_min + (c_max - c_min)*i/R0.n_max
        weights = np.ones(len(ok))*MP/h100
        plt.hist(r[ok], weights=weights, cumulative=True, histtype="step",
                 lw=2.5, color=pc("k", c_val), bins=bins)

    random.seed(4)
    colors = [pc("r"), pc("o"), pc("g"), pc("b"), pc("p")]
    n_example = len(colors)
    
    r_half_model = star_tagging.Jiang2019RHalf()
    m_star_model = star_tagging.UniverseMachineMStar()
    profile_model = star_tagging.PlummerProfile()

    rvir = h0["rvir"]*scale[snap]/h100
    r_half = [r_half_model.r_half(cvir=h0["cvir"], z=z[snap], rvir=rvir)
              for _ in range(n_example)]
    mpeak = np.max(m[h_idx,:]["mvir"])
    m_star = [m_star_model.m_star(mpeak=mpeak/h100, z=z[snap])
              for _ in range(n_example)]
    
    for i in range(n_example):
        plt.figure(1)
        prof = profile_model.m_enc(m_star[i], r_half[i], bins)        
        plt.plot(bins, prof, c=colors[i])

        mp_star = R0.set_mp_star(rvir, profile_model, r_half[i], m_star[i])

        plt.figure(2)
        E_mid = (R0.E_edges[1:] + R0.E_edges[:-1]) / 2
        ok = R0.mp_star_table > 1e-1
        plt.plot(E_mid[ok], R0.mp_star_table[ok], c=colors[i])
        plt.plot(E_mid[ok], R0.mp_star_table[ok], "o", c=colors[i])

        plt.figure(3)
        #plt.hist(r, weights=mp_star[idx], bins=bins, cumulative=True,
        #         color=colors[i], histtype="step",lw=3)
        m_star_diff, edges = np.histogram(r, weights=mp_star[idx], bins=bins)
        m_star_sum = np.zeros(len(edges))
        m_star_sum[1:] = np.cumsum(m_star_diff)
        plt.plot(edges, m_star_sum, color=colors[i])
        plt.plot(bins, prof, "--", c="k", lw=2)

        plt.figure(4)

        m_star_diff_final_1, edges = np.histogram(
            r_final_1, weights=mp_star[idx_final_1], bins=bins)
        m_star_sum_final_1 = np.zeros(len(edges))
        m_star_sum_final_1[1:] = np.cumsum(m_star_diff_final_1)

        m_star_diff_final_2, edges = np.histogram(
            r_final_2, weights=mp_star[idx_final_2], bins=bins)
        m_star_sum_final_2 = np.zeros(len(edges))
        m_star_sum_final_2[1:] = np.cumsum(m_star_diff_final_2)
        
        plt.plot(edges, m_star_sum_final_1, c=colors[i])
        plt.plot(edges, m_star_sum_final_2, "--", c=colors[i])
        plt.plot(edges, m_star_sum, "--", c="k", lw=2)
        
        
    plt.figure(1)
    plt.ylim(MP/h100, None)
    plt.xlim(np.min(r), np.max(r))

    plt.xlabel(r"$r\ ({\rm kpc})$")
    plt.ylabel(r"$M(<r)\ (M_\odot)$")
    plt.xscale("log")
    plt.yscale("log")

    plt.savefig("../plots/star_tagging_plots/raw_mass_profiles.%03d.%d.png" %
                (BASE_HALO, h_idx))

    plt.figure(2)
    plt.xlabel(r"$E_0/(m_{p,{\rm DM}}\, V_{\rm max}^2)$")
    plt.ylabel(r"$m_{\rm p,\star}(E_0)$")
    plt.yscale("log")
    plt.ylim(3e-1, 3e5)

    plt.savefig("../plots/star_tagging_plots/mp_star.%03d.%d.png" %
                (BASE_HALO, h_idx))
    
    plt.figure(3)
    plt.xlim(np.min(r), np.max(r))

    plt.xlabel(r"$r\ ({\rm kpc})$")
    plt.ylabel(r"$M_\star(<r)\ (M_\odot)$")
    plt.xlim(0, 6)

    plt.savefig("../plots/star_tagging_plots/profile_matching.%03d.%d.png" %
                (BASE_HALO, h_idx))

    plt.figure(4)
    plt.xlabel(r"$r\ ({\rm kpc})$")
    plt.ylabel(r"$M_\star(<r)\ (M_\odot)$")
    plt.xlim(0, 6)
    plt.plot([], [], pc("r"), label=r"$2\cdot t_{\rm orbit,vir}$")
    plt.plot([], [], "--", c=pc("r"), label=r"$4\cdot t_{\rm orbit,vir}$")
    plt.legend(loc="lower right", fontsize=17)

    plt.savefig("../plots/star_tagging_plots/profile_evolution.%03d.%d.png" %
                (BASE_HALO, h_idx))

    
def main():
    palette.configure(False)

    base_dir, h_idx, snap_start, comp_snaps = subhalo_indices(BASE_HALO)

    #plot_energy_levels(base_dir, h_idx, snap_start, comp_snaps)
    #plot_relaxation(base_dir, h_idx, snap_start, comp_snaps)
    plot_tagging_steps(base_dir, h_idx, snap_start, comp_snaps)
    
    plt.show()
    
if __name__ == "__main__":
    main()
