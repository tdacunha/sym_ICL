import symlib
import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import os.path as path
from colossus.cosmology import cosmology
import scipy.interpolate as interpolate
import read_mw_sats

#invalid_hosts = [36]
invalid_hosts = []

USE_TEX = True

SUITE = "SymphonyMilkyWay"
OUT_DIR = "../plots/core_tracking"
BASE_DIR = "/sdf/home/p/phil1/ZoomIns"
invalid_hosts = []
#SUITE = "SymphonyMilkyWay"
#OUT_DIR = "../plots/core_plots/MilkyWay"
#invalid_hosts = []

def plot_mass_loss_old():
    palette.configure(USE_TEX)

    out_dir = OUT_DIR
    suite = SUITE
    sim_dir = symlib.get_host_directory(base_dir, suite, 0)

    param = symlib.simulation_parameters(suite)
    a = symlib.scale_factors(sim_dir)
    
    h, hist = symlib.read_subhalos(param, sim_dir)
    h = symlib.set_units_halos(h, a, param)
    c = symlib.read_cores(sim_dir)

    targets = np.arange(1, len(h), dtype=int)

    fig, ax = plt.subplots(2, figsize=(7, 14), sharex=True)
    
    for i_sub in targets:
        ax[0].cla()
        ax[1].cla()

        ax[0].plot(a, h["mvir"][0], "--", c="k")
        
        okh, mvir = h["ok"][i_sub], h["mvir"][i_sub]
        rh = np.sqrt(np.sum(h["x"][i_sub]**2, axis=1))
    
        ax[0].plot(a[okh], mvir[okh], pc("r"), label=r"$M_{\rm vir,RS}$")
        ax[1].plot(a[okh], rh[okh], pc("r"), label=r"$r_{\rm Rockstar}$")
        okc = c["ok"][i_sub]
        m_bound, m_tidal = c["m_bound"][i_sub], c["m_tidal_bound"][i_sub]
        r_half = c["r50_bound"][i_sub]
        rc = np.sqrt(np.sum(c["x"][i_sub]**2, axis=1))
        ax[0].plot(a[okc], m_bound[okc], pc("o"), label=r"$M_{\rm bound}$")
        ax[0].plot(a[okc], m_tidal[okc], pc("b"), label=r"$M_{\rm tidal}$")
        ax[1].plot(a[okc], rc[okc], pc("b"), label=r"$r_{\rm core}$")
        ax[1].plot(a[okh], rh[okh], ":", c=pc("r"))
        ax[1].plot(a, h[0]["rvir"], "--", c="k")
        ax[1].plot(a[okc], r_half[okc], "--", c=pc("a"))

        mp = param["mp"]/param["h100"]
        ax[0].set_ylim(mp, None)
        ax[0].fill_between([0, 1], [mp]*2, [25*mp]*2, color="k", alpha=0.2)
        ax[0].set_yscale("log")
        ax[1].set_yscale("log")
        ax[1].set_xlim(0, 1)
        ax[1].set_ylim(np.min(h[0]["rvir"]), None)
        
        ax[1].set_xlabel(r"$a(t)$")
        ax[1].set_ylabel(r"$r\ ({\rm kpc})$")
        ax[0].set_ylabel(r"$M\ (M_\odot)$")

        ax[0].legend(loc="upper right", fontsize=16)
        ax[1].legend(loc="lower left", fontsize=16)

        fig.savefig(path.join(out_dir, "ml_%04d.png" % i_sub))

def plot_mass_loss():
    palette.configure(USE_TEX)

    out_dir = path.join(OUT_DIR, "sub_mass_loss")
    suite = SUITE
    base_dir = BASE_DIR
    sim_dir = symlib.get_host_directory(base_dir, suite, 0)

    param = symlib.simulation_parameters(suite)
    a = symlib.scale_factors(sim_dir)
    
    h, hist = symlib.read_subhalos(sim_dir)
    c = symlib.read_cores(sim_dir)

    targets = np.arange(1, len(h), dtype=int)

    # fig, ax = plt.subplots(2, figsize=(7, 14), sharex=True)
    fig, ax = plt.subplots(figsize=(6,6))
    
    for i_sub in targets:
        if i_sub != 14: continue
        print("i_sub =", i_sub)
        """
        ax[0].cla()
        ax[1].cla()

        ax[0].plot(a, h["mvir"][0], "--", c="k")
        
        okh, mvir = h["ok"][i_sub], h["mvir"][i_sub]
        rh = np.sqrt(np.sum(h["x"][i_sub]**2, axis=1))
    
        ax[0].plot(a[okh], mvir[okh], pc("r"), label=r"$M_{\rm vir,RS}$")
        ax[1].plot(a[okh], rh[okh], pc("r"), label=r"$r_{\rm Rockstar}$")
        okc = c["ok"][i_sub]
        m_bound, m_tidal = c["m_bound"][i_sub], c["m_tidal_bound"][i_sub]
        r_half = c["r50_bound"][i_sub]
        rc = np.sqrt(np.sum(c["x"][i_sub]**2, axis=1))
        ax[0].plot(a[okc], m_bound[okc], pc("o"), label=r"$M_{\rm bound}$")
        ax[0].plot(a[okc], m_tidal[okc], pc("b"), label=r"$M_{\rm tidal}$")
        ax[1].plot(a[okc], rc[okc], pc("b"), label=r"$r_{\rm core}$")
        ax[1].plot(a[okh], rh[okh], ":", c=pc("r"))
        ax[1].plot(a, h[0]["rvir"], "--", c="k")
        ax[1].plot(a[okc], r_half[okc], "--", c=pc("a"))

        mp = param["mp"]/param["h100"]
        ax[0].set_ylim(mp, None)
        ax[0].fill_between([0, 1], [mp]*2, [25*mp]*2, color="k", alpha=0.2)
        ax[0].set_yscale("log")
        ax[1].set_yscale("log")
        ax[1].set_xlim(0, 1)
        ax[1].set_ylim(np.min(h[0]["rvir"]), None)
        
        ax[1].set_xlabel(r"$a(z)$")
        ax[1].set_ylabel(r"$r\ ({\rm kpc})$")
        ax[0].set_ylabel(r"$M\ (M_\odot)$")

        ax[0].legend(loc="upper right", fontsize=16)
        ax[1].legend(loc="lower left", fontsize=16)

        fig.savefig(path.join(out_dir, "ml_%04d.png" % i_sub))
        """

        ax.cla()
        
        okh, mvir = h["ok"][i_sub], h["mvir"][i_sub]
        rh = np.sqrt(np.sum(h["x"][i_sub]**2, axis=1))
    
        okc = c["ok"][i_sub]
        m_bound, m_tidal = c["m_bound"][i_sub], c["m_tidal_bound"][i_sub]
        r_half = c["r50_bound"][i_sub]
        rc = np.sqrt(np.sum(c["x"][i_sub]**2, axis=1))

        ax.plot(a[okh], rh[okh], pc("r"), label=r"${\rm Rockstar}$")
        ax.plot(a[okc], rc[okc], pc("b"), label=r"${\rm Particle-tracking}$")
        ax.plot(a[okh], rh[okh], ":", c=pc("r"))
        ax.plot(a, h[0]["rvir"], "--", c="k", label=r"$R_{\rm vir,host}$")

        ax.set_yscale("log")
        ax.set_xlim(0, 1)
        ax.set_ylim(np.min(h[0]["rvir"]), None)
        
        ax.set_xlabel(r"$a(t)$")
        ax.set_ylabel(r"$r\ ({\rm kpc})$")

        ax.legend(loc="lower right", fontsize=18)

        fig.savefig(path.join(out_dir, "ml_%04d.pdf" % i_sub))


def mass_function():
    print("mass_function_compare")
    palette.configure(USE_TEX)

    base_dir = BASE_DIR

    gs_kw = {"height_ratios": [3, 1.5]}
    fig_main, axs_main = plt.subplots(nrows=2, sharex=True,
                                      gridspec_kw=gs_kw)
    ax_main, rat_main = axs_main[0], axs_main[1]
    fig_c, axs_c = plt.subplots(nrows=2, sharex=True,
                                      gridspec_kw=gs_kw)
    ax_c, rat_c = axs_c[0], axs_c[1]

    fig_rs, axs_rs = plt.subplots(nrows=2, sharex=True,
                                      gridspec_kw=gs_kw)
    ax_rs, rat_rs = axs_rs[0], axs_rs[1]

    fig_hr, axs_hr = plt.subplots(nrows=2, sharex=True,
                                      gridspec_kw=gs_kw)
    ax_hr, rat_hr = axs_hr[0], axs_hr[1]


    #fig_hr, ax_hr = plt.subplots(sharex=True)
    #fig_c, ax_c = plt.subplots(sharex=True)
    #fig_rs, ax_rs = plt.subplots(sharex=True)

    figs = [fig_main, fig_hr, fig_c, fig_rs]
    axes = [ax_main, ax_hr, ax_c, ax_rs]
    rats = [rat_main, rat_hr, rat_c, rat_rs]

    #radii_mults = [1, 1/5]
    #radii_labels = [r"$r < R_{\rm vir}$", r"$r < R_{\rm vir}/5$"]
    #radii_colors = [pc("r"), pc("b")]
    
    radii_mults = [1, 1/2, 1/4, 1/8][:3]
    radii_labels = [r"$r < R_{\rm vir}$", r"$r < R_{\rm vir}/2$",
                    r"$r < R_{\rm vir}/4$", r"$r < R_{\rm vir}/8$"][:3]
    radii_colors = [pc("r"), pc("o"), pc("b"), pc("p")]
    
    #radii_mults = [1, 1/10]
    #radii_labels = [r"$r < R_{\rm vir}$", r"$r < R_{\rm vir}/10$"]
    #radii_colors = [pc("r"), pc("b")]

    suites = ["SymphonyMilkyWay", "SymphonyMilkyWayHR", "SymphonyMilkyWayLR"]
    
    saved_ratio = [[None, None] for _ in range(4)]
    saved_r = [[None, None] for _ in range(4)]

    fig_suites = [[0, 0], [1, 1], [1, 2], [1, 2]]
    fig_types = [[1, 0], [1, 0], [1, 1], [0, 0]]
    fig_suffix = ["main", "hr", "c", "rs"]
    fig_labels = [[r"${\rm Particle{-}tracking}$", r"${\rm Rockstar}$"],
                  [r"${\rm Particle{-}tracking}$", r"${\rm Rockstar}$"],
                  [r"${\rm High{-}res}$", r"${\rm Fiducial}$"], 
                  [r"${\rm High{-}res}$", r"${\rm Fiducial}$"]]

    for i_suite in range(len(suites)):
        suite = suites[i_suite]
        param = symlib.simulation_parameters(suite)

        cs = [None]*symlib.n_hosts(suite)
        hs, hists = [None]*len(cs), [None]*len(cs)

        for halo_type in range(2):
            ratios = []
            rads = []
            
            n_host = 0
            
            for i_host in range(symlib.n_hosts(suite)):    
                print("%s: %2d/%2d" % (suite, i_host, symlib.n_hosts(suite) - 1))

                sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

                if halo_type == 0:
                    hs[i_host], hists[i_host] = symlib.read_subhalos(sim_dir)
                    cs[i_host] = symlib.read_cores(sim_dir)
                h, hist, c = hs[i_host], hists[i_host], cs[i_host]

                if halo_type == 0:
                    mvir, x = h["mvir"], h["x"]
                    ok = c["ok_rs"] #rockstar_ok(h, c, param)
                elif halo_type == 1:
                    mvir, x = c["m_bound"], c["x"]
                    ok = c["ok"]
                    #mvir, x, ok = combine_cores_rockstar(h, c, param)
                else:
                    assert(0)

                ratio = hist["mpeak_pre"] / h["mvir"][0,-1]
                r = np.sqrt(np.sum(x[:,-1,:]**2, axis=1)) / h["rvir"][0,-1]

                ok[0,-1] = False

                ratios.append(ratio[ok[:,-1]])
                rads.append(r[ok[:,-1]])
                n_host += 1

            r = np.hstack(rads)
            ratio = np.hstack(ratios)

            for i_fig in range(len(figs)):
                suites_to_plot = fig_suites[i_fig]
                types_to_plot = fig_types[i_fig]

                for i_curve in range(len(suites_to_plot)):
                    if suites_to_plot[i_curve] != i_suite: continue
                    if types_to_plot[i_curve] != halo_type: continue

                    saved_r[i_fig][i_curve] = r
                    saved_ratio[i_fig][i_curve] = ratio

    for i_fig in range(len(figs)):
        fig, ax, rat = figs[i_fig], axes[i_fig], rats[i_fig]
        suites_to_plot = fig_suites[i_fig]
        types_to_plot = fig_types[i_fig]

        funcs = [[None, None] for _ in range(len(radii_mults))]

        for i_curve in range(len(suites_to_plot)):
            
            r = saved_r[i_fig][i_curve]
            ratio = saved_ratio[i_fig][i_curve]

            for i_mult in range(len(radii_mults)):
                ok = r < radii_mults[i_mult]

                ls = "-" if i_curve == 0 else "--"
                color = radii_colors[i_mult]

                ax.plot(np.sort(ratio[ok]),
                        np.arange(np.sum(ok))[::-1]/n_host,
                        ls=ls, c=color)                            

                funcs[i_mult][i_curve] = interpolate.interp1d(
                    np.log10(np.sort(ratio[ok])),
                    np.arange(np.sum(ok))[::-1],
                    bounds_error=False,
                    fill_value=-1
                )

                if i_curve == 0: continue

                rat_eval = np.linspace(-5, -1, 200)

                n_min = 5
                n_num = funcs[i_mult][0](rat_eval)
                n_den = funcs[i_mult][1](rat_eval)
                ok = (n_num > n_min) & (n_den > n_min)
                rat.plot([1e-5, 1], [1, 1], "--", lw=1, c="k")
                rat.plot(10**rat_eval[ok], n_num[ok]/n_den[ok], c=color)


    for i_fig in range(len(figs)):
        fig, ax, rat = figs[i_fig], axes[i_fig], rats[i_fig]
        suffix = fig_suffix[i_fig]

        if i_fig == 0 or i_fig == 1:
            for i_rad in range(len(radii_mults)):
                ax.plot([], [], label=radii_labels[i_rad],
                        c=radii_colors[i_rad])
        ax.plot([], [], "-", c=pc("a"), label=fig_labels[i_fig][0])
        ax.plot([], [], "--", c=pc("a"), label=fig_labels[i_fig][1])
        #ax.plot([], [], c=pc("r"), label=r"${\rm Rockstar}$")
        #ax.plot([], [], c=pc("b"), label=r"${\rm Particle-Tracking}$")
        #ax.plot([], [], "-", c=pc("a"), label=r"$r<300\ {\rm kpc}$")
        #ax.plot([], [], "--", c=pc("a"), label=r"$r<30\ {\rm kpc}$")

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e-5, 0.2)
        ax.set_ylim(1e-1, 2e3)
        rat.set_xlabel(r"$\mu=m_{\rm peak,sub}/M_{\rm host}$")
        if i_fig == 0 or i_fig == 1:
            rat.set_ylabel(r"$N_{\rm pt}/N_{\rm Rockstar}$")
        else:
            rat.set_ylabel(r"$N_{\rm HR}/N_{\rm fid}$")
        ax.set_ylabel(r"$N(>\mu)$")
        if i_fig == 1:
            rat.set_ylim(0.6, 4)
        if i_fig > 0: 
            rat.set_xlim(1e-5, 4e-2)

        ax.legend(loc="lower left", fontsize=17)
    
        fig.savefig(path.join(OUT_DIR, "core_shmf_%s.pdf" % suffix))

def mass_function_resolution():
    print("mass_function_resolution")
    palette.configure(USE_TEX)

    base_dir = BASE_DIR

    fig_vir, ax_vir = plt.subplots()
    fig_peak, ax_peak = plt.subplots()

    radii_lims = [50, 100, 250]
    radii_labels = [r"$r < 50\ {\rm kpc}$", r"$r < 100\ {\rm kpc}$",
                    r"$r < 250\ {\rm kpc}$",]
    radii_colors = [pc("r"), pc("o"), pc("b")]

    suites = ["SymphonyMilkyWayHR", "SymphonyMilkyWayLR", "SymphonyMilkyWay"]
    suite_ls = ["-", "--", ":"]
    
    for i_suite in range(len(suites)):
        suite = suites[i_suite]
        param = symlib.simulation_parameters(suite)

        nvirs_c, npeaks_c, rads_c = [], [], []
        nvirs_h, npeaks_h, rads_h = [], [], []

        for i_host in range(symlib.n_hosts(suite)):

            print("%s: %2d/%2d" % (suite, i_host, symlib.n_hosts(suite) - 1))
            sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
            
            h, hist = symlib.read_subhalos(sim_dir)
            c = symlib.read_cores(sim_dir)

            mvir, x, ok = combine_cores_rockstar(h, c, param)
            r = np.sqrt(np.sum(x[:,-1,:]**2, axis=1))
            ok[0,-1] = False
            h["ok"][0,-1] = False
            ok_h =h["ok"]

            mp = param["mp"] / param["h100"]
            nvirs_c.append(mvir[ok[:,-1],-1] / mp)
            npeaks_c.append(hist["mpeak_pre"][ok[:,-1]] / mp)
            rads_c.append(r[ok[:,-1]])
            nvirs_h.append(mvir[ok_h[:,-1],-1] / mp)
            npeaks_h.append(hist["mpeak_pre"][ok_h[:,-1]] / mp)
            rads_h.append(r[ok_h[:,-1]])

        nvir_c = np.hstack(nvirs_c)
        npeak_c = np.hstack(npeaks_c)
        r_c = np.hstack(rads_c)
        nvir_h = np.hstack(nvirs_h)
        npeak_h = np.hstack(npeaks_h)
        r_h = np.hstack(rads_h)

        n_max = max(np.max(npeak_c), np.max(npeak_h))
        mass_bins = 10**np.linspace(0, np.log10(n_max)*1.01, 100)
        
        for i_rad in range(len(radii_lims)):
            ok_h = r_h < radii_lims[i_rad]
            ok_c = r_c < radii_lims[i_rad]

            N_vir_c, N_edges = np.histogram(nvir_c[ok_c], bins=mass_bins)
            N_vir_c = np.cumsum(N_vir_c[::-1])[::-1]
            N_vir_h, _ = np.histogram(nvir_h[ok_h], bins=mass_bins)
            N_vir_h = np.cumsum(N_vir_h[::-1])[::-1]
            N_peak_c, _ = np.histogram(npeak_c[ok_c], bins=mass_bins)
            N_peak_c = np.cumsum(N_peak_c[::-1])[::-1]
            N_peak_h, _ = np.histogram(npeak_h[ok_h], bins=mass_bins)
            N_peak_h = np.cumsum(N_peak_h[::-1])[::-1]

            N_min = symlib.n_hosts(suite)*1.5
            N_max = 1e11/mp
            ok_vir = ((N_vir_c > N_min) & (N_vir_h > N_min) &
                      (N_edges[:-1] < N_max))
            ok_peak = ((N_peak_c > N_min) & (N_peak_h > N_min) & 
                       (N_edges[:-1] < N_max) & (N_edges[:-1] > 300))

            N_mids = np.sqrt(N_edges[:-1] * N_edges[1:])

            c = radii_colors[i_rad]
            ls = suite_ls[i_suite]
            ax_vir.plot(N_mids[ok_vir], 1 - N_vir_h[ok_vir]/N_vir_c[ok_vir],
                        ls=ls, c=c)
            ax_peak.plot(N_mids[ok_peak],
                         1 - N_peak_h[ok_peak]/N_peak_c[ok_peak],
                         ls=ls, c=c)

    for i_rad in range(len(radii_labels)):
        ax_vir.plot([], [], c=radii_colors[i_rad], label=radii_labels[i_rad])
        ax_peak.plot([], [], c=radii_colors[i_rad], label=radii_labels[i_rad])

    ax_vir.set_xscale("log")
    ax_peak.set_xscale("log")
    ax_vir.set_xlim(30, None)
    ax_peak.set_xlim(30, None)
    ax_vir.set_ylim(0, 1)
    ax_vir.set_ylim(0, 1)

    ax_vir.set_xlabel(r"$N_{\rm vir,sub}$")
    ax_peak.set_xlabel(r"$N_{\rm peak,sub}$")
    ax_vir.set_ylabel(r"$f_{\rm missing}$")
    ax_peak.set_ylabel(r"$f_{\rm missing}$")

    ax_vir.legend(loc="upper left", fontsize=17)
    ax_peak.legend(loc="upper left", fontsize=17)
    fig_vir.savefig("../plots/core_plots/f_missing_resolution_vir.png")
    fig_peak.savefig("../plots/core_plots/f_missing_resolution_peak.png")



def t_survive(cosmo, a, infall_snap, ok):
    z = 1/a - 1
    t = cosmo.age(z)

    t_first = t[infall_snap]
    #print(infall_snap, np.arange(len(ok), dtype=int)[ok])
    t_last = np.max(t[ok])


    return t_last - t_first

def survival_time():
    palette.configure(USE_TEX)

    base_dir = BASE_DIR
    suite = SUITE

    param = symlib.simulation_parameters(suite)
    mp = param["mp"]/param["h100"]

    cut_edges = 32*mp*np.array([10, 100, 1e3, 1e4, 1e5])
    n_merger = np.zeros(len(cut_edges)-1, dtype=int)
    n_disrupt = np.zeros(len(cut_edges)-1, dtype=int)
    t_ratio_merger = [[] for _ in range(len(cut_edges)-1)]
    t_ratio_disrupt = [[] for _ in range(len(cut_edges)-1)]

    c_param = symlib.colossus_parameters(param)
    cosmo = cosmology.setCosmology(suite, c_param)

    for i_halo in range(symlib.n_hosts(suite)):
        if i_halo in invalid_hosts:
            continue

        sim_dir = symlib.get_host_directory(base_dir, suite, i_halo)

        a = symlib.scale_factors(sim_dir)
        h, hist = symlib.read_subhalos(sim_dir)
        c = symlib.read_cores(sim_dir)

        m_vir = h["mvir"][1:,:]
        m_peak = np.max(h["mvir"], axis=1)[1:]
        okh = h["ok"][1:,:] & (m_vir > 32*mp)
        
        r_h = np.sqrt(np.sum(h["x"][1:,:]**2, axis=2))
        r_c = np.sqrt(np.sum(c["x"][1:,:]**2, axis=2))
        r_half = c["r50_bound"][1:,:]
        m_bound = c["m_bound"][1:,:]
        m_tidal = c["m_tidal"][1:,:]
        okc = c["ok"][1:,:] & (m_bound > 32*mp) & (r_c > r_half)
        merger_snap = hist["merger_snap"][1:]

        n_okh = np.sum(okh, axis=1)
        n_okc = np.sum(okc, axis=1)
        is_gone = ((~okc) & (~okh))[:,-1]

        is_merger = hist["merger_ratio"][1:] > 0.05

        print("%2d" %  i_halo)

        for i in range(len(cut_edges)-1):
            in_cut = ((cut_edges[i] < m_peak) & (m_peak <= cut_edges[i+1]) & 
                      (n_okc > 1) & (n_okh > 1))
            cut_mergers = np.where(in_cut & is_merger & is_gone)[0]
            cut_disrupts = np.where(in_cut & (~is_merger) & is_gone)[0]

            t_ratio_merger_i = np.zeros(len(cut_mergers))
            t_ratio_disrupt_i = np.zeros(len(cut_disrupts))
            for j in range(len(cut_disrupts)):
                k = cut_disrupts[j]
                t_rockstar = t_survive(cosmo, a, merger_snap[k], okh[k,:])
                t_core = t_survive(cosmo, a, merger_snap[k], okc[k,:])
                if t_rockstar == 0 or t_core == 0: continue
                t_ratio_disrupt_i[j] = t_core/t_rockstar

            t_ratio_disrupt[i].append(t_ratio_disrupt_i)


            n_merger[i] += len(cut_mergers)
            n_disrupt[i] += np.sum(cut_disrupts)

    fig, ax = plt.subplots()

    colors = [pc("r"), pc("o"), pc("b"), pc("p")]
    Np = [r"3\times 10^2", r"3\times 10^3",
          r"3\times 10^4", r"3\times 10^5"]
    for i in range(len(cut_edges)-1):
        if n_disrupt[i] < 5: continue
        t_ratio_disrupt[i] = np.hstack(t_ratio_disrupt[i])
        print(n_disrupt[i], np.mean(t_ratio_disrupt[i]), np.std(t_ratio_disrupt[i]))
        ax.hist(np.log10(t_ratio_disrupt[i]), bins=300, lw=3, color=colors[i],
                histtype="step", density=True, cumulative=True, range=(-1, 2),
                label=r"$N_{\rm peak} \approx %s$" % Np[i])
        
    ax.set_xlim(-1, 2)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$\log_{10}(\tau = t_{\rm core-tracking}/t_{\rm Rockstar})$")
    ax.set_ylabel(r"$N(<\tau)$")
    ax.legend(loc="lower right", fontsize=17)

    ax.grid()
    fig.savefig("../plots/core_plots/t_ratio_hist.png")
    fig.savefig(path.join(OUT_DIR, "t_ratio_hist.png"))

    print(n_merger)
    print(n_disrupt)

def stitching_errors():
    print("stitching_errors")
    palette.configure(USE_TEX)

    base_dir = BASE_DIR
    suite = SUITE

    if suite == "SymphonyMilkyWay":
        bins = [(10**2.5, 10**3), (10**3.5, 10**4), (10**4.5, 10**5)]
        colors = [pc("b"), pc("o"), pc("r")]
        labels = [r"$10^{2.5} < N_{\rm peak} < 10^{3}$",
                  r"$10^{3.5} < N_{\rm peak} < 10^{4}$",
                  r"$10^{4.5} < N_{\rm peak} < 10^{5}$"]
    else:
        bins = [(10**2.5, 10**3), (10**3.5, 10**4), (10**4.5, 10**5)]
        colors = [pc("b"), pc("o"), pc("r")]
        labels = [r"$10^{2.5} < N_{\rm peak} < 10^{3}$",
                  r"$10^{3.5} < N_{\rm peak} < 10^{4}$",
                  r"$10^{4.5} < N_{\rm peak} < 10^{5}$"]        

    n_mass_bins = len(bins)

    param = symlib.simulation_parameters(suite)
    
    n_host = 0
    n_bins = 200
    r_bins = np.logspace(-3, 0, n_bins+1)

    N_tot_now = [np.zeros(n_bins) for _ in range(n_mass_bins)]
    N_err_now = [np.zeros(n_bins) for _ in range(n_mass_bins)]

    N_tot_disrupt = [np.zeros(n_bins) for _ in range(n_mass_bins)]
    N_err_disrupt = [np.zeros(n_bins) for _ in range(n_mass_bins)]

    for i_halo in range(symlib.n_hosts(suite)):
        if i_halo in invalid_hosts:
            continue
        print(i_halo)

        sim_dir = symlib.get_host_directory(base_dir, suite, i_halo)

        a = symlib.scale_factors(sim_dir)
        h, hist = symlib.read_subhalos(sim_dir)
        c = symlib.read_cores(sim_dir)

        mp = param["mp"]/param["h100"]

        host = h[0]
        h, c, hist = h[1:], c[1:], hist[1:]

        m_vir, m_peak = h["mvir"][:,-1], hist["mpeak_pre"]
        
        for j in range(n_mass_bins):
            r_h = np.sqrt(np.sum(h["x"][:,-1]**2, axis=1))
            r_c = np.sqrt(np.sum(c["x"][:,-1]**2, axis=1))
            m_bound = c["m_bound"][:,-1]
            
            ok = (h["ok"][:,-1] & (m_peak/mp >= bins[j][0]) &
                  (m_peak/mp >= bins[j][0]))
            is_err = h["ok"][:,-1] & (~c["ok_rs"][:,-1])

            r0 = host["rvir"][-1]

            n_tot, _ = np.histogram(r_h[ok]/r0, bins=r_bins)
            n_err, _ = np.histogram(r_h[ok & is_err]/r0, bins=r_bins)

            N_tot_now[j] += n_tot
            N_err_now[j] += n_err

            r_c = np.sqrt(np.sum(c["x"]**2, axis=2))
            r_h = np.sqrt(np.sum(h["x"]**2, axis=2))
            snap = np.arange(h.shape[1])

            r_peri = np.zeros(len(h))
            for k in range(len(h)):
                if np.sum(c["ok_rs"][k]) == 0:
                    r_peri[k] = np.min(r_h[k,h["ok"][k]])
                else:
                    r_peri[k] = np.min(r_h[k,c["ok_rs"][k]])

            ok = h["ok"]
            mass_ok = (m_peak/mp >= bins[j][0]) & (m_peak/mp < bins[j][1])

            last_snap = np.zeros(len(h), dtype=int)
            for i in range(len(last_snap)):
                if np.sum(ok[i]) == 0 or h["ok"][i,-1]:
                    last_snap[i] = -1
                else:
                    last_snap[i] = np.max(snap[ok[i]])

            is_disrupt = (last_snap != -1) & mass_ok
            is_err = ~c["ok_rs"][np.arange(len(c), dtype=int),last_snap]

            r0 = host["rvir"][last_snap]
            
            n_tot_disrupt, _ = np.histogram(
                r_peri[is_disrupt]/r0[is_disrupt], bins=r_bins)
            n_err_disrupt, _ = np.histogram(
                r_peri[is_disrupt & is_err]/r0[is_disrupt & is_err],
                bins=r_bins)
            N_tot_disrupt[j] += n_tot_disrupt
            N_err_disrupt[j] += n_err_disrupt

    N_tot_now = [np.cumsum(N_tot_now[j]) for j in range(n_mass_bins)]
    N_err_now = [np.cumsum(N_err_now[j]) for j in range(n_mass_bins)]
    N_tot_disrupt = [np.cumsum(N_tot_disrupt[j]) for j in range(n_mass_bins)]
    N_err_disrupt = [np.cumsum(N_err_disrupt[j]) for j in range(n_mass_bins)]
    
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()

    for j in range(n_mass_bins):
        ok_now = N_tot_now[j] > 5
        ok_disrupt = N_tot_disrupt[j] > 5

        right_bin_now = r_bins[1:][ok_now]
        right_bin_disrupt = r_bins[1:][ok_disrupt]
        
        f_err_now = N_err_now[j][ok_now]/N_tot_now[j][ok_now]
        f_err_disrupt = N_err_disrupt[j][ok_disrupt]/N_tot_disrupt[j][ok_disrupt]

        ax1.plot(right_bin_now, f_err_now, c=colors[j], label=labels[j])
        ax2.plot(right_bin_disrupt, f_err_disrupt, c=colors[j], label=labels[j])

    ax1.set_xscale("log")
    ax1.set_xlabel(r"$r/R_{\rm vir}$")
    ax1.set_ylabel(r"$f_{\rm err,RS}(< r/R_{\rm vir})$")
    ax1.legend(loc="upper right", fontsize=17)
    fig1.savefig(path.join(OUT_DIR, "core_stitching_now.pdf"))

    ax2.set_xscale("log")
    ax2.set_xlabel(r"$r_{\rm peri}/R_{\rm vir}$")
    ax2.set_ylabel(r"$f_{\rm err,disrupt,RS}(< r_{\rm peri}/R_{\rm vir})$")
    ax2.legend(loc="upper right", fontsize=17)
    fig2.savefig(path.join(OUT_DIR, "core_stitching_disrupt.pdf"))

def core_ok(c, param):
    mp = param["mp"]/param["h100"]
    r = np.sqrt(np.sum(c["x"]**2, axis=2))
    ok = c["ok"] & (c["m_bound"] > 32*mp) & (r > c["r50_bound"])
    ok[0,:] = False
    return ok
    
def rockstar_stitching_error(h, c, param):
    c_ok = core_ok(c, param)
    dr =  np.sqrt(np.sum((h["x"] - c["x"])**2, axis=2))
    return dr > c["r50_bound"]
    
def previous_stitching_err(s_err, c_ok, h_ok):
    # Couldn't figure out how to do this one without the loop, ha ha
    out = np.zeros(h_ok.shape, dtype=bool)
    for hi in range(len(h_ok)):
        curr_s_err = False
        for si in range(len(h_ok[0])):
            if c_ok[hi,si] and h_ok[hi,si]:
                curr_s_err = s_err[hi,si]
            out[hi,si] = curr_s_err
    return out

def rockstar_ok(h, c, param):
    c_ok = core_ok(c, param)
    s_err = rockstar_stitching_error(h, c, param)
    h_ok = h["ok"]
    curr_s_err = c_ok & s_err & h_ok
    h_ok[curr_s_err] = False
    prev_s_err = previous_stitching_err(s_err, c_ok, h_ok)
    h_ok[prev_s_err] = False
    return h_ok

def combine_cores_rockstar(h, c, param):
    ok_h = rockstar_ok(h, c, param)
    ok_c = core_ok(c, param)
    
    ok = ok_h | ok_c
    x = np.zeros((h.shape[0], h.shape[1], 3))
    mvir = np.zeros(h.shape)

    for snap in range(h.shape[1]):
        mvir[ok_h[:,snap],snap] = h["mvir"][ok_h[:,snap],snap]
        x[ok_h[:,snap],snap] = h["x"][ok_h[:,snap],snap]

        mvir[ok_c[:,snap],snap] = c["m_bound"][ok_c[:,snap],snap]
        x[ok_c[:,snap],snap] = c["x"][ok_c[:,snap],snap]

    return mvir, x, ok
        

def ufd_frequency():
    palette.configure(USE_TEX)

    base_dir = BASE_DIR
    suite = SUITE

    param = symlib.simulation_parameters(suite)
    mp = param["mp"]/param["h100"]
    n_hosts = symlib.n_hosts(suite)

    mw_sats = read_mw_sats.read()

    r_bins = np.linspace(0, 50, 200)
    
    fig, ax = plt.subplots()
    names = ["1e80", "1e85", "1e90", "1e95"]
    exponents = ["8", "8.5", "9", "9.5"]
    colors = [pc("r"), pc("o"), pc("b"), pc("p")]
    mpeak_cuts = [10**8, 10**8.5, 10**9, 10**9.5]
    
    c_arrays = [[] for _ in range(len(names))]
    h_arrays = [[] for _ in range(len(names))]

    for i in range(n_hosts):
        if i in invalid_hosts: continue
        print("%2d/%d" % (i, n_hosts))

        sim_dir = symlib.get_host_directory(base_dir, suite, i)
        h, hist = symlib.read_subhalos(sim_dir)
        c = symlib.read_cores(sim_dir)

        h_r = np.sqrt(np.sum(h["x"][1:,-1,:]**2, axis=1))
        h_ok = rockstar_ok(h, c, param)[1:,-1]
        c_r = np.sqrt(np.sum(c["x"][1:,-1,:]**2, axis=1))
        c_ok = core_ok(c, param)[1:,-1]

        for j in range(len(mpeak_cuts)):
            mpeak_ok = hist["mpeak_pre"][1:] > mpeak_cuts[j]
        
            h_n = np.cumsum(np.histogram(h_r[h_ok & mpeak_ok],
                                         bins=r_bins)[0])

            c_n = np.cumsum(np.histogram(c_r[c_ok & mpeak_ok],
                                         bins=r_bins)[0])
            all_r = np.hstack([c_r[c_ok & mpeak_ok], 
                               h_r[h_ok & (~c_ok) & mpeak_ok]])
            #c_n = np.cumsum(np.histogram(all_r, bins=r_bins)[0])

            c_arrays[j].append(c_n)
            h_arrays[j].append(h_n)

    ok_1 = mw_sats["class"] >= 4
    ok_2 = mw_sats["class"] >= 3
    n_1 = np.cumsum(np.histogram(mw_sats["r"][ok_1], bins=r_bins)[0])
    n_2 = np.cumsum(np.histogram(mw_sats["r"][ok_2], bins=r_bins)[0])
    ax.plot(r_bins[:-1], n_2, c=pc("a"),
               label=r"${\rm probable\ dwarfs}$")
    ax.plot(r_bins[:-1], n_1, c=pc("k"),
               label=r"${\rm confirmed\ dwarfs}$")

    for i in range(len(mpeak_cuts)):
        ax.plot(r_bins[:-1], np.mean(c_arrays[i], axis=0),
                "-", c=colors[i],
                label=r"$M_{lim} = 10^{%s}\ M_\odot$" % exponents[i])
        ax.plot(r_bins[:-1], np.mean(h_arrays[i], axis=0),
                "--", c=colors[i])

    ax.set_xlabel(r"$r\ ({\rm kpc})$")
    ax.set_ylabel(r"$N(<r)$")
    ax.legend(loc="upper left", fontsize=16)
    ax.set_yscale("log")
    ax.set_ylim(0.2, 100)
    ax.set_xlim(0, 50)
    fig.savefig(path.join(OUT_DIR, "ufd_freq_%s.png") % names[i])

def nfw_mass_frac(r_rvir, cvir):
    r_rs = r_rvir * cvir
    def M_enc(x): return np.log(1 + x) + 1/(1 + x) - 1
    return M_enc(r_rs) / M_enc(cvir)

def radius_cdf():
    print("radius_cdf")
    palette.configure(USE_TEX)

    base_dir = BASE_DIR

    gs_kw = {"height_ratios": [3, 1.5]}
    fig_main, axs_main = plt.subplots(nrows=2, sharex=True,
                                      gridspec_kw=gs_kw)
    ax_main, rat_main = axs_main[0], axs_main[1]
    fig_c, axs_c = plt.subplots(nrows=2, sharex=True,
                                      gridspec_kw=gs_kw)
    ax_c, rat_c = axs_c[0], axs_c[1]

    fig_rs, axs_rs = plt.subplots(nrows=2, sharex=True,
                                      gridspec_kw=gs_kw)
    ax_rs, rat_rs = axs_rs[0], axs_rs[1]

    fig_hr, axs_hr = plt.subplots(nrows=2, sharex=True,
                                      gridspec_kw=gs_kw)
    ax_hr, rat_hr = axs_hr[0], axs_hr[1]

    figs = [fig_main, fig_hr, fig_c, fig_rs]
    axes = [ax_main, ax_hr, ax_c, ax_rs]
    rats = [rat_main, rat_hr, rat_c, rat_rs]

    npeak_cuts = [(10**2.5, 10**3), (10**3.5, 10**4), (10**4.5, 10**5)]
    #npeak_cuts = [0.3e3, 3e3, 30e3]
    cut_labels = [
        r"$10^{2.5} < N_{\rm peak} < 10^3$",
        r"$10^{3.5} < N_{\rm peak} < 10^4$",
        r"$10^{4.5} < N_{\rm peak} < 10^5$",
        #r"$N_{\rm peak} > 3\times 10^2$",
        #r"$N_{\rm peak} > 3\times 10^3$",
        #r"$N_{\rm peak} > 3\times 10^4$",
    ]
    cut_colors = [pc("r"), pc("o"), pc("b")]

    suites = ["SymphonyMilkyWay", "SymphonyMilkyWayHR", "SymphonyMilkyWayLR"]
    
    fig_suites = [[0, 0], [1, 1], [1, 2], [1, 2]]
    fig_types = [[1, 0], [1, 0], [1, 1], [0, 0]]
    fig_suffix = ["main", "hr", "c", "rs"]
    fig_labels = [[r"${\rm Particle{-}tracking}$", r"${\rm Rockstar}$"],
                  [r"${\rm Particle{-}tracking}$", r"${\rm Rockstar}$"],
                  [r"${\rm High{-}res}$", r"${\rm Low{-}res}$"], 
                  [r"${\rm High{-}res}$", r"${\rm Low{-}res}$"]]

    for i_suite in range(len(suites)):
        suite = suites[i_suite]
        param = symlib.simulation_parameters(suite)
        mp = param["mp"]/param["h100"]

        cs = [None]*symlib.n_hosts(suite)
        hs, hists = [None]*len(cs), [None]*len(cs)

        for halo_type in range(2):
            npeaks = []
            rads = []
            
            n_host = 0
            
            for i_host in range(symlib.n_hosts(suite)):
                print("%s: %2d/%2d" % (suite, i_host, symlib.n_hosts(suite) - 1))

                sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

                if halo_type == 0:
                    hs[i_host], hists[i_host] = symlib.read_subhalos(sim_dir)
                    cs[i_host] = symlib.read_cores(sim_dir)
                h, hist, c = hs[i_host], hists[i_host], cs[i_host]

                if halo_type == 0:
                    mvir, x = h["mvir"], h["x"]
                    ok = c["ok_rs"] #rockstar_ok(h, c, param)
                elif halo_type == 1:
                    mvir, x = c["m_bound"], c["x"]
                    ok = c["ok"]
                    #mvir, x, ok = combine_cores_rockstar(h, c, param)
                else:
                    assert(0)

                r = np.sqrt(np.sum(x[:,-1,:]**2, axis=1)) / h["rvir"][0,-1]
                npeak = hist["mpeak_pre"]/mp

                ok[0,-1] = False

                npeaks.append(npeak[ok[:,-1]])
                rads.append(r[ok[:,-1]])
                n_host += 1

            r = np.hstack(rads)
            npeak = np.hstack(npeaks)

            for i_fig in range(len(figs)):
                fig, ax = figs[i_fig], axes[i_fig]
                suites_to_plot = fig_suites[i_fig]
                types_to_plot = fig_types[i_fig]

                for i_curve in range(len(suites_to_plot)):
                    if suites_to_plot[i_curve] != i_suite: continue
                    if types_to_plot[i_curve] != halo_type: continue

                    for i_mult in range(len(npeak_cuts)):
                        ok = (npeak > npeak_cuts[i_mult][0]) & (r < 1) & (npeak< npeak_cuts[i_mult][1])

                        ls = "-" if i_curve == 0 else "--"
                        color = cut_colors[i_mult]

                        ax.plot(np.sort(r[ok]),
                                np.arange(np.sum(ok))/np.sum(ok),
                                ls=ls, c=color)

    for i_fig in range(len(figs)):
        fig, ax, rat = figs[i_fig], axes[i_fig], rats[i_fig]
        suffix = fig_suffix[i_fig]

        for i_rad in range(len(npeak_cuts)):
            ax.plot([], [], label=cut_labels[i_rad],
                    c=cut_colors[i_rad])

        ax.plot([], [], "-", c=pc("a"), label=fig_labels[i_fig][0])
        ax.plot([], [], "--", c=pc("a"), label=fig_labels[i_fig][1])
        
        r_part, m_part = particle_mass_profile()

        ax.plot(r_part, m_part, ":", c="k",
                label=r"$M(<r/R_{\rm vir})/M_{\rm vir}$")


        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        rat.set_xlabel(r"$r/R_{\rm vir}$")
        ax.set_ylabel(r"$N(<r/R_{\rm vir})$")

        ax.legend(loc="upper left", fontsize=17)
    
        fig.savefig(path.join(OUT_DIR, "radius_cdf_%s.pdf" % suffix))

def particle_mass_profile():
    r_hi, scaled_rho = np.loadtxt("tables/SymphonyMilkyWay_density_profile.txt").T
    rho = scaled_rho/r_hi**2
    
    r_edge = np.zeros(len(r_hi)+1)
    r_edge[1:] = r_hi
    vol = r_edge[1:]**3 - r_edge[:-1]**3
    mass = np.zeros(len(r_edge))
    mass[1:] = vol*rho
    mass = np.cumsum(mass)

    f_mass = interpolate.interp1d(np.log10(r_edge), mass)
    m_vir = f_mass(0)

    #print(r_edge)
    #print(mass/m_vir)

    return r_edge, mass / m_vir

def mass_loss():
    print("mass_loss")
    palette.configure(USE_TEX)
    
    base_dir = BASE_DIR
    suite = SUITE

    if suite == "SymphonyMilkyWay":
        bins = [1e8, 1e9, 1e10, 1e11]
        colors = [pc("b"), pc("o"), pc("r")]
        labels = [r"$10^{8} < M_{\rm peak}/M_\odot < 10^{9}$",
                  r"$10^{9} < M_{\rm peak}/M_\odot < 10^{10}$",
                  r"$10^{10} < M_{\rm peak}/M_\odot < 10^{11}$"]

    else:
        bins = [1e7, 1e8, 1e9, 1e10, 1e11]
        colors = [pc("p"), pc("b"), pc("o"), pc("r")]
        labels = [r"$10^{7} < M_{\rm peak}/M_\odot < 10^{8}$",
                  r"$10^{8} < M_{\rm peak}/M_\odot < 10^{9}$",
                  r"$10^{9} < M_{\rm peak}/M_\odot < 10^{10}$",
                  r"$10^{10} < M_{\rm peak}/M_\odot < 10^{11}$"]
    n_bins = len(bins) - 1

    r_bins = np.linspace(0, 1, 101)
    h_n = np.zeros((len(bins)-1, len(r_bins)-1))
    c_n = np.zeros((len(bins)-1, len(r_bins)-1))

    
    n = len(bins) - 1
    ratio_h, ratio_c = [[] for _ in range(n)], [[] for _ in range(n)]
    r_h, r_c = [[] for _ in range(n)], [[] for _ in range(n)]

    for i in range(symlib.n_hosts(suite)):
        if i in invalid_hosts: continue
        print("%2d/%d" % (i, symlib.n_hosts(suite)))

        sim_dir = symlib.get_host_directory(base_dir, suite, i)
        param = symlib.simulation_parameters(sim_dir)
        h, hist = symlib.read_subhalos(sim_dir)
        c = symlib.read_cores(sim_dir)

        h_r = np.sqrt(np.sum(h["x"][1:,-1,:]**2, axis=1))
        h_ok = rockstar_ok(h, c, param)[1:,-1]
        c_r = np.sqrt(np.sum(c["x"][1:,-1,:]**2, axis=1))
        c_ok = core_ok(c, param)[1:,-1]
        
        for j in range(len(bins) - 1):
            mpeak_ok = ((hist["mpeak_pre"] > bins[j]) &
                        (hist["mpeak_pre"] < bins[j+1]))[1:]

            ratio = (h["mvir"][:,-1]/hist["mpeak_pre"])[1:]
            ratio_h[j].append(ratio[mpeak_ok & h_ok])
            r_h[j].append(h_r[mpeak_ok & h_ok] / h["rvir"][0,-1])

            r = np.copy(h_r)
            ratio[c_ok] = (c["m_bound"][:,-1]/hist["mpeak_pre"])[1:][c_ok]
            r[c_ok] = c_r[c_ok]

            ratio_c[j].append(ratio[mpeak_ok & (c_ok | h_ok)])
            r_c[j].append(r[mpeak_ok & (c_ok | h_ok)] / h["rvir"][0,-1])

    for i in range(len(ratio_h)):
        ratio_h[i] = np.hstack(ratio_h[i])
        ratio_c[i] = np.hstack(ratio_c[i])
        r_h[i] = np.hstack(r_h[i])
        r_c[i] = np.hstack(r_c[i])

    r_lim = [0.2, 0.5, 1]
    rlim_name = ["02", "05", "10"]
    r_bins = 10**np.linspace(-4, 0, 301)
    low, med, high = np.zeros((n_bins, 300)), np.zeros((n_bins, 300)), np.zeros((n_bins, 300))

    n_min = 8

    for i in range(len(r_lim)):
        fig, ax = plt.subplots()

        for j in range(len(bins) - 1):
            ok_h = r_h[j] < r_lim[i]
            ok_c = r_c[j] < r_lim[i]

            ax.hist(ratio_h[j][ok_h], color=colors[j], cumulative=True,
                    density=True, histtype="step", lw=3, bins=r_bins, ls="--")
            ax.hist(ratio_c[j][ok_c], color=colors[j], cumulative=True,
                    density=True, histtype="step", lw=3, bins=r_bins,
                    label=labels[j])
        
            ax.set_xlabel(r"$m=M(z=0)/M_{\rm peak}$")
            ax.set_ylabel(r"$N(<m)$")
            ax.legend(loc="upper left", fontsize=18)

            """
            for k in range(len(r_bins) - 1):
                if i != 0: break
                ok = (r[j] < r_bins[k]) & ok_c
                n = np.sum(ok)
                if n < n_min:
                    low[k], med[k], high[k] = -1, -1, -1
                else:
                    low[k] = np.percentile(r[j][ok], 50 - 68/2)
                    med[k] = np.percentile(r[j][ok], 50)
                    high[k] = np.percentile(r[j][ok], 50 + 68/2)
            """

            
        ax.set_xscale("log")
        ax.set_xlim(1e-3, 1)
        fig.savefig(path.join(OUT_DIR, "mass_loss_high_r%s.png") %
                    rlim_name[i])

def mass_frac_cdf():
    print("radius_cdf")
    palette.configure(USE_TEX)

    base_dir = BASE_DIR

    fig_main, ax_main = plt.subplots()
    fig_hr, ax_hr = plt.subplots()
    fig_c, ax_c = plt.subplots()
    fig_rs, ax_rs = plt.subplots()

    figs = [fig_main, fig_hr, fig_c, fig_rs]
    axes = [ax_main, ax_hr, ax_c, ax_rs]

    npeak_cuts = [0.3e3, 3e3, 30e3]
    cut_labels = [
        r"$N_{\rm peak} > 3\times 10^2$",
        r"$N_{\rm peak} > 3\times 10^3$",
        r"$N_{\rm peak} > 3\times 10^4$",
    ]
    cut_colors = [pc("r"), pc("o"), pc("b")]

    suites = ["SymphonyMilkyWay", "SymphonyMilkyWayHR", "SymphonyMilkyWayLR"]
    
    fig_suites = [[0, 0], [1, 1], [1, 2], [1, 2]]
    fig_types = [[1, 0], [1, 0], [1, 1], [0, 0]]
    fig_suffix = ["main", "hr", "c", "rs"]
    fig_labels = [[r"${\rm Core-tracking}$", r"${\rm Rockstar}$"],
                  [r"${\rm Core-tracking}$", r"${\rm Rockstar}$"],
                  [r"${\rm High-res}$", r"${\rm Low-res}$"], 
                  [r"${\rm High-res}$", r"${\rm Low-res}$"]]

    for i_suite in range(len(suites)):
        suite = suites[i_suite]
        param = symlib.simulation_parameters(suite)
        mp = param["mp"]/param["h100"]

        cs = [None]*symlib.n_hosts(suite)
        hs, hists = [None]*len(cs), [None]*len(cs)

        for halo_type in range(2):
            npeaks = []
            rads = []
            
            n_host = 0
            
            for i_host in range(symlib.n_hosts(suite)):
                print("%s: %2d/%2d" % (suite, i_host, symlib.n_hosts(suite) - 1))

                sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

                if halo_type == 0:
                    hs[i_host], hists[i_host] = symlib.read_subhalos(sim_dir)
                    cs[i_host] = symlib.read_cores(sim_dir)
                h, hist, c = hs[i_host], hists[i_host], cs[i_host]

                if halo_type == 0:
                    mvir, x = h["mvir"], h["x"]
                    ok = rockstar_ok(h, c, param)
                elif halo_type == 1:
                    mvir, x, ok = combine_cores_rockstar(h, c, param)
                else:
                    assert(0)

                r = np.sqrt(np.sum(x[:,-1,:]**2, axis=1)) / h["rvir"][0,-1]
                npeak = hist["mpeak_pre"]/mp

                ok[0,-1] = False

                npeaks.append(npeak[ok[:,-1]])
                rads.append(r[ok[:,-1]])
                n_host += 1

            r = np.hstack(rads)
            npeak = np.hstack(npeaks)

            for i_fig in range(len(figs)):
                fig, ax = figs[i_fig], axes[i_fig]
                suites_to_plot = fig_suites[i_fig]
                types_to_plot = fig_types[i_fig]

                for i_curve in range(len(suites_to_plot)):
                    if suites_to_plot[i_curve] != i_suite: continue
                    if types_to_plot[i_curve] != halo_type: continue

                    for i_mult in range(len(npeak_cuts)):
                        ok = (npeak > npeak_cuts[i_mult]) & (r < 1)

                        ls = "-" if i_curve == 0 else "--"
                        color = cut_colors[i_mult]

                        ax.plot(np.sort(r[ok]),
                                np.arange(np.sum(ok))/np.sum(ok),
                                ls=ls, c=color)

    for i_fig in range(len(figs)):
        fig, ax = figs[i_fig], axes[i_fig]
        suffix = fig_suffix[i_fig]

        for i_rad in range(len(npeak_cuts)):
            ax.plot([], [], label=cut_labels[i_rad],
                    c=cut_colors[i_rad])

        ax.plot([], [], "-", c=pc("a"), label=fig_labels[i_fig][0])
        ax.plot([], [], "--", c=pc("a"), label=fig_labels[i_fig][1])
        
        r = np.linspace(0, 1, 200)
        ax.plot(r, nfw_mass_frac(r, 10), ":", c="k",
                label=r"${\rm NFW,\,c_{\rm vir}=10}$")


        ax.set_xlim(1e-3, 1)
        ax.set_ylim(0, 1)
        ax.set_xscale("log")
        ax.set_xlabel(r"$m = M_{\rm vir}/M_{\rm peak}$")
        ax.set_ylabel(r"$N(<m)$")

        ax.legend(loc="lower right", fontsize=17)
    
        fig.savefig(path.join(OUT_DIR, "mass_frac_cdf_%s.png" % suffix))


def main():
    #plot_mass_loss()
    #ufd_frequency()
    
    #survival_time()
    #stitching_errors()
    #mass_loss()

    #mass_function()
    #mass_function_resolution()
    radius_cdf()
    #radius_cdf()
    #mass_frac_cdf()
    #stitching_errors()

if __name__ == "__main__": main()
