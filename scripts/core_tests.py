import symlib
import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import os.path as path
from colossus.cosmology import cosmology
import read_mw_sats

invalid_hosts = [6, 9, 10, 16, 17, 31, 36, 37, 40, 42, 43]

def plot_mass_loss():
    palette.configure(False)

    out_dir = "../plots/core_tests/mass_loss"
    base_dir = "../tmp_data"
    suite = "SymphonyMilkyWay"
    sim_dir = symlib.get_host_directory(base_dir, suite, 0)

    param = symlib.simulation_parameters(suite)
    a = symlib.scale_factors(sim_dir)
    
    h, hist = symlib.read_subhalos(param, sim_dir)
    h = symlib.set_units_halos(h, a, param)
    c = symlib.read_cores(sim_dir)

    targets = np.arange(1, len(h), dtype=int)

    fig, ax = plt.subplots(2, figsize=(7, 14), sharex=True)
    
    for i_sub in targets:
        print(i_sub)
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

def mass_function():
    palette.configure(False)

    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    suite = "SymphonyMilkyWay"

    param = symlib.simulation_parameters(suite)
    
    ratios_rs, ratios_core = [], []
    r_rs, r_core = [], []
    n_host = 0

    for i_halo in range(symlib.n_hosts(suite)):
        if i_halo in invalid_hosts:
            continue
        print(i_halo)

        sim_dir = symlib.get_host_directory(base_dir, suite, i_halo)

        a = symlib.scale_factors(sim_dir)
        h, hist = symlib.read_subhalos(param, sim_dir)
        h = symlib.set_units_halos(h, a, param)
        c = symlib.read_cores(sim_dir)

        mp = param["mp"]/param["h100"]

        m_vir = h["mvir"][1:,-1]
        m_peak = np.max(h["mvir"], axis=1)[1:]
        okh = h["ok"][1:,-1] & (m_vir > 32*mp)
        
        r_h = np.sqrt(np.sum(h["x"][1:,-1]**2, axis=1))
        r_c = np.sqrt(np.sum(c["x"][1:,-1]**2, axis=1))
        r_half = c["r50_bound"][1:,-1]
        m_bound = c["m_bound"][1:,-1]
        m_tidal = c["m_tidal"][1:,-1]
        okc = c["ok"][1:,-1] & (m_bound > 32*mp) & (r_c > r_half)
        m0, r0 = h["mvir"][0,-1], h["rvir"][0,-1]
        
        ratios_rs.append(m_peak[okh]/m0)
        ratios_core.append(m_peak[okc]/m0)
        r_rs.append(r_h[okh]/r0)
        r_core.append(r_c[okc]/r0)
        n_host += 1

    ratios_rs = np.hstack(ratios_rs)
    ratios_core = np.hstack(ratios_core)
    r_rs = np.hstack(r_rs)
    r_core = np.hstack(r_core)

    r0_rs, r0_core = r_rs < 1, r_core < 1
    r1_rs, r1_core = r_rs < 0.25, r_core < 0.25
    r2_rs, r2_core = r_rs < 0.1, r_core < 0.1

    fig, ax = plt.subplots()
    ax.plot(np.sort(ratios_rs[r0_rs]),
            np.arange(np.sum(r0_rs))[::-1] / n_host,
            c=pc("r"))
    ax.plot(np.sort(ratios_core[r0_core]),
            np.arange(np.sum(r0_core))[::-1] / n_host, "--",
            c=pc("r"))
    ax.plot(np.sort(ratios_rs[r1_rs]),
            np.arange(np.sum(r1_rs))[::-1] / n_host,
            c=pc("o"))
    ax.plot(np.sort(ratios_core[r1_core]),
            np.arange(np.sum(r1_core))[::-1] / n_host, "--",
            c=pc("o"))
    ax.plot(np.sort(ratios_rs[r2_rs]),
            np.arange(np.sum(r2_rs))[::-1] / n_host,
            c=pc("b"))
    ax.plot(np.sort(ratios_core[r2_core]),
            np.arange(np.sum(r2_core))[::-1] / n_host, "--",
            c=pc("b"))
    
    ax.plot([], [], c=pc("r"), label=r"${\rm <R_{\rm vir}}$")
    ax.plot([], [], c=pc("o"), label=r"${\rm <R_{\rm vir}/4}$")
    ax.plot([], [], c=pc("b"), label=r"${\rm <R_{\rm vir}/10}$")
    ax.plot([], [], c=pc("k"), label=r"${\rm Rockstar}$")
    ax.plot([], [], "--", c=pc("k"), label=r"${\rm core-tracking}$")

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(r"$m=M_{\rm peak,sub}/M_{\rm host}$")
    ax.set_ylabel(r"$N(>m)$")

    ax.legend(loc="upper right", fontsize=17)
    
    fig.savefig("../plots/core_plots/core_shmf.png")

def t_survive(cosmo, a, infall_snap, ok):
    z = 1/a - 1
    t = cosmo.age(z)

    t_first = t[infall_snap]
    #print(infall_snap, np.arange(len(ok), dtype=int)[ok])
    t_last = np.max(t[ok])


    return t_last - t_first

def survival_time():
    palette.configure(False)

    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    suite = "SymphonyMilkyWay"

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
        h, hist = symlib.read_subhalos(param, sim_dir)
        h = symlib.set_units_halos(h, a, param)
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

    print(n_merger)
    print(n_disrupt)

def stitching_errors():
    palette.configure(False)

    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    suite = "SymphonyMilkyWay"

    param = symlib.simulation_parameters(suite)
    
    n_host = 0

    n_bins = 200
    r_bins = np.logspace(-3, 0, n_bins+1)

    N_tot_now = np.zeros(n_bins)
    N_err_now = np.zeros(n_bins)

    N_tot_disrupt = np.zeros(n_bins)
    N_err_disrupt = np.zeros(n_bins)

    for i_halo in range(symlib.n_hosts(suite)):
        if i_halo in invalid_hosts:
            continue
        print(i_halo)

        sim_dir = symlib.get_host_directory(base_dir, suite, i_halo)

        a = symlib.scale_factors(sim_dir)
        h, hist = symlib.read_subhalos(param, sim_dir)
        h = symlib.set_units_halos(h, a, param)
        c = symlib.read_cores(sim_dir)

        mp = param["mp"]/param["h100"]

        host = h[0]
        h, c, hist = h[1:], c[1:], hist[1:]

        m_vir, m_peak = h["mvir"][:,-1], hist["mpeak"]

        dr_hc = np.sqrt(np.sum((h["x"][:,-1] - c["x"][:,-1])**2, axis=1))
        r_h = np.sqrt(np.sum(h["x"][:,-1]**2, axis=1))
        r_c = np.sqrt(np.sum(c["x"][:,-1]**2, axis=1))
        r_half = c["r50_bound"][:,-1]
        m_bound = c["m_bound"][:,-1]
        m_tidal = c["m_tidal"][:,-1]

        ok = (h["ok"][:,-1] & (m_vir > 32*mp) &
              c["ok"][:,-1] & (m_bound > 32*mp) &
              (r_c > r_half))
        is_err = dr_hc > r_half

        r0 = host["rvir"][-1]

        n_tot, _ = np.histogram(r_h[ok]/r0, bins=r_bins)
        n_err, _ = np.histogram(r_h[ok & is_err]/r0, bins=r_bins)
        
        N_tot_now += n_tot
        N_err_now += n_err

        r_peri = np.min(np.sum(h["x"]**2, axis=2), axis=1)
        r_c = np.sqrt(np.sum(c["x"]**2, axis=2))
        r_h = np.sqrt(np.sum(h["x"]**2, axis=2))
        snap = np.arange(h.shape[1])

        ok = (h["ok"] & (h["mvir"] > 32*mp) & 
              c["ok"] & (c["m_bound"] > 32*mp) & 
              (r_c > c["r50_bound"]))

        last_snap = np.zeros(len(h), dtype=int)
        rd_hc = np.zeros(len(h))
        r_half = np.zeros(len(h))
        rh_last = np.zeros(len(h))
        for i in range(len(last_snap)):
            if np.sum(ok[i]) == 0 or h["ok"][i,-1]:
                last_snap[i] = -1
            else:
                last_snap[i] = np.max(snap[ok[i]])
            dr_hc[i] = np.sqrt(np.sum((
                h["x"][i,last_snap[i],:] - c["x"][i,last_snap[i],:])**2,
            ))
            r_half = c["r50_bound"][i,last_snap[i]]
            rh_last = r_h[i,last_snap[i]]

        is_valid = last_snap != -1
        is_err = r_half < dr_hc
        r0 = host["rvir"][last_snap]

        n_tot_disrupt, _ = np.histogram(
            r_peri[is_valid]/r0[is_valid], bins=r_bins)
        n_err_disrupt, _ = np.histogram(
            r_peri[is_valid & is_err]/r0[is_valid & is_err], bins=r_bins)
        N_tot_disrupt += n_tot_disrupt
        N_err_disrupt += n_err_disrupt

    N_tot_now = np.cumsum(N_tot_now)
    N_err_now = np.cumsum(N_err_now)
    N_tot_disrupt = np.cumsum(N_tot_disrupt)
    N_err_disrupt = np.cumsum(N_err_disrupt)

    print("now:    ", N_tot_now[-1], N_err_now[-1])
    print("disrupt:", N_tot_disrupt[-1], N_err_disrupt[-1])

    ok_now = N_tot_now > 2
    ok_disrupt = N_tot_disrupt > 10

    f_err_now = N_err_now[ok_now]/N_tot_now[ok_now]
    right_bin_now = r_bins[1:][ok_now]
    f_err_disrupt = N_err_disrupt[ok_disrupt]/N_tot_disrupt[ok_disrupt]
    right_bin_disrupt = r_bins[1:][ok_disrupt]

    fig, ax = plt.subplots()
    ax.plot(right_bin_now, f_err_now, c=pc("r"))
    ax.set_xscale("log")
    ax.set_xlabel(r"$r/R_{\rm vir}$")
    ax.set_ylabel(r"$f_{\rm error}(< r/R_{\rm vir})$")
    fig.savefig("../plots/core_plots/core_stitching_now.png")

    fig, ax = plt.subplots()
    ax.plot(right_bin_disrupt, f_err_disrupt, c=pc("r"))
    ax.set_xscale("log")
    ax.set_xlabel(r"$r_{\rm peri}/R_{\rm vir}$")
    ax.set_ylabel(r"$f_{\rm error}(< r_{\rm peri}/R_{\rm vir})$")
    fig.savefig("../plots/core_plots/core_stitching_disrupt.png")

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

def ufd_frequency():
    palette.configure(False)

    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    suite = "SymphonyMilkyWay"

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
            mpeak_ok = hist["mpeak"][1:] > mpeak_cuts[j]
        
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
    fig.savefig("../plots/core_plots/ufd_freq_%s.png" % names[i])

def nfw_mass_frac(r_rvir, cvir):
    r_rs = r_rvir * cvir
    def M_enc(x): return np.log(1 + x) + 1/(1 + x) - 1
    return M_enc(r_rs) / M_enc(cvir)


def radial_distribution():
    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    suite = "SymphonyMilkyWay"
    
    bins = [1e8, 1e9, 1e10, 1e11]
    colors = [pc("r"), pc("o"), pc("b")]

    r_bins = np.linspace(0, 1, 101)
    h_n = np.zeros((len(bins)-1, len(r_bins)-1))
    c_n = np.zeros((len(bins)-1, len(r_bins)-1))

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

        for i in range(len(bins) - 1):
            mpeak_ok = ((hist["mpeak"] > bins[i]) &
                        (hist["mpeak"] < bins[i+1]))[1:]
            h_n[i,:] += np.histogram(h_r[h_ok & mpeak_ok]/h["rvir"][0,-1],
                                     bins=r_bins)[0]
            c_n[i,:] += np.histogram(c_r[c_ok & mpeak_ok]/h["rvir"][0,-1],
                                     bins=r_bins)[0]
    
    palette.configure(True)

    fig, ax = plt.subplots()
    r_mid = (r_bins[1:] + r_bins[:-1]) / 2

    for i in range(len(bins) - 1):
        h_n[i] = np.cumsum(h_n[i]) / np.sum(h_n[i])
        c_n[i] = np.cumsum(c_n[i]) / np.sum(c_n[i])
        ax.plot(r_mid, h_n[i,:], "--", c=colors[i])
        ax.plot(r_mid, c_n[i,:], "-", c=colors[i])

    ax.plot(r_mid, nfw_mass_frac(r_mid, 10), ":", c="k",
            label=r"${\rm NFW,\,c_{\rm vir}=10}$")
    ax.plot([], [], c=colors[0],
            label=r"$10^8 < M_{\rm peak}/M_\odot < 10^9$")
    ax.plot([], [], c=colors[1],
            label=r"$10^9 < M_{\rm peak}/M_\odot < 10^{10}$")
    ax.plot([], [], c=colors[2],
            label=r"$10^{10} < M_{\rm peak}/M_\odot < 10^{11}$")

    ax.set_xlabel(r"$r/R_{\rm vir}$")
    ax.set_ylabel(r"$N(<r/R_{\rm vir})/N(<R_{\rm vir})$")
    ax.legend(loc="lower right", fontsize=17)

    fig.savefig("../plots/core_plots/radial_cdf.png")

def main():
    #plot_mass_loss()
    #mass_function()
    #survival_time()
    #stitching_errors()
    #ufd_frequency()
    radial_distribution()

if __name__ == "__main__": main()
