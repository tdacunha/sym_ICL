import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import symlib
import os.path as path
from colossus.cosmology import cosmology
from colossus.halo import mass_so
import scipy.signal as signal
import scipy.stats as stats

ratio_limit = 0.1
suite = "SymphonyMilkyWay"
n_hr_max = 5
suite_lr = "SymphonyMilkyWayLR"
suite_hr = "SymphonyMilkyWayHR"
base_dir = "/sdf/home/p/phil1/ZoomIns"
plot_dir_fmt = "/scratch/phil1/plots/tracking/mah/h%d"

def plot_percentiles(ax, color, x, y, label=None):
    bins, n_min = 20, 100
    n, edges, _ = stats.binned_statistic(np.log10(x), y, "count", bins=bins)
    hi, _, _ = stats.binned_statistic(
        np.log10(x), y, lambda xx: np.percentile(xx, 50+68/2), bins=bins)
    mi, _, _ = stats.binned_statistic(
        np.log10(x), y, lambda xx: np.percentile(xx, 50), bins=bins)
    lo, _, _ = stats.binned_statistic(
        np.log10(x), y, lambda xx: np.percentile(xx, 50-68/2), bins=bins)
    ok = n > n_min

    mid = 10**((edges[1:] + edges[:-1])/2)

    ax.plot(mid[ok], mi[ok], c=color, label=label)
    #ax.plot(mid[ok], hi[ok], "--", c=color, lw=1.5)
    #ax.plot(mid[ok], lo[ok], "--", c=color, lw=1.5)

def mass_loss_vars():
    n_method = 2
    method_colors = [pc("r"), pc("b")]

    dm_dt = [[] for _ in range(n_method)]
    z_sub = [[] for _ in range(n_method)]
    msub_mhost = [[] for _ in range(n_method)]
    msub_mpeak = [[] for _ in range(n_method)]
    csub_chost = [[] for _ in range(n_method)]
    rsub_rhost = [[] for _ in range(n_method)]

    for i_host in range(symlib.n_hosts(suite)):
        if "HR" in suite:
            cut_mult = 8
        else:
            cut_mult = 1
        print(i_host, symlib.n_hosts(suite))
        plot_dir = plot_dir_fmt % i_host
        host_dir = symlib.get_host_directory(base_dir, suite, i_host)
        param = symlib.simulation_parameters(host_dir)
        cosmo = cosmology.setCosmology('', symlib.colossus_parameters(param))
        h, hist = symlib.read_subhalos(host_dir)

        c = symlib.read_cores(host_dir)
        scale = symlib.scale_factors(host_dir)
        z = 1/scale - 1
        t = cosmo.age(z)
        T = mass_so.dynamicalTime(z, "vir", "crossing")
        
        mp = param["mp"]/param["h100"]
        ok = hist["mpeak"]/mp > 300*cut_mult
        c, h, hist = c[ok], h[ok], hist[ok]
        
        smooth_order = 4
        smooth_window = 15

        for i in range(1, len(h)):
            R_host, m_host, ok_host = h["rvir"][0], h["mvir"][0], h["ok"][0]
            m_rs, m_c = h["mvir"][i], c["m_bound"][i]
            ok1_rs, ok2_rs, ok_c = h["ok"][i], c["ok_rs"][i], c["ok"][i]
            ok1_rs, ok2_rs = ok1_rs & h["ok"][0,:], ok2_rs & h["ok"][0,:]
            r_rs = np.sqrt(np.sum(h["x"][i]**2, axis=1))
            r_c = np.sqrt(np.sum(c["x"][i]**2, axis=1))

            snap = np.arange(h.shape[1], dtype=int)
            i_start = hist["merger_snap"][i]
            is_sub_rs = ok2_rs & (snap >= i_start)

            is_sub_method = [
                ok2_rs & (snap >= i_start),
                ok_c,
            ]
            m_method = [m_rs, m_c]
            r_method = [r_rs, r_c]

            for i_method in range(n_method):
                r, m = r_method[i_method], m_method[i_method]
                is_sub = is_sub_method[i_method]

                if (np.sum(is_sub) < smooth_window or
                    hist["merger_ratio"][i] > ratio_limit): continue
            
                dlnm_smooth = signal.savgol_filter(
                    np.log(m[is_sub]), smooth_window, smooth_order,
                    deriv=1, delta=1.0, mode="interp")[:-1]
                dt = t[is_sub][1:] - t[is_sub][:-1]
                dm_dt_i = dlnm_smooth/dt * T[is_sub][:-1]

                z_sub_i = z[is_sub][:-1]
                msub_mhost_i = (m[is_sub]/m_host[is_sub])[:-1]
                msub_mpeak_i = (m[is_sub]/np.max(m_rs[:i_start]))[:-1]
                csub_chost_i = (h[i,i_start]["cvir"]/h[0,is_sub]["cvir"])[:-1]
                rsub_rhost_i = (r[is_sub]/h[0,is_sub]["rvir"])[:-1]

                dm_dt[i_method].append(dm_dt_i)
                z_sub[i_method].append(z_sub_i)
                msub_mhost[i_method].append(msub_mhost_i)
                msub_mpeak[i_method].append(msub_mpeak_i)
                csub_chost[i_method].append(csub_chost_i)
                rsub_rhost[i_method].append(rsub_rhost_i)

    for i in range(n_method):
        dm_dt[i] = np.hstack(dm_dt[i])
        z_sub[i] = np.hstack(z_sub[i])
        msub_mhost[i] = np.hstack(msub_mhost[i])
        msub_mpeak[i] = np.hstack(msub_mpeak[i])
        csub_chost[i] = np.hstack(csub_chost[i])
        rsub_rhost[i] = np.hstack(rsub_rhost[i])

        ok = (msub_mhost[i] < 0.1) & (msub_mpeak[i] < 1.0) #& (z_sub[i] < 6)
        print(np.sum(ok))
        dm_dt[i] = dm_dt[i][ok]
        z_sub[i] = z_sub[i][ok]
        msub_mhost[i] = msub_mhost[i][ok]
        msub_mpeak[i] = msub_mpeak[i][ok]
        csub_chost[i] = csub_chost[i][ok]
        rsub_rhost[i] = rsub_rhost[i][ok]

    fig, axs = plt.subplots(2, 3, figsize=(24, 16))
    ax_z, ax_msub_mhost, ax_msub_mpeak = axs[0,0], axs[0,1], axs[0,2]
    ax_csub_chost, ax_rsub_rhost = axs[1,0], axs[1,1]
    fig.delaxes(axs[1,2])

    for i in range(n_method):
        c = method_colors[i]
        plot_percentiles(ax_z, c, z_sub[i], dm_dt[i])
        plot_percentiles(ax_msub_mhost, c, msub_mhost[i], dm_dt[i])
        plot_percentiles(ax_msub_mpeak, c, msub_mpeak[i], dm_dt[i])
        
        plot_percentiles(ax_csub_chost, c, csub_chost[i], dm_dt[i])
        plot_percentiles(ax_rsub_rhost, c, rsub_rhost[i], dm_dt[i])

    axs[0,0].set_ylabel(r"$\dot{m}_{\rm sub}/m_{\rm sub}$")
    axs[0,1].set_ylabel(r"$\dot{m}_{\rm sub}/m_{\rm sub}$")

    ax_z.set_xlabel(r"$1+z$")
    ax_msub_mhost.set_xlabel(r"$m_{\rm sub}/M_{\rm host}$")
    ax_msub_mpeak.set_xlabel(r"$m_{\rm sub}/m_{\rm peak}$")
    ax_csub_chost.set_xlabel(r"$c_{\rm sub,infall}/c_{\rm host}$")
    ax_rsub_rhost.set_xlabel(r"$r_{\rm sub}/R_{\rm vir,host}$")
    ax_msub_mpeak.grid()

    ax_z.set_xscale("log")
    ax_msub_mhost.set_xscale("log")
    ax_msub_mpeak.set_xscale("log")
    ax_csub_chost.set_xscale("log")
    ax_rsub_rhost.set_xscale("log")

    ax_z.set_xscale("log")
    ax_msub_mhost.set_ylim((-3, 0))
    ax_msub_mpeak.set_ylim((-3, 0))
    ax_csub_chost.set_ylim((-3, 0))
    ax_rsub_rhost.set_ylim((-3, 0))

    fig.savefig("../plots/core_tracking/mass_loss_vars_%s.pdf" % suite)

def mass_loss_msub_mpeak():
    n_method = 4
    colors = [pc("r"), pc("b"), pc("o"), pc("b")]
    panel = [0, 0, 1, 1]
    labels = [r"${\rm Rockstar}$", r"${\rm Particle{-}tracking}$",
              r"{\rm Fiducial}", r"${\rm High{-}res}$"]
    dm_dt = [[] for _ in range(n_method)]
    msub_mpeak = [[] for _ in range(n_method)]
    msub_mhost = [[] for _ in range(n_method)]

    for i_host in range(symlib.n_hosts(suite)):
        print(i_host, "out of", symlib.n_hosts(suite))
        plot_dir = plot_dir_fmt % i_host
        host_dir = symlib.get_host_directory(base_dir, suite, i_host)
        param = symlib.simulation_parameters(host_dir)
        cosmo = cosmology.setCosmology('', symlib.colossus_parameters(param))
        h, hist = symlib.read_subhalos(host_dir)

        c = symlib.read_cores(host_dir)
        scale = symlib.scale_factors(host_dir)
        z = 1/scale - 1
        t = cosmo.age(z)
        T = mass_so.dynamicalTime(z, "vir", "crossing")
        
        mp = param["mp"]/param["h100"]
        ok = hist["mpeak"]/mp > 300*8
        c, h, hist = c[ok], h[ok], hist[ok]

        if i_host < n_hr_max:
            host_dir = symlib.get_host_directory(base_dir, suite_hr, i_host)
            param = symlib.simulation_parameters(host_dir)
            cosmo = cosmology.setCosmology(
                '', symlib.colossus_parameters(param))
            h_hr, hist_hr = symlib.read_subhalos(host_dir)
            
            c_hr = symlib.read_cores(host_dir)
            scale = symlib.scale_factors(host_dir)

            mp = param["mp"]/param["h100"]
            ok = hist_hr["mpeak"]/mp > 300*8
            c_hr, h_hr, hist_hr = c_hr[ok], h_hr[ok], hist_hr[ok]

        smooth_order = 4
        smooth_window = 15

        for i in range(1, len(h)):
            R_host, m_host, ok_host = h["rvir"][0], h["mvir"][0], h["ok"][0]
            m_rs, m_c = h["mvir"][i], c["m_bound"][i]
            ok1_rs, ok2_rs, ok_c = h["ok"][i], c["ok_rs"][i], c["ok"][i]
            ok1_rs, ok2_rs = ok1_rs & h["ok"][0,:], ok2_rs & h["ok"][0,:]

            if i_host < n_hr_max:
                ok_c_hr = c_hr["ok"][i]
                m_hr = c_hr["m_bound"][i]
            else:
                ok_c_hr = None
                m_hr = None

            snap = np.arange(h.shape[1], dtype=int)
            i_start = hist["merger_snap"][i]

            is_sub_method = [
                ok2_rs & (snap >= i_start),
                ok_c, ok_c, ok_c_hr
            ]
            m_method = [m_rs, m_c, m_c, m_hr]

            for i_method in range(n_method):
                if ((i_method == 3 or i_method == 2) and
                    i_host >= n_hr_max): continue

                m = m_method[i_method]
                is_sub = is_sub_method[i_method]

                if (np.sum(is_sub) < smooth_window or
                    hist["merger_ratio"][i] > ratio_limit): continue
            
                dlnm_smooth = signal.savgol_filter(
                    np.log(m[is_sub]), smooth_window, smooth_order,
                    deriv=1, delta=1.0, mode="interp")[:-1]
                dlnm_smooth = signal.savgol_filter(
                    np.log(m[is_sub]), smooth_window, smooth_order,
                    deriv=1, delta=1.0)[:-1]
                dt = t[is_sub][1:] - t[is_sub][:-1]
                dm_dt_i = dlnm_smooth/dt * T[is_sub][:-1]
                msub_mpeak_i = (m[is_sub]/np.max(m_rs[:i_start]))[:-1]
                msub_mhost_i = (m[is_sub]/m_host[is_sub])[:-1]

                dm_dt[i_method].append(dm_dt_i)
                msub_mhost[i_method].append(msub_mhost_i)
                msub_mpeak[i_method].append(msub_mpeak_i)

        #break

    for i in range(n_method):
        dm_dt[i] = np.hstack(dm_dt[i])
        msub_mpeak[i] = np.hstack(msub_mpeak[i])
        msub_mhost[i] = np.hstack(msub_mhost[i])

        ok = (msub_mhost[i] < 0.1) & (msub_mpeak[i] < 1.0)
        dm_dt[i] = dm_dt[i][ok]
        msub_mpeak[i] = msub_mpeak[i][ok]
        msub_mhost[i] = msub_mhost[i][ok]

        print(len(msub_mpeak[i]))

    fig, ax = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

    plot_percentiles(ax[0], pc("r"), msub_mpeak[0], dm_dt[0], labels[0])
    plot_percentiles(ax[0], pc("b"), msub_mpeak[1], dm_dt[1], labels[1])
    plot_percentiles(ax[1], pc("r"), msub_mpeak[2], dm_dt[2], labels[2])
    plot_percentiles(ax[1], pc("b"), msub_mpeak[3], dm_dt[3], labels[3])
    
    ax[0].set_ylabel(r"$\dot{m}_{\rm sub}/(m_{\rm sub}/\tau_{\rm crossing})$")
    ax[0].set_xscale("log")
    ax[1].set_xscale("log")
    ax[0].set_xlabel(r"$m_{\rm sub}/m_{\rm peak}$")
    ax[1].set_xlabel(r"$m_{\rm sub}/m_{\rm peak}$")

    ax[0].set_xlim(1e-3, 1)
    ax[1].set_xlim(1e-3, 1)
    ax[0].set_ylim(-3, 0)
    ax[1].set_ylim(-3, 0)

    ax[0].set_ylim((-3, 0))
    ax[1].set_ylim((-3, 0))

    ax[0].legend(loc="lower right", fontsize=17)
    ax[1].legend(loc="lower right", fontsize=17)

    fig.savefig("../plots/core_tracking/mass_loss_msub_mpeak.pdf")


if __name__ == "__main__":
    palette.configure(True)
    #mass_loss_vars()
    mass_loss_msub_mpeak()
