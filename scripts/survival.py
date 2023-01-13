import numpy as np
import scipy.interpolate as interpolate
import symlib
from colossus.cosmology import cosmology
from colossus.halo import mass_so
import palette
from palette import pc
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

max_ratio = 0.1
base_dir = "/sdf/home/p/phil1/ZoomIns"

def kaplan_meier(dt, censored, t_eval, decreasing=False):
    if decreasing:
        dt, t_eval = -dt, -t_eval
    dt_death = dt[~censored]
    dt = np.sort(dt)
    ti, di = np.unique(dt_death, return_counts=True)
    ni = len(dt) - np.searchsorted(dt, ti, side="left")

    hi = 1 - di/ni
    Si = np.cumprod(hi)

    # Greenwood forumla
    diff = ni-di
    diff[(ni-di) <= 0] = 1
    err = Si * np.sqrt(np.cumsum(di / (ni*diff)))
    
    f_S = interpolate.interp1d(ti, Si, fill_value=np.nan,
                               bounds_error=False)
    f_err = interpolate.interp1d(ti, err, fill_value=np.nan,
                                 bounds_error=False)

    return f_S(t_eval), f_err(t_eval)
    
def calc_mpeak_pre(h, hist):
    mpeak_pre = np.zeros(len(h))
    for i in range(1, len(h)):
        mpeak_pre[i] = np.max(h["mvir"][i,:hist["first_infall_snap"][i]+1])
    mpeak_pre[0] = np.max(h["mvir"][0,:])
    return mpeak_pre

def survival_time():
    suite = "SymphonyMilkyWay"
    suite_lr = "SymphonyMilkyWayLR"
    suite_hr = "SymphonyMilkyWayHR"

    npeak = []
    t_c = []
    t_rs = []
    m_c = []
    m_rs = []
    censored_c = []
    censored_rs = []

    npeak_lr = []
    t_lr = []
    m_lr = []
    censored_lr = []

    npeak_hr = []
    t_hr = []
    m_hr = []
    censored_hr = []

    for suite_i in [suite, suite_lr, suite_hr]:

        n_hosts = symlib.n_hosts(suite_i)

        for i_host in range(n_hosts):
            print(suite_i, i_host)
            sim_dir = symlib.get_host_directory(base_dir, suite_i, i_host)
        
            param = symlib.simulation_parameters(sim_dir)
            mp = param["mp"]/param["h100"]
            cosmo = cosmology.setCosmology('', symlib.colossus_parameters(param))
            a = symlib.scale_factors(sim_dir)
            z = 1/a - 1
            t = cosmo.age(z)
            T_orbit = mass_so.dynamicalTime(z, "vir", "crossing")

            h, hist = symlib.read_subhalos(sim_dir)
            c = symlib.read_cores(sim_dir)
            mpeak_pre = calc_mpeak_pre(h, hist)

            last_snap_rs = np.zeros(len(h), dtype=int)
            last_snap_c = np.zeros(len(h), dtype=int)
            snap = np.arange(h.shape[1], dtype=int)
            for i in range(1, len(h)):
                ok_rs, ok_c = c["ok_rs"][i], c["ok"][i]
                if np.sum(ok_c) == 0:
                    last_snap_c[i] = hist["first_infall_snap"][i]
                else:
                    last_snap_c[i] = np.max(snap[ok_c])
                if np.sum(ok_rs) == 0:
                    last_snap_rs[i] = hist["first_infall_snap"][i]
                else:
                    last_snap_rs[i] = np.max(snap[ok_rs])
            last_snap_c[0], last_snap_rs[0] = len(snap)-1, len(snap)-1
            dt_rs = t[last_snap_rs] - t[hist["first_infall_snap"]]
            dt_c = t[last_snap_c] - t[hist["first_infall_snap"]]
            dt_rs /= T_orbit[hist["first_infall_snap"]]
            dt_c /= T_orbit[hist["first_infall_snap"]]
            
            idx = np.arange(len(c), dtype=int)
            dm_c = c["m_bound"][idx,last_snap_c]/mpeak_pre
            dm_rs = h["mvir"][idx,last_snap_rs]/mpeak_pre
    
            ratio_ok = ((hist["merger_ratio"] < max_ratio) &
                        (hist["first_infall_snap"] < last_snap_c) &
                        (hist["first_infall_snap"] < last_snap_rs))
            ratio_ok[0] = False

            if suite_i == suite:
                censored_rs.append(c["ok_rs"][ratio_ok,-1])
                censored_c.append(c["ok"][ratio_ok,-1])
                t_rs.append(dt_rs[ratio_ok])
                t_c.append(dt_c[ratio_ok])
                m_rs.append(dm_rs[ratio_ok])
                m_c.append(dm_c[ratio_ok])
                npeak.append((mpeak_pre/mp)[ratio_ok])
            elif suite_i == suite_lr:
                censored_lr.append(c["ok"][ratio_ok,-1])
                t_lr.append(dt_c[ratio_ok])
                m_lr.append(dm_c[ratio_ok])
                npeak_lr.append((mpeak_pre/mp)[ratio_ok])
            elif suite_i == suite_hr:
                censored_hr.append(c["ok"][ratio_ok,-1])
                t_hr.append(dt_c[ratio_ok])
                m_hr.append(dm_c[ratio_ok])
                npeak_hr.append((mpeak_pre/mp)[ratio_ok])
                
    censored_rs = np.hstack(censored_rs)
    censored_c = np.hstack(censored_c)
    t_rs = np.hstack(t_rs)
    t_c = np.hstack(t_c)
    m_rs = np.hstack(m_rs)
    m_c = np.hstack(m_c)
    npeak = np.hstack(npeak)

    censored_lr = np.hstack(censored_lr)
    t_lr = np.hstack(t_lr)
    m_lr = np.hstack(m_lr)
    npeak_lr = np.hstack(npeak_lr)

    censored_hr = np.hstack(censored_hr)
    t_hr = np.hstack(t_hr)
    m_hr = np.hstack(m_hr)
    npeak_hr = np.hstack(npeak_hr)

    n_eval = 100
    t_eval = 10**np.linspace(np.log10(0.1), np.log10(50), n_eval)
    m_eval = 10**np.linspace(-4, 1, n_eval)

    lims = [(10**2.5, 10**3), (10**3.5, 10**4), (10**4.5, 10**5)]
    labels = [r"$10^{2.5}<N_{\rm peak}<10^3$",
              r"$10^{3.5}<N_{\rm peak}<10^4$",
              r"$10^{4.5}<N_{\rm peak}<10^5$"]
    colors = [pc("r"), pc("o"), pc("b")]

    fig1, (ax_m, ax_hr) = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

    for i in range(len(lims)):
        low, high = lims[i]
        ok = (npeak > low) & (npeak < high)
        med_npeak = np.median(npeak[ok])

        S_t_c, err_t_c = kaplan_meier(t_c[ok], censored_c[ok], t_eval)
        S_t_rs, err_t_rs = kaplan_meier(t_rs[ok], censored_rs[ok], t_eval)
        S_m_c, err_m_c = kaplan_meier(m_c[ok], censored_c[ok], m_eval,
                                      decreasing=True)
        S_m_rs, err_m_rs = kaplan_meier(m_rs[ok], censored_rs[ok], m_eval,
                                        decreasing=True)
        
        ok_lr = (npeak_lr > low) & (npeak_lr < high)
        S_m_lr, err_m_lr = kaplan_meier(m_lr[ok_lr], censored_lr[ok_lr],
                                        m_eval, decreasing=True)
        ok_hr = (npeak_hr > low) & (npeak_hr < high)
        S_m_hr, err_m_hr = kaplan_meier(m_hr[ok_hr], censored_hr[ok_hr],
                                        m_eval, decreasing=True)

        #ok = ~np.isnan(S_t_rs)
        #ax_t.plot(t_eval[ok], S_t_rs[ok], "--", c=colors[i])
        #ax_t.fill_between(t_eval[ok], (S_t_rs+err_t_rs)[ok],
        #                  (S_t_rs-err_t_rs)[ok], color=colors[i], alpha=0.2)
        #ok = ~np.isnan(S_t_c)
        #ax_t.plot(t_eval[ok], S_t_c[ok], "-", c=colors[i])
        #ax_t.fill_between(t_eval[ok], (S_t_c+err_t_c)[ok],
        #                  (S_t_c-err_t_c)[ok], color=colors[i], alpha=0.2)
        ok = ~np.isnan(S_m_rs)
        ax_m.plot(m_eval[ok], S_m_rs[ok], "--", c=colors[i])
        ax_m.fill_between(m_eval[ok], (S_m_rs+err_m_rs)[ok],
                          (S_m_rs-err_m_rs)[ok], color=colors[i], alpha=0.2)

        ok = ~np.isnan(S_m_c)
        ax_m.plot(m_eval[ok], S_m_c[ok], "-", c=colors[i])
        ax_m.fill_between(m_eval[ok], (S_m_c+err_m_c)[ok],
                          (S_m_c-err_m_c)[ok], color=colors[i], alpha=0.2)

        func_S_m_c = interpolate.interp1d(np.log10(m_eval), S_m_c)
        ref_n = 50
        ax_m.plot(ref_n/med_npeak, func_S_m_c(np.log10(ref_n/med_npeak)),
                  "o", c=colors[i], ms=12)

        ok = ~np.isnan(S_m_lr)
        ax_hr.plot(m_eval[ok], S_m_lr[ok], "-", c=colors[i])
        ax_hr.fill_between(m_eval[ok], (S_m_lr+err_m_lr)[ok],
                           (S_m_lr-err_m_lr)[ok], color=colors[i], alpha=0.2)

        ok = ~np.isnan(S_m_hr)
        ax_hr.plot(m_eval[ok], S_m_hr[ok], "--", c=colors[i])
        ax_hr.fill_between(m_eval[ok], (S_m_hr+err_m_hr)[ok],
                           (S_m_hr-err_m_hr)[ok], color=colors[i], alpha=0.2)

        ax_m.plot([], [], c=colors[i], label=labels[i])

    ax_m.plot([], [], "--", c="k", label=r"$\textsc{Rockstar}$")
    ax_hr.plot([], [], "-", c="k", label=r"${\rm fiducial\ resolution}$")
    ax_hr.plot([], [], "--", c="k", label=r"${\rm high\ resolution}$")

    #ax_t.set_ylim(0, 1.2)
    ax_m.legend(loc="upper left", fontsize=17)
    #ax_t.set_xlabel(r"$t_{\rm disrupt}/T_{\rm crossing,infall}$")
    #ax_t.set_ylabel(r"${\rm Pr}(>t_{\rm disrupt}/T_{\rm crossing,infall})$")
    #ax_t.set_xscale("log")
    #ax_t.set_xlim(0, 30)
    ax_m.set_ylim(0, 1.25)
    ax_m.set_xlabel(r"$M_{\rm disrupt}/M_{\rm peak}$")
    ax_m.set_ylabel(r"${\rm Pr}(<M_{\rm disrupt}/M_{\rm peak})$")
    ax_m.set_xscale("log")
    ax_m.set_xlim((1e-4, 1))

    #ax_hr.set_ylim(0, 1.2)
    ax_hr.legend(loc="upper left", fontsize=17)
    ax_hr.set_xlabel(r"$M_{\rm disrupt}/M_{\rm peak}$")
    #ax_hr.set_ylabel(r"${\rm Pr}(<M_{\rm disrupt}/M_{\rm peak})$")
    ax_hr.set_xscale("log")
    ax_hr.set_xlim((1e-4, 1))

    fig1.savefig("../plots/core_tracking/survival.pdf")
    #fig2.savefig("../plots/core_tracking/survival_hr.pdf")

def main():
    palette.configure(True)

    survival_time()

if __name__ == "__main__": main()
