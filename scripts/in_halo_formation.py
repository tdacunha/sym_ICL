import symlib
import palette
from palette import pc
import scipy
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

def calc_mpeaks(idx, b, mvir):
    mpeak = np.zeros(len(idx))
    for i, j in enumerate(idx):
        mpeak[i] = np.max(mvir[b["start"][j]: b["end"][j]])
    return mpeak

def main():
    palette.configure(True)

    base_dir = "/sdf/home/p/phil1/ZoomIns"
    #suites = ["SymphonyLMC", "SymphonyMilkyWay", "SymphonyGroup",
    #          "SymphonyLCluster"]
    #suites = ["SymphonyCluster"]
    suites = ["SymphonyMilkyWay"]
    
    n_bins = 200
    bins = 10**np.linspace(0, 6, n_bins + 1)
    hists = [np.zeros(n_bins, dtype=int) for _ in range(6)] 

    fig_many, ax_many = plt.subplots()

    for suite in suites:
        print(suite)
        scale = symlib.scale_factors(suite)
        param = symlib.simulation_parameters(suite)
        n_host = symlib.n_hosts(suite)

        mp = param["mp"]/param["h100"]

        for i_host in range(n_host):
            print("%s %d/%d" % (suite, i_host, n_host-1))
            sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
            
            h, hist = symlib.read_subhalos(sim_dir)
            max_snap = len(h[0]) - 1

            b = symlib.read_branches(sim_dir)
            mvir, snap = symlib.read_tree(sim_dir, ["mvir", "snap"])

            is_sub = b["is_main_sub"]
            is_surv = snap[b["start"]] == max_snap
            is_ok = b["is_real"] & (~b["is_disappear"])

            flags = [
                is_sub,              is_sub & is_ok,
                is_sub & is_surv,    is_sub & is_surv & is_ok,
                is_sub & (~is_surv), is_sub & (~is_surv) & is_ok,
                is_sub & (~is_surv) & (~is_ok)
            ]

            idxs = [np.where(flag)[0] for flag in flags]
            npeaks = [calc_mpeaks(idx, b, mvir)/mp for idx in idxs]

            for i in range(len(hists)):
                hists[i] += np.cumsum(np.histogram(
                    npeaks[i], bins=bins)[0][::-1])[::-1]

            n_ok =  np.cumsum(np.histogram(
                npeaks[5], bins=bins)[0][::-1])[::-1]
            n_all =  np.cumsum(np.histogram(
                npeaks[4], bins=bins)[0][::-1])[::-1]
            ok = n_all > 10

            n4, n5 = np.sum(npeaks[4] > 1e4), np.sum(npeaks[5] > 1e4)
            frac_1e4 = 1 - n5/n4

            if frac_1e4 > 0.2:
                print("%s %d/%d" % (suite, i_host, n_host-1))
                print(sim_dir)
                print("%.3f" % frac_1e4)

                is_big = npeaks[6] > 1e4
                big_idx = idxs[6][is_big]

                h_cmov, hist_cmov = symlib.read_subhalos(sim_dir, comoving=True)

                print(h_cmov["mvir"][1][197:])
                print(h_cmov["x"][1][197:])
                print(h_cmov["v"][1][197:])

                print("branch_idx", big_idx)
                print("start", b["start"][big_idx])
                print("n_snap", b["end"][big_idx] - b["start"][big_idx])
                
                print("last_snap", snap[b["start"][big_idx]])
                print(np.max(snap))
                print()
                    

            ax_many.plot(bins[1:][ok], 1 - n_ok[ok]/n_all[ok],
                         lw=1, c=pc("b"))


    fig, ax = plt.subplots()
    colors = [pc("r"), pc("o"), pc("b")]
    labels = [r"${\rm All}$", r"${\rm Surviving},\ z=0$", r"${\rm Disrupted}$"]

    for i in range(len(colors)):
        ok = (hists[2*i] > 0) & (hists[2*i + 1] > 0)
        frac = 1 - hists[2*i + 1][ok]/hists[2*i][ok]
        ax.plot(bins[1:][ok], frac, c=colors[i], label=labels[i])

    ylo, yhi = ax.set_ylim(0, None)
    ax.set_xlim(1, 1e4)
    ax.legend(loc="upper right", fontsize=17)
    ax.set_xscale("log")
    ax.set_ylabel(r"$f_{\rm error,RS,leaf}(> N_{\rm peak})$")
    ax.set_xlabel(r"$N_{\rm peak}$")
    if len(suites) > 1:
        fig.savefig("../plots/core_tracking/in_halo_frac.pdf")
    else:
        fig.savefig("../plots/core_tracking/in_halo_frac_%s.pdf" % suites[0])

    ylo, yhi = ax_many.set_ylim(0, None)
    ax_many.set_xlim(1, 1e4)
    ax_many.set_xscale("log")
    ax_many.set_ylabel(r"$f_{\rm error,RS,leaf}(> N_{\rm peak})$")
    ax_many.set_xlabel(r"$N_{\rm peak}$")
    if len(suites) > 1:
        fig_many.savefig("../plots/core_tracking/in_halo_frac_all.pdf")
    else:
        fig_many.savefig("../plots/core_tracking/in_halo_frac_all_%s.pdf" % suites[0])


if __name__ == "__main__": main()
