import numpy as np
import symlib
import palette
from palette import pc
import matplotlib.pyplot as plt

param_sets = [
    (64, 32, pc("r"), "-", 2),
    (64, 16, pc("o"), "-", 2),
    (64, 8, pc("g"), "-", 2),
    (64, 4, pc("b"), "-", 2),
    (64, 1, pc("p"), "-", 2),

    (32, 32, pc("r"), "--", 2),
    (32, 16, pc("o"), "--", 2),
    (32, 8, pc("g"), "--", 2),
    (32, 4, pc("b"), "--", 2),
    (32, 1, pc("p"), "--", 2),

    (16, 32, pc("r"), ":", 2),
    (16, 32, pc("r"), ":", 2),
    (16, 16, pc("o"), ":", 2),
    (16, 8, pc("g"), ":", 2),
    (16, 4, pc("b"), ":", 2),
    (16, 1, pc("p"), ":", 2),

    (8, 8, pc("g"), ".-", 2),
    (8, 4, pc("b"), ".-", 2),
    (8, 1, pc("p"), ".-", 2)

]

base_dir = "/sdf/home/p/phil1/ZoomIns"
suite = "SymphonyMilkyWay"
plot_dir = "../plots/core_tracking/param_tests"
suffixes = ["k%02d_n%02d" % (p[0], p[1]) for p in param_sets]
suffixes[10] = "K16_N32"

def mpeak_pre(h, hist):
    m_pre = np.zeros(len(h))
    snap = np.arange(h.shape[1], dtype=int)
    m_pre[0] = hist["mpeak"][0]
    for i in range(1, len(h)):
        ok = np.where(h[i,:]["ok"] & (snap < hist["first_infall_snap"][i]))[0]
        if len(ok) == 0:
            m_pre[i] = 0
        else:
            m_pre[i] = np.max(h["mvir"][i,ok])
    return m_pre

def plot_valid_z0_fractions():    
    param = symlib.simulation_parameters(suite)
    mp = param["mp"]/param["h100"]

    n_bins = 200
    n_peak_bins = 10**np.linspace(np.log10(3e2), np.log10(1e5), n_bins+1)
    is_ok_rs = []
    is_err_rs = []
    n_peak = []
    r_rs = []
    f_core_rs = []
    
    is_ok_c = [[] for _ in range(len(suffixes))]
    r_c = [[] for _ in range(len(suffixes))]
    is_h_c_err = [[] for _ in range(len(suffixes))]
    f_core_c = [[] for _ in range(len(suffixes))]

    n_hosts = symlib.n_hosts(suite)
    for i_host in range(n_hosts):
        print("%d/%d" % (i_host+1, n_hosts))
        sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
        h, hist = symlib.read_subhalos(sim_dir)
        c = symlib.read_cores(sim_dir, suffix=suffixes[10])

        m_pre = mpeak_pre(h, hist)
        is_ok_rs_i = h["ok"][:,-1]
        n_peak_i = m_pre/mp
        r = np.sqrt(np.sum(h["x"][:,-1]**2, axis=1))
        r_rs_i = r/h["rvir"][0,-1]
        is_err_rs_i = c["f_core_rs"][:,-1] <= 0

        is_ok_rs_i = is_ok_rs_i & (r > c["r50_bound_rs"][:,-1])

        is_ok_rs.append(is_ok_rs_i[1:])
        is_err_rs.append(is_err_rs_i[1:])
        r_rs.append(r_rs_i[1:])
        n_peak.append(n_peak_i[1:])
        f_core_rs.append(c["f_core_rs"][1:,-1])

        for i_suf in range(len(suffixes)):
            if i_suf > 0:
                c = symlib.read_cores(sim_dir, suffix=suffixes[i_suf])
            r = np.sqrt(np.sum(c["x"][:,-1]**2, axis=1))
            is_ok_c_i = c["ok"][:,-1] & (c["f_core"][:,-1] > 0) & (r > c["r50_bound"][:,-1])
            dx_h_c = c["x"][:,-1] - h["x"][:,-1]
            dr_h_c = np.sqrt(np.sum(dx_h_c**2, axis=1))
            is_h_c_err_i = dr_h_c > c["r50_bound"][:,-1]

            is_ok_c[i_suf].append(is_ok_c_i[1:])
            is_h_c_err[i_suf].append(is_h_c_err_i[1:])
            f_core_c[i_suf].append(c["f_core"][1:,-1])

    is_ok_rs = np.hstack(is_ok_rs)
    is_err_rs = np.hstack(is_err_rs)
    r_rs = np.hstack(r_rs)
    n_peak = np.hstack(n_peak)
    f_core_rs = np.hstack(f_core_rs)

    for cut in [3e2, 1e3, 1e4, 1e5]:
        print("%1g particles" % cut)
        for i_suf in range(len(suffixes)):
            is_ok_c[i_suf] = np.hstack(is_ok_c[i_suf])
            is_h_c_err[i_suf] = np.hstack(is_h_c_err[i_suf])
            f_core_c[i_suf] = np.hstack(f_core_c[i_suf])

            n_peak_ok = n_peak >= cut
            true_ok_h = is_ok_rs & (~is_err_rs)
            both_ok = is_ok_c[i_suf] & true_ok_h & n_peak_ok
            either_ok = (is_ok_c[i_suf] | true_ok_h) & n_peak_ok
            rs_better = both_ok & is_h_c_err[i_suf] & (f_core_c[i_suf] < f_core_rs)
            c_better = both_ok & is_h_c_err[i_suf] & (f_core_c[i_suf] > f_core_rs)
            rs_missing = either_ok & (~true_ok_h)
            c_missing = either_ok & (~is_ok_c[i_suf])

            print("%s %d %d %3d %3d %3d %3d" % (
                suffixes[i_suf], np.sum(either_ok), np.sum(both_ok),
                np.sum(rs_better), np.sum(c_better), np.sum(rs_missing),
                np.sum(c_missing)
            ))
        print()

    lim_colors = [pc("r"), pc("o"), pc("b")]
    r_lims = [1, 0.5, 0.25]

    fig, ax = plt.subplots()
    for i_lim in range(len(r_lims)):
        r_lim = r_lims[i_lim]
        in_r = r_rs < r_lim

        N_tot, edges = np.histogram(n_peak[is_ok_rs & in_r], bins=n_peak_bins)
        N_tot = np.cumsum(N_tot[::-1])[::-1]
        N_tot += np.sum(n_peak[in_r & is_ok_rs] > n_peak_bins[-1])

        N_err, _ = np.histogram(n_peak[is_ok_rs & is_err_rs & in_r], bins=n_peak_bins)
        N_err = np.cumsum(N_err[::-1])[::-1]
        N_err += np.sum(n_peak[is_ok_rs & is_err_rs & in_r] > n_peak_bins[-1])
    
        N_ok = N_err > 2

        f_err = N_err / N_tot
        ax.plot(edges[:-1][N_ok], f_err[N_ok], c=lim_colors[i_lim],
                label=r"$<%.2f\cdot R_{\rm vir}$")
    ax.set_xscale("log")
    ax.set_xlim(n_peak_bins[0], n_peak_bins[-1])
    ax.set_xlabel(r"$N_{\rm peak,pre}$")
    ax.set_ylabel(r"$f_{\rm err,RS}(>N_{\rm peak,pre};z=0)$")
    ax.legend(loc="upper left")
    fig.savefig("%s/f_p_err_rs.png" % plot_dir)

    fig, ax = plt.subplots()

def main():
    palette.configure(False)
    plot_valid_z0_fractions()
    
 
if __name__ == "__main__": main()
