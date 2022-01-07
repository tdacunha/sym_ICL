import numpy as np
import matplotlib.pyplot as plt
import lib
import os.path as path
import palette
from palette import pc
import scipy.interpolate as interpolate

OMEGA_M = 0.286
#DIR_FORMAT = "../tmp_data/Halo%03d"
#HALO_NUMS = [4]
DIR_FORMAT = "/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo%03d"
HALO_NUMS = [4, 113, 169, 170, 222, 229, 282, 327, 349, 407, 453, 523, 625,
             659, 666, 719, 747, 756, 788, 858, 953, 975, 983]
DIR_NAMES = [DIR_FORMAT % n for n in HALO_NUMS]
MP = 2.8e5
MVIR_CONV = MP * 300
INDIVIDUAL_SUBS = 10
colors = [pc("k")] + [pc("r", 0.2 + (0.8 - 0.2)*p/INDIVIDUAL_SUBS) 
                      for p in range(INDIVIDUAL_SUBS)] + [pc("r"), pc("b")]
N_BINS = 200

def main():
    palette.configure(False)

    r_subs, prog_idxs = [None]*len(DIR_NAMES), [None]*len(DIR_NAMES)
    for i in range(len(DIR_NAMES)):
        print(DIR_NAMES[i])
        r_subs[i], prog_idxs[i] = sub_info(DIR_NAMES[i])
    r_sub, prog_idx = lib.flatten(r_subs), lib.flatten(prog_idxs)
        
    n_tot = interpolate.interp1d(r_sub, np.arange(len(r_sub)))

    r_sub_pre = [None]*(INDIVIDUAL_SUBS + 1 + 2)
    r_sub_pre[0] = np.sort(r_sub[prog_idx > -1])
    for i in range(1, INDIVIDUAL_SUBS + 1 + 2):
        if i <= INDIVIDUAL_SUBS:
            r_sub_pre[i] = np.sort(r_sub[(prog_idx >= 1) & (prog_idx <= i)])
        else:
            r_sub_pre[i] = np.sort(r_sub[prog_idx == i])

    # Plotting

    r_range = (np.min(r_sub)+1e-4, 1.5)
    r, n_tot = n_contain(r_sub, r_range, N_BINS)
    
    fig, ax = plt.subplots()
    for i in range(len(r_sub_pre) - 2):
        _, n_sub = n_contain(r_sub_pre[i], r_range, N_BINS)
        ax.plot(r, n_sub / n_tot, c=colors[i])

    ax.set_xlim(r_range[0], r_range[1])
    ax.set_xlabel(r"$r/R_{\rm vir}$")
    ax.set_ylabel(r"$N_{\rm preprocessed}(<r)/N_{\rm tot}(<r)$")

    fig.savefig("../plots/preprocess_fraction.png")

    fig, ax = plt.subplots()
    for i in range(len(r_sub_pre) - 2, len(r_sub_pre)):
        _, n_sub = n_contain(r_sub_pre[i], r_range, N_BINS)
        ax.plot(r, n_sub / n_tot, c=colors[i])

    plt.plot([], [], pc("r"), label=r"${\rm LMC}$")
    plt.plot([], [], pc("b"), label=r"${\rm GSE}$")

    plt.legend(loc="upper right")

    ax.set_xlim(r_range[0], r_range[1])
    ax.set_xlabel(r"$r/R_{\rm vir}$")
    ax.set_ylabel(r"$N_{\rm preprocessed}(<r)/N_{\rm tot}(<r)$")

    fig.savefig("../plots/preprocess_fraction_2.png")

def n_contain(r, r_range, n_bins):
    n = np.zeros(n_bins+1) 
    n[1:], edges = np.histogram(r, range=r_range, bins=n_bins)
    n[0] = np.sum(r < r_range[0])
    return edges, np.cumsum(n)
    
def sub_info(dir_name):
    a = lib.scale_factors()
    
    lmc_idx, gse_idx = lib.read_merger_idxs(dir_name)
    m_idx, m = lib.read_mergers(dir_name)
    b = lib.read_branches(dir_name)
    x, mvir, snap = lib.read_tree(dir_name, ["X", "Mvir", "Snap"])
    ci = m_idx[0]
    
    surv_sub = np.where((snap[b["start"]] == 235) & b["is_main_sub"] &
                        b["is_real"])[0]
    conv = is_converged(b["start"][surv_sub], b["end"][surv_sub],
                        mvir, cut_var="mpeak")
    surv_sub = surv_sub[conv]
        
    x0 = x[b["start"][ci]]
    rvir0 = lib.mvir_to_rvir(mvir[b["start"][ci]], a[-1], OMEGA_M)

    pre_sub = b["preprocess"][surv_sub]
    start_sub, end_sub = b["start"][surv_sub], b["end"][surv_sub]
        
    x_sub = x[b["start"]][surv_sub]
    r_sub = distance(x0, x_sub)/rvir0

    prog_idx = np.ones(len(r_sub)) * -1
    prog_idx[pre_sub != -1] = 0
    for i in range(1, len(m_idx)):
        prog_idx[pre_sub == m_idx[i]] = i

    if lmc_idx != -1:
        prog_idx[pre_sub == lmc_idx] = len(m_idx)
    if gse_idx != -1:
        prog_idx[pre_sub == gse_idx] = len(m_idx) + 1

    return r_sub, prog_idx

    
def is_converged(start, end, mvir, cut_var="mvir"):
    if cut_var == "mvir":
        return mvir[start] > MVIR_CONV
    elif cut_var == "mpeak":
        mpeak = np.zeros(len(start))
        for i in range(len(start)):
            mpeak[i] = np.max(mvir[start[i]: end[i]])
        return mpeak > MVIR_CONV
    else:
        assert 0
    
def distance(x0, x):
    return np.sqrt(np.sum((x-x0)**2, axis=1))
    
        
if __name__ == "__main__": main()
