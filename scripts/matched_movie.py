import symlib
from palette import pc
import palette
import numpy as np
import matplotlib.pyplot as plt
import os.path as path
import scipy.spatial as spatial

def pre_infall_mpeak(h, hist):
    n_snap = h.shape[1]
    snap = np.arange(n_snap, dtype=int)

    out = np.zeros(len(hist))
    for i in range(len(h)):
        infall_snap = hist["first_infall_snap"][i]
        if infall_snap == -1: infall_snap = n_snap
        out[i] = np.max(h[i,snap<=infall_snap]["mvir"])

    return out

def plot_matched_haloes(ax, h_1, h_2, match_1, match_2, snap, r_max,
                        plot_set_1=None, plot_set_2=None):
    if plot_set_1 is None: plot_set_1 = np.arange(len(h_1), dtype=int)
    if plot_set_2 is None: plot_set_2 = np.arange(len(h_2), dtype=int)

    ls_1, ls_2 = "-", "--"
    main_c = pc("a")
    surv_c, dis_c, unmatch_c = pc("b"), pc("r"), pc("k")
    print(snap)

    for i in plot_set_1:
        if not h_1["ok"][i,snap]: continue

        m = match_1[i]

        if i == 0:
            c = main_c
        elif m == -1:
            c = unmatch_c
        elif not h_2[m,snap]["ok"]:
            c = dis_c
        else:
            c = surv_c
            
        symlib.plot_circle(ax, h_1["x"][i,snap,0],
                           h_1["x"][i,snap,1],
                           h_1["rvir"][i,snap],
                           c=c, ls=ls_1, lw=1.5)
        plt.plot([h_1["x"][i,snap,0]], [h_1["x"][i,snap,1]], ".", c=main_c)

        if c == surv_c:
            symlib.plot_circle(ax, h_2["x"][m,snap,0],
                               h_2["x"][m,snap,1],
                               h_2["rvir"][m,snap],
                               c=c, ls=ls_2, lw=1.5)
            plt.plot([h_2["x"][m,snap,0]], [h_2["x"][m,snap,1]], ".", c=main_c)
            plt.plot([h_2["x"][m,snap,0], h_1["x"][i,snap,0]],
                     [h_2["x"][m,snap,1], h_1["x"][i,snap,1]],
                     lw=1, c=main_c)

    for i in plot_set_2:
        if not h_2["ok"][i,snap]: continue

        m = match_2[i]

        if i == 0:
            c = main_c
        elif m == -1:
            c = unmatch_c
        elif not h_1[m,snap]["ok"]:
            c = dis_c            
        else:
            c = surv_c

        symlib.plot_circle(ax, h_2["x"][i,snap,0],
                           h_2["x"][i,snap,1],
                           h_2["rvir"][i,snap],
                           c=c, ls=ls_2, lw=1.5)
        plt.plot([h_2["x"][i,snap,0]], [h_2["x"][i,snap,1]], ".", c=main_c)

        if c == surv_c:
            symlib.plot_circle(ax, h_1["x"][m,snap,0],
                               h_1["x"][m,snap,1],
                               h_1["rvir"][m,snap],
                               c=c, ls=ls_1, lw=1.5)
            ax.plot([h_1["x"][m,snap,0]], [h_1["x"][m,snap,1]], ".", c=main_c)
            ax.plot([h_1["x"][m,snap,0], h_2["x"][i,snap,0]],
                    [h_1["x"][m,snap,1], h_2["x"][i,snap,1]],
                    lw=1, c=main_c)


    ax.set_xlim(-r_max, +r_max)
    ax.set_ylim(-r_max, +r_max)
    ax.set_xlabel(r"$X\ ({\rm kpc})$")
    ax.set_ylabel(r"$Y\ ({\rm kpc})$")
    
class BinnedKDTrees(object):
    def __init__(self, x, mpeak, ok, bin_factor=1.25):
        log_mpeak = np.log10(mpeak)
        log_m_min, log_m_max = np.min(log_mpeak), np.max(log_mpeak)
        d_log_m = np.log10(bin_factor)

        n_bins = int(np.ceil((log_m_max - log_m_min) / d_log_m))

        self.bins = np.linspace(log_m_min, log_m_min+(n_bins+1)*d_log_m, n_bins+1)
        self.trees = [None]*n_bins
        self.idxs = [None]*n_bins

        for i in range(n_bins):
            ok_i = ok & (log_mpeak >= self.bins[i]) & (log_mpeak < self.bins[i+1])
            if np.sum(ok_i) > 0:
                self.idxs[i] = np.where(ok_i)[0]
                self.trees[i] = spatial.KDTree(x[ok_i])

    def best_match(self, x, mpeak, search_factor=2):
        m_low = np.log10(mpeak/search_factor)
        m_high = np.log10(mpeak*search_factor)
        
        idx_high = min(np.searchsorted(self.bins, m_high), len(self.bins)-1)
        idx_low = max(np.searchsorted(self.bins, m_low), 0)

        match_idxs, match_rs = [], []

        for i in range(idx_low, idx_high):
            if self.trees[i] is None: continue
            r, j = self.trees[i].query(x)
            match_rs.append(r)
            match_idxs.append(self.idxs[i][j])

        if len(match_rs) == 0:
            return -1

        return match_idxs[np.argmin(match_rs)]
        

def match_subhaloes(h_1, h_2, hist_1, hist_2, min_votes=4):
    n_snap = h_1.shape[1]
    mpeak_1 = pre_infall_mpeak(h_1, hist_1)
    mpeak_2 = pre_infall_mpeak(h_2, hist_2)

    votes_1, votes_2 = [], []

    for snap in range(n_snap):
        match_1 = np.ones(len(h_1), dtype=int)*-1
        match_2 = np.ones(len(h_2), dtype=int)*-1

        ok_1 = (snap < hist_1["first_infall_snap"]) & h_1["ok"][:,snap]
        ok_2 = (snap < hist_2["first_infall_snap"]) & h_2["ok"][:,snap]
        trees_1 = BinnedKDTrees(h_1["x"][:,snap,:], mpeak_1, ok_1)
        trees_2 = BinnedKDTrees(h_2["x"][:,snap,:], mpeak_2, ok_2)

        for i in np.where(ok_1)[0]:
            match_1[i] = trees_2.best_match(h_1["x"][i,snap], mpeak_1[i])
        for i in np.where(ok_2)[0]:
            match_2[i] = trees_1.best_match(h_2["x"][i,snap], mpeak_2[i])
                    
        for i in range(len(match_1)):
            if match_2[match_1[i]] == i:
                votes_1.append(i)
                votes_2.append(match_1[i])

    votes_1 = np.array(votes_1, dtype=int)
    votes_2 = np.array(votes_2, dtype=int)        

    vote_count = { }

    for i in range(len(votes_1)):
        pair = (votes_1[i], votes_2[i])
        vote_count[pair] = 1 + vote_count.get(pair, 0)

    idx_1, idx_2, n_votes = process_vote_count(vote_count)
    order = np.argsort(n_votes)[::-1]
    idx_1, idx_2, n_votes = idx_1[order], idx_2[order], n_votes[order]

    match_1 = np.ones(len(h_1), dtype=int)*-1
    match_2 = np.ones(len(h_2), dtype=int)*-1
    in_use_1, in_use_2 = {}, {}
    for i in range(len(n_votes)):
        if idx_1[i] in in_use_1 or idx_2[i] in in_use_2 or n_votes[i] < min_votes:
            continue
        i1, i2 = idx_1[i], idx_2[i]
        match_1[i1] = i2
        match_2[i2] = i1
        in_use_1[i1] = None
        in_use_2[i2] = None
    
    return match_1, match_2
    
def process_vote_count(vote_count):
    pairs = vote_count.keys()
    idx_1, idx_2 = list(zip(*pairs))
    idx_1, idx_2 = np.array(idx_1, dtype=int), np.array(idx_2, dtype=int)
    n_votes = np.array([vote_count[pair] for pair in pairs], dtype=int)
    return idx_1, idx_2, n_votes

def main():
    palette.configure(False)

    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    out_dir = "/home/users/phil1/code/src/github.com/phil-mansfield/symphony_pipeline/plots/core_plots/matches/"
    suite_hr, suite_lr = "SymphonyMilkyWayHR", "SymphonyMilkyWayLR"
    n_host = symlib.n_hosts(suite_hr)

    colors = [pc("k"), pc("r", 0.4), pc("o", 0.4), pc("g", 0.4),
              pc("b", 0.4), pc("p", 0.4), pc("r", 0.8), pc("o", 0.8),
              pc("g", 0.8), pc("b", 0.8)]

    f = open("tables/lr_hr_match.txt", "w+")
    print("""# 0 - host index
    1 - suite (0 - LR, 1 - HR)
2 - LR index
3 - HR index""", file=f)

    for i_host in range(n_host):
        sim_dir_lr = symlib.get_host_directory(base_dir, suite_lr, i_host)
        sim_dir_hr = symlib.get_host_directory(base_dir, suite_hr, i_host)
        h_lr, hist_lr = symlib.read_subhalos(sim_dir_lr, comoving=True)
        h_hr, hist_hr = symlib.read_subhalos(sim_dir_hr, comoving=True)
        hist_lr[0]["first_infall_snap"] = h_lr.shape[1]
        hist_hr[0]["first_infall_snap"] = h_hr.shape[1]
        h_phys_lr, _ = symlib.read_subhalos(sim_dir_lr)
        h_phys_hr, _ = symlib.read_subhalos(sim_dir_hr)

        n_max_lr, n_max_hr = len(h_lr), len(h_hr)
        h_lr, h_hr = h_lr[:n_max_lr], h_hr[:n_max_hr]
        h_phys_lr, h_phys_hr = h_phys_lr[:n_max_lr], h_phys_hr[:n_max_hr]
        hist_lr, hist_hr = hist_lr[:n_max_lr], hist_hr[:n_max_hr]

        targets_lr = np.where(hist_lr["mpeak"]/hist_lr["mpeak"][0] > 1e-3)[0]
        targets_hr = np.where(hist_hr["mpeak"]/hist_hr["mpeak"][0] > 1e-3)[0]

        match_lr, match_hr = match_subhaloes(
            h_lr, h_hr, hist_lr, hist_hr
        )

        for i_lr, i_hr in enumerate(match_lr):
            print("%d %d %4d %4d" % (i_host, 0, i_lr, i_hr), file=f)
        for i_hr, i_lr in enumerate(match_hr):
            print("%d %d %4d %4d" % (i_host, 1, i_lr, i_hr), file=f)

        pre_mpeak_lr = pre_infall_mpeak(h_lr, hist_lr)
        pre_mpeak_hr = pre_infall_mpeak(h_hr, hist_hr)
        
        fig, ax = plt.subplots()
        snap = np.arange(h_hr.shape[1], dtype=int)
        scale = symlib.scale_factors(suite_hr)
        #print(i_host, len(h_hr), len(h_lr))
        for j in range(len(colors)):
            j_lr, j_hr = match_hr[j], j
            #print("   ", j_lr, j_hr)

            ok1_hr = h_hr["ok"][j_hr,:]
            ok2_hr = ok1_hr & (snap < hist_hr["first_infall_snap"][j_hr])

            ax.plot(scale[ok1_hr], h_hr["mvir"][j_hr,ok1_hr], c=colors[j], lw=1.5)
            ax.plot(scale[ok2_hr], h_hr["mvir"][j_hr,ok2_hr], c=colors[j])

            if j_lr != -1:
                ok1_lr = h_lr["ok"][j_lr,:]
                ok2_lr = ok1_lr & (snap < hist_lr["first_infall_snap"][j_lr])
                ax.plot(scale[ok1_lr], h_lr["mvir"][j_lr,ok1_lr],
                        "--", c=colors[j], lw=1.5)
                ax.plot(scale[ok2_lr], h_lr["mvir"][j_lr,ok2_lr],
                        "--", c=colors[j])

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$a(t)$")
        ax.set_ylabel(r"$M_{\rm vir}$")
        fig.savefig(path.join(out_dir, "mass_match_%d.png" % i_host))

        continue

        r_max = 1.25*h_phys_hr["rvir"][0,-1]
        n_snap = h_phys_hr.shape[1]

        fig, ax = plt.subplots()

        for snap in range(n_snap):
            if snap % 50 == 0: print("   snap %d" % snap)
            ax.cla()

            plot_matched_haloes(
                ax, h_phys_hr, h_phys_lr, match_hr, match_lr,
                snap, r_max, plot_set_1=targets_hr, plot_set_2=targets_lr
            )

            fig.savefig(path.join(out_dir, "h%d" % i_host, 
                                  "z0_match_%03d.png" % snap))
            
            if snap == n_snap-1:
                fig.savefig(path.join(out_dir, "h%d" % i_host, 
                                      "z0_match_%03d.png" % i_host))
        

if __name__ == "__main__": main()
