import numpy as np
import lib
import matplotlib.pyplot as plt
import palette
from palette import pc
from colossus.cosmology import cosmology
from colossus.halo import mass_so
import scipy.interpolate as interpolate
import scipy.signal as signal
import scipy.stats as stats
import binning

cosmo = cosmology.setCosmology("chinchilla",
                               {"flat": True, "H0": 70, "Om0": 0.286,
                                'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.96})

#DIR_FORMAT = "/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo%03d"
#HALO_NUMS = [4, 113, 169, 170, 222, 229, 282, 327, 349, 407, 453, 523, 625,
#             659, 666, 719, 747, 756, 788, 858, 953, 975, 983]
DIR_FORMAT = "../tmp_data/Halo%03d"
HALO_NUMS = [4, 282]
DIR_NAMES = [DIR_FORMAT % n for n in HALO_NUMS]
MP = 2.8e5
MVIR_CONV = MP * 300
HIST_BINS = 50
TIME_BINS = 40

TARGET_IDX = -1

def calc_lopsidedness(dx):
    if len(dx) <= 1: return 0
    i, j = unique_index_pairs(len(dx))
    theta = angle_between(dx[i], dx[j])
    return np.mean(theta)

def pair_counts(dx, bins):
    if len(dx) < 2: return np.zeros(bins), np.zeros(bins)

    i, j = unique_index_pairs(len(dx))
    theta = angle_between(dx[i], dx[j])
    bin_idx = np.array(bins*theta/np.pi, dtype=int)
    angle_edges = np.linspace(0, np.pi, bins+1)
    angles = (angle_edges[1:] + angle_edges[:-1])/2

    bin_idx[bin_idx < 0] = 0
    bin_idx[bin_idx >= bins] = bins - 1  
    
    return (angles, np.bincount(bin_idx, minlength=bins))

def cr_pair_counts(dx, dv, bins):
    if len(dx) < 2: return np.zeros(bins), np.zeros(bins)

    J = np.cross(dx, dv)

    i, j = unique_index_pairs(len(dx))

    J_dot_J = multi_dot(J[i], J[j])
    theta = angle_between(dx[i], dx[j])

    is_corot = J_dot_J > 0

    bin_idx = np.array(bins*theta/np.pi, dtype=int)
    angle_edges = np.linspace(0, np.pi, bins+1)
    angles = (angle_edges[1:] + angle_edges[:-1])/2
    
    return (angles, np.bincount(bin_idx, minlength=bins),
            np.bincount(bin_idx[is_corot], minlength=bins))

def normalize_pair_counts(theta, n):
    d_theta = theta[1] - theta[0]
    return (2*n) / (np.sum(n) * d_theta * np.sin(theta))

def unique_index_pairs(n):
    i, j = np.arange(n, dtype=int), np.arange(n, dtype=int)
    i, j = np.meshgrid(i, j)
    i, j = i.flatten(), j.flatten()
    return i[i < j], j[i < j]

def angle_between(x1, x2):
    scaled_dot = multi_dot(x1, x2) / np.sqrt(r_squared(x1)*r_squared(x2))
    scaled_dot[scaled_dot > 1] = 1
    scaled_dot[scaled_dot < -1] = -1
    return np.arccos(scaled_dot)

def r_squared(dx):
    r2 = np.zeros(len(dx))
    for dim in range(3):
        r2 += dx[:,dim]*dx[:,dim]
    return r2

def multi_dot(a, b):
    dims = a.shape[1]
    out = np.zeros(a.shape[0])
    for dim in range(dims):
        out += a[:,dim]*b[:,dim]
    return out

def is_contained_in(a, b):
    """ reutrn a boolean array of length |a| which returns true is a is in b.
    """
    if len(b) == 0: return np.array([], dtype=bool)
    
    b = np.sort(b)
    idx = np.searchsorted(b, a)
    idx[idx >= len(b)] = -1
    return (idx >= 0) & (b[idx] == a)

def flatten_nested_list(xs):
    out = []
    for x in xs:
        out += list(x)
    return out

def get_mpeak(b, mvir):
    mpeak = np.zeros(len(b["start"]))
    for i in range(len(b["start"])):
        mpeak[i] = np.max(mvir[b["start"][i]: b["end"][i]])
    return mpeak

def collect_indices(idx, cent_idx):
    idx = np.copy(idx)
    idx[idx == -1] = cent_idx
    bins = binning.bin_ints(idx, np.max(idx)+1)
    bins[cent_idx] = np.array([], dtype=np.int64)

    out_bins = []
    for i in range(len(bins)):
        if len(bins[i]) > 0:
            out_bins.append(np.append([i], bins[i]))

    return out_bins

def collect_by_snap(snap, start, end, idx, snap_start):    
    snap_low = np.zeros(len(idx), dtype=int)
    snap_high = np.zeros(len(idx), dtype=int)
    for i in range(len(idx)):
        start_i, end_i = start[idx[i]], end[idx[i]]
        snap_low[i] = np.min(snap[start_i: end_i])
        snap_high[i] = np.max(snap[start_i: end_i])

    out_snap = np.arange(max(np.min(snap_low), snap_start),
                         np.max(snap_high)+1, dtype=int)

    # branch index and variable index
    b_idx, v_idx = [], []
    for i in range(len(out_snap)):
        snap_i = out_snap[i]
        ok = (snap_i >= snap_low) & (snap_i <= snap_high)

        b_idx.append(idx[ok])
        v_idx.append(np.asarray(start[idx] + (snap[start[idx]] - snap_i),
                                dtype=int)[ok])
        
    return out_snap, b_idx, v_idx

def surviving_fraction(snap, b, m_snap, m_lookup, group_idxs, np_cutoffs, mvir):
    t_min, t_max, t_bins_num = 0, 10,TIME_BINS
    t_bins = np.linspace(t_min, t_max, t_bins_num+1)
    t_hists = np.zeros((len(np_cutoffs), t_bins_num))
    
    scale = lib.scale_factors()
    z = 1/scale - 1
    t_orbit = mass_so.dynamicalTime(z[m_snap], "vir", "orbit")
    age = cosmo.age(z)

    n_tots = np.zeros(len(np_cutoffs))
    
    for i in range(len(group_idxs)):
        m_idx = np.searchsorted(m_lookup, group_idxs[i][0])
        if TARGET_IDX != -1 and group_idxs[i][0] != TARGET_IDX: continue
        m_snap_i = m_snap[m_idx]
        
        group_snap, group_b_idx, group_v_idx = collect_by_snap(
            snap, b["start"], b["end"], group_idxs[i], m_snap_i)
        if len(group_snap) == 0: continue

        dt = (age[group_snap] - age[m_snap_i]) / t_orbit[m_snap_i]
        n = np.array(list(map(len, group_b_idx)))

        for j in range(len(np_cutoffs)):
            mvir_cutoff = MP * np_cutoffs[j]
            n_cutoff = np.zeros(len(group_v_idx))
            for k in range(len(group_v_idx)):
                n_cutoff[k] = np.sum(mvir[group_v_idx[k]] > mvir_cutoff)
                
            if n_cutoff[0] == 0: continue

            def safe_mean(x):
                if len(x) == 0: return 0
                return np.mean(x)
            t_count, t_edges, _ = stats.binned_statistic(
                dt, n_cutoff / n_cutoff[0], safe_mean, bins=t_bins)
            t_count[np.isnan(t_count)] = 0
            
            n_tots[j] += n_cutoff[0]
            t_hists[j] += t_count*n_cutoff[0]


    for j in range(len(t_hists)):
        t_hists[j] /= n_tots[j]
        
    t_centers = (t_edges[1:] + t_edges[:-1]) / 2
    return t_centers, t_hists, n_tots

def pair_hist(t_low, t_high, snap, b, m_snap, m_lookup, group_idxs, dx, dv):
    hist_bins = HIST_BINS
    x_hists = np.zeros((len(t_low), hist_bins))
    L_hists = np.zeros((len(t_low), hist_bins))
    
    scale = lib.scale_factors()
    z = 1/scale - 1
    t_orbit = mass_so.dynamicalTime(z[m_snap], "vir", "orbit")
    age = cosmo.age(z)
    
    for i in range(len(group_idxs)):
        m_idx = np.searchsorted(m_lookup, group_idxs[i][0])
        if TARGET_IDX != -1 and group_idxs[i][0] != TARGET_IDX: continue
        m_snap_i = m_snap[m_idx]

        group_snap, group_b_idx, group_v_idx = collect_by_snap(
            snap, b["start"], b["end"], group_idxs[i], m_snap_i)
        if len(group_snap) == 0: continue

        n = np.array(list(map(len, group_b_idx)))
        n_pairs_i = (n*(n-1))/2
        dt = (age[group_snap] - age[m_snap_i]) / t_orbit[m_snap_i]
        
        x_hists_i = np.zeros((len(dt), hist_bins))
        L_hists_i = np.zeros((len(dt), hist_bins))
        
        for j in range(len(dt)):
            if len(group_v_idx[j]) < 2: continue

            dx_j, dv_j = dx[group_v_idx[j]], dv[group_v_idx[j]]
            L_j = np.cross(dx_j, dv_j)
            theta, x_hists_i[j,:] = pair_counts(dx_j, hist_bins)
            _, L_hists_i[j,:] = pair_counts(L_j, hist_bins)
    
        for k in range(len(t_low)):
            in_t_range = np.where((dt > t_low[k]) & (dt < t_high[k]))[0]
            x_hists[k,:] += np.sum(x_hists_i[in_t_range], axis=0)
            L_hists[k,:] += np.sum(L_hists_i[in_t_range], axis=0)

    return theta, x_hists, L_hists

def lopsidedness(id, snap, b, m_snap, m_lookup, group_idxs, dx, dv):
    t_min, t_max, t_bins_num = 0, 10, TIME_BINS
    t_bins = np.linspace(t_min, t_max, t_bins_num+1)
    t_hist = np.zeros(t_bins_num)
    
    scale = lib.scale_factors()
    z = 1/scale - 1
    t_orbit = mass_so.dynamicalTime(z[m_snap], "vir", "orbit")
    age = cosmo.age(z)

    x_lop, L_lop = np.zeros(t_bins_num), np.zeros(t_bins_num)
    n_pairs = np.zeros(t_bins_num)
    
    for i in range(len(group_idxs)):
        m_idx = np.searchsorted(m_lookup, group_idxs[i][0])
        if TARGET_IDX != -1 and  group_idxs[i][0] != TARGET_IDX: continue
        m_snap_i = m_snap[m_idx]

        group_snap, group_b_idx, group_v_idx = collect_by_snap(
            snap, b["start"], b["end"], group_idxs[i], m_snap_i)
        if len(group_snap) == 0: continue

        n = np.array(list(map(len, group_b_idx)))
        n_pairs_i = (n*(n-1))/2
        dt = (age[group_snap] - age[m_snap_i]) / t_orbit[m_snap_i]

        x_lop_i, L_lop_i = np.zeros(len(dt)), np.zeros(len(dt))
        for j in range(len(x_lop_i)):
            dx_j, dv_j = dx[group_v_idx[j]], dv[group_v_idx[j]]
            id_j = id[group_v_idx[j]]
            L_j = np.cross(dx_j, dv_j)
            x_lop_i[j] = calc_lopsidedness(dx_j)
            L_lop_i[j] = calc_lopsidedness(L_j)
            
        n_pairs_hist, t_edges, _ = stats.binned_statistic(
            dt, n_pairs_i, "sum", bins=t_bins)
        
        def safe_mean(x):
            if len(x) == 0: return 0
            return np.mean(x)
        x_lop_hist, t_edges, _ = stats.binned_statistic(
            dt, x_lop_i, safe_mean, bins=t_bins)
        L_lop_hist, t_edges, _ = stats.binned_statistic(
            dt, L_lop_i, safe_mean, bins=t_bins)

        n_pairs_hist[np.isnan(n_pairs_hist)] = 0
        x_lop_hist[np.isnan(x_lop_hist)] = 0
        L_lop_hist[np.isnan(L_lop_hist)] = 0

        n_pairs += n_pairs_hist
        x_lop += x_lop_hist*n_pairs_hist
        L_lop += L_lop_hist*n_pairs_hist

    t_centers = (t_edges[1:] + t_edges[:-1]) / 2
    n_pairs[n_pairs == 0] = 1
    x_lop /= n_pairs
    L_lop /= n_pairs
    n_pairs[n_pairs == 1] = 0
    
    return t_centers, x_lop, L_lop, n_pairs

def propagate_parent_idxs(idx):
    # If you wanted to do this fast, you'd do union-find on it, but I don't
    # want to.
    
    n_changed = -1
    while n_changed != 0:
        n_changed = 0
        for i in range(len(idx)):
            if idx[i] == i: idx[i] = -1
            if idx[i] != -1 and idx[[idx[i]]] != -1:
                n_changed += 1
                idx[i] = idx[idx[i]]
  
def filter_group_idxs(pristine, group_idxs):
    first_idxs = np.array([group_idxs[i][0] for i in range(len(group_idxs))],
                           dtype=int)
    ok = is_contained_in(first_idxs, pristine)

    out = []
    for i in range(len(ok)):
        if ok[i]:
            out.append(group_idxs[i])
    return out
            
def main():
    palette.configure(False)

    np_cutoffs = [0, 30, 100, 300, 1000]
    colors = [pc("r"), pc("o"), pc("g"), pc("b"), pc("p")]

    t_names = ["0", "2", "4", "6", "8"]
    t_low =  np.array([-0.25, 1.75, 3.75, 5.75, 7.75])
    t_high = np.array([ 0.25, 2.25, 4.25, 6.25, 8.25]) 

    n_sat_tot = 0
    n_pair_tot = 0
    f_surv = np.zeros((len(t_low), TIME_BINS))
    x_lop = np.zeros(TIME_BINS)
    L_lop = np.zeros(TIME_BINS)
    x_hists = np.zeros((len(t_low), HIST_BINS))
    L_hists = np.zeros((len(t_low), HIST_BINS))
    
    for dir_i, dir_name in enumerate(DIR_NAMES):
        halo_name = dir_name.split("/")[-1]
        print(halo_name)

        m_idx, m = lib.read_mergers(dir_name)
        b = lib.read_branches(dir_name)
        id, x, v, mvir, snap = lib.read_tree(
            dir_name, ["ID", "X", "V", "Mvir", "Snap"])
        mpeak = get_mpeak(b, mvir)
        mw_idx, mw = m_idx[0], m[0]

        print(m_idx)
        
        dx = x - m[0]["x"][snap]
        dv = v - m[0]["v"][snap]

        propagate_parent_idxs(b["preprocess"])
        pristine = lib.pristine_merger_indices(b)
        
        # Remember, this only returns results for pristine haloes.
        _, m_snap, ratio = lib.merger_stats(b, m, x, mvir, snap)
        
        group_idxs = collect_indices(
            np.asarray(b["preprocess"], dtype=int), mw_idx)
        group_idxs = filter_group_idxs(pristine, group_idxs)        
        
        t_centers, f_surv_i, n_sat_tot_i = surviving_fraction(
            snap, b, m_snap, pristine, group_idxs, np_cutoffs, mvir)
        
        t_centers, x_lop_i, L_lop_i, n_pair_tot_i = lopsidedness(
            id, snap, b, m_snap, pristine, group_idxs, dx, dv)

        theta, x_hists_i, L_hists_i = pair_hist(
            t_low, t_high, snap, b, m_snap, pristine, group_idxs, dx, dv)

        for k in range(len(n_sat_tot_i)):
            f_surv[k,:] += f_surv_i[k,:]*n_sat_tot_i[k]
        x_lop += x_lop_i*n_pair_tot_i
        L_lop += L_lop_i*n_pair_tot_i
        x_hists += x_hists_i
        L_hists += L_hists_i
        n_sat_tot += n_sat_tot_i
        n_pair_tot += n_pair_tot_i

        plt.figure(0)
        for i in range(len(np_cutoffs)):
            plt.plot(t_centers, f_surv_i[i,:], colors[i],
                     label=r"$N_{\rm vir,min} = %d$" % np_cutoffs[i])

        plt.figure(1)
        plt.plot(t_centers, x_lop_i, pc("r"), label=r"$\vec{x}$")
        plt.plot(t_centers, L_lop_i, pc("b"), label=r"$\vec{L}$")

        plt.figure(2)
        for i in range(len(x_hists)):
            x_norm = normalize_pair_counts(theta, x_hists_i[i,:])
            plt.plot(theta, x_norm, colors[i],
                     label=(r"$t - t_{\rm infall} = %s\,t_{\rm orbit}$" %
                            t_names[i]))

        plt.figure(3)
        for i in range(len(L_hists)):
            L_norm = normalize_pair_counts(theta, L_hists_i[i,:])
            plt.plot(theta, L_norm, colors[i],
                     label=(r"$t - t_{\rm infall} = %s\,t_{\rm orbit}$" %
                            t_names[i]))

        plt.figure(0)
        plt.ylabel(r"$f_{\rm survive} (N_{\rm vir} > N_{\rm vir,min})$")
        plt.xlabel(r"$(t - t_{\rm infall}) / t_{\rm orbit}$")
        plt.legend(loc="upper right", fontsize=16)
        plt.yscale("log")
        plt.savefig("../plots/lopsidedness/f_survive_vs_t.%s.png" % halo_name)
        plt.clf()

        plt.figure(1)
        plt.ylabel(r"$\langle\theta\rangle\ ({\rm radians})$")
        plt.xlabel(r"$(t - t_{\rm infall}) / t_{\rm orbit}$")
        plt.legend(loc="upper right", fontsize=16)
        lo, hi = plt.xlim()
        plt.xlim(lo, hi)
        plt.plot([lo, hi], [np.pi/2]*2, "--", c="k", lw=1.5)
        plt.savefig("../plots/lopsidedness/lop_vs_t.%s.png" % halo_name)
        plt.clf()

        plt.figure(2)
        plt.xlabel(r"$\theta_x$")
        plt.ylabel(r"$N(\theta_x) / N_{\rm isotropic}(\theta_x)$")
        plt.legend(loc="upper right", fontsize=16)
        lo, hi = plt.xlim()
        plt.xlim(lo, hi)
        plt.plot([lo, hi], [1]*2, "--", c="k", lw=1.5)
        plt.ylim(0, None)
        plt.savefig("../plots/lopsidedness/theta_x_t_hist.%s.png" % halo_name)
        plt.clf()

        plt.figure(3)
        plt.xlabel(r"$\theta_L$")
        plt.ylabel(r"$N(\theta_L) / N_{\rm isotropic}(\theta_L)$")
        plt.legend(loc="upper right")
        plt.legend(loc="upper right", fontsize=16)
        lo, hi = plt.xlim()
        plt.xlim(lo, hi)
        plt.plot([lo, hi], [1]*2, "--", c="k", lw=1.5)
        plt.ylim(0, None)
        plt.savefig("../plots/lopsidedness/theta_L_t_hist.%s.png" % halo_name)
        plt.clf()
            
        
    for k in range(len(n_sat_tot_i)):
        f_surv[k,:] /= n_sat_tot[k]
    x_lop /= n_pair_tot
    L_lop /= n_pair_tot
        
    plt.figure(0)
    for i in range(len(np_cutoffs)):
        plt.plot(t_centers, f_surv[i,:], colors[i],
                 label=r"$N_{\rm vir,min} = %d$" % np_cutoffs[i])

    plt.figure(1)
    plt.plot(t_centers, x_lop, pc("r"), label=r"$\vec{x}$")
    plt.plot(t_centers, L_lop, pc("b"), label=r"$\vec{L}$")

    plt.figure(2)
    for i in range(len(x_hists)):
        x_norm = normalize_pair_counts(theta, x_hists[i,:])
        plt.plot(theta, x_norm, colors[i],
                 label=(r"$t - t_{\rm infall} = %s\,t_{\rm orbit}$" %
                        t_names[i]))

    plt.figure(3)
    for i in range(len(L_hists)):
        L_norm = normalize_pair_counts(theta, L_hists[i,:])
        plt.plot(theta, L_norm, colors[i],
                 label=(r"$t - t_{\rm infall} = %s\,t_{\rm orbit}$" %
                        t_names[i]))
            
    plt.figure(0)
    plt.ylabel(r"$f_{\rm survive} (N_{\rm vir} > N_{\rm vir,min})$")
    plt.xlabel(r"$(t - t_{\rm infall}) / t_{\rm orbit}$")
    plt.legend(loc="upper right", fontsize=16)
    plt.yscale("log")
    plt.savefig("../plots/lopsidedness/f_survive_vs_t.png")
    
    plt.figure(1)
    plt.ylabel(r"$\langle\theta\rangle\ ({\rm radians})$")
    plt.xlabel(r"$(t - t_{\rm infall}) / t_{\rm orbit}$")
    plt.legend(loc="upper right", fontsize=16)
    lo, hi = plt.xlim()
    plt.xlim(lo, hi)
    plt.plot([lo, hi], [np.pi/2]*2, "--", c="k", lw=1.5)
    plt.savefig("../plots/lopsidedness/lop_vs_t.png")
    
    plt.figure(2)
    plt.xlabel(r"$\theta_x$")
    plt.ylabel(r"$N(\theta_x) / N_{\rm isotropic}(\theta_x)$")
    plt.legend(loc="upper right", fontsize=16)
    lo, hi = plt.xlim()
    plt.xlim(lo, hi)
    plt.plot([lo, hi], [1]*2, "--", c="k", lw=1.5)
    plt.ylim(0, None)
    plt.savefig("../plots/lopsidedness/theta_x_t_hist.png")
    
    plt.figure(3)
    plt.xlabel(r"$\theta_L$")
    plt.ylabel(r"$N(\theta_L) / N_{\rm isotropic}(\theta_L)$")
    plt.legend(loc="upper right")
    plt.legend(loc="upper right", fontsize=16)
    lo, hi = plt.xlim()
    plt.xlim(lo, hi)
    plt.plot([lo, hi], [1]*2, "--", c="k", lw=1.5)
    plt.ylim(0, None)
    plt.savefig("../plots/lopsidedness/theta_L_t_hist.png")

    plt.show()
    
if __name__ == "__main__": main()
