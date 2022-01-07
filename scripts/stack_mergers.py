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

cosmo = cosmology.setCosmology("chinchilla",
                               {"flat": True, "H0": 70, "Om0": 0.286,
                                'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.96})

DIR_FORMAT = "/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo%03d"
HALO_NUMS = [4, 113, 169, 170, 222, 229, 282, 327, 349, 407, 453, 523, 625,
             659, 666, 719, 747, 756, 788, 858, 953, 975, 983]
#DIR_FORMAT = "../tmp_data/Halo%03d"
#HALO_NUMS = [4]
DIR_NAMES = [DIR_FORMAT % n for n in HALO_NUMS]
MP = 2.8e5
MVIR_CONV = MP * 300

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

def main():
    palette.configure(False)

    log_ratio_min, log_ratio_max = -1.5, 0.0
    log_ratio_bins = 6
    log_ratio_edges = np.linspace(log_ratio_min, log_ratio_max,
                                  log_ratio_bins+1)
    dlog_ratio = (log_ratio_max - log_ratio_min) / log_ratio_bins
    scale = lib.scale_factors()
    
    merger_counts = np.zeros(log_ratio_bins)
    sub_counts = np.zeros(log_ratio_bins)

    trj_ms = [[] for _ in range(log_ratio_bins)]
    trj_rs = [[] for _ in range(log_ratio_bins)]
    trj_ts = [[] for _ in range(log_ratio_bins)]
    
    for dir_i, dir_name in enumerate(DIR_NAMES):
        halo_name = dir_name.split("/")[-1]
        print(halo_name)

        m_idx, m = lib.read_mergers(dir_name)
        b = lib.read_branches(dir_name)
        x, mvir, snap = lib.read_tree(dir_name, ["X", "Mvir", "Snap"])
        mpeak = get_mpeak(b, mvir)
        mw = m[0]

        p_sub_idx = lib.pristine_merger_indices(b)
        _, m_snap, ratio = lib.merger_stats(b, m, x, mvir, snap)
        log_ratio = np.log10(ratio)
        ok = mpeak[p_sub_idx] > MVIR_CONV

        merger_counts += np.histogram(
            np.log10(ratio[ok]), bins=log_ratio_edges)[0]
        
        for i in range(log_ratio_bins):
            low, high = log_ratio_edges[i], log_ratio_edges[i+1]
            in_range = (log_ratio >= low) & (log_ratio < high)
            p_sub_idx_i = p_sub_idx[in_range & ok]
            m_snap_i = m_snap[in_range & ok]

            trj_t, trj_r, trj_m = trj_ts[i], trj_rs[i], trj_ms[i]
            
            for j in range(len(p_sub_idx_i)):
                trj_t_i, trj_r_i, trj_m_i = sub_sub_trajectories(
                    mw, scale, b, x, mvir, snap,
                    mpeak, p_sub_idx_i[j], m_snap_i[j])
                trj_t.append(trj_t_i)
                trj_r.append(trj_r_i)
                trj_m.append(trj_m_i)

    for i in range(log_ratio_bins):
        low, high = log_ratio_edges[i], log_ratio_edges[i+1]
        trj_t, trj_r, trj_m = trj_ts[i], trj_rs[i], trj_ms[i]
        
        trj_t = flatten_nested_list(trj_t)
        trj_r = flatten_nested_list(trj_r)
        trj_m = flatten_nested_list(trj_m)
            
        t_all = np.array(flatten_nested_list(trj_t))
        r_all = np.array(flatten_nested_list(trj_r))
        m_all = np.array(flatten_nested_list(trj_m))
            
        sub_counts[i] = len(trj_r)
            
        plt.figure(2*i)
        for j in range(len(trj_t)):
            plt.plot(trj_t[j], trj_m[j], lw=1, alpha=0.25, c="k")

        plot_contour(t_all, m_all, pc("r"), (-1.5, 1.5), 20)
        
        plt.title(r"$%.2f< \log_{10}(M_{\rm sub}/M_{\rm host}) <%.2f$" %
                  (low, high))
        plt.xlabel(r"$(t - t_{\rm infall})/t_{\rm orbit}$")
        plt.ylabel(r"$M_{\rm sub}/M_{\rm sub,peak}$")
        plt.xlim(-1.5, 1.5)
        plt.ylim(-2, 0)
        plt.plot([0, 0], [-2, 0], "--", c="k")

        plt.savefig("../plots/merger_stack_m_%d.png")
        
        plt.figure(2*i + 1)
        for j in range(len(trj_t)):
            plt.plot(trj_t[j], trj_r[j], lw=1, alpha=0.25, c="k")

        plot_contour(t_all, r_all, pc("b"), (-1.5, 1.5), 20)
        
        plt.title(r"$%.2f< r/R_{\rm host} <%.2f$" %
                  (low, high))
        plt.xlabel(r"$(t - t_{\rm infall})/t_{\rm orbit}$")
        plt.ylabel(r"$r/R_{\rm host}$")
        plt.xlim(-1.5, 1.5)
        plt.ylim(0, 5)
        plt.plot([0, 0], [0, 5], "--", c="k")

        plt.savefig("../plots/merger_stack_m_%d.png")

    #plt.show()

def plot_contour(x, y, color, x_range, bins):
    n, edges, _ = stats.binned_statistic(
        x, y, "count", range=x_range, bins=bins)
    mid = (edges[1:] + edges[:-1]) / 2
    
    med, _, _ = stats.binned_statistic(
        x, y, "median", range=x_range, bins=bins)
    lo, _, _ = stats.binned_statistic(
        x, y, lambda xx: np.percentile(xx, 50-68/2), range=x_range, bins=bins)
    hi, _, _ = stats.binned_statistic(
        x, y, lambda xx: np.percentile(xx, 50+68/2), range=x_range, bins=bins)

    
    ok = n > 5
    plt.plot(mid[ok], med[ok], color)
    plt.fill_between(mid[ok], lo[ok], hi[ok], alpha=0.25, color=color)
    
def sub_sub_trajectories(mw, scale, b, x, mvir, snap,
                         mpeak, grp_idx, m_snap):
    z = 1/scale - 1
    t_orbit = mass_so.dynamicalTime(z[m_snap], "vir", "orbit")
    rvir = mass_so.M_to_R(mw["mvir"][m_snap], z[m_snap], "vir")/1e3
    
    age = cosmo.age(z)
    dt = age - age[m_snap]

    sub_idx = np.where((mpeak > MVIR_CONV) & (b["preprocess"] == grp_idx))[0]

    n = len(sub_idx)
    trj_t, trj_r, trj_m = [None]*n, [None]*n, [None]*n

    for i in range(n):
        j = sub_idx[i]

        start, end = b["start"][j], b["end"][j]
        
        trj_t[i] = dt[snap[start: end]] / t_orbit
        trj_m[i] = np.log10(mvir[start: end] / mpeak[j])
        dx = x[start: end] - mw["x"][snap[start: end]]
        trj_r[i] = np.sqrt(np.sum(dx**2, axis=1)) / rvir
        
    return trj_t, trj_r, trj_m
            
if __name__ == "__main__": main()
