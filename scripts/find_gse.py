import numpy as np
import matplotlib.pyplot as plt
import lib
import palette
from palette import pc

DIR_FORMAT = "../tmp_data/Halo%03d"
HALO_NUMS = [4]
DIR_NAMES = [DIR_FORMAT % n for n in HALO_NUMS]
#DIR_FORMAT = "/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo%03d"
#HALO_NUMS = [4, 113, 169, 170, 222, 229, 282, 327, 349, 407, 453, 523, 625,
#             659, 666, 719, 747, 756, 788, 858, 953, 975, 983]
MP = 2.8e5
MVIR_CONV = MP * 300

def merger_snap(h, x_sub, snap_sub):    
    dist = np.sqrt(np.sum((h["x"][snap_sub] - x_sub)**2, axis=1))
    within = dist < h["rvir"][snap_sub]
    merger = h["ok"][snap_sub] & within
    if np.sum(merger) == 0:
        return -1
    else:
        return np.min(snap_sub[merger])
    
def merger_stats(b, m, x, mvir, snap):
    a = lib.scale_factors()
    sub_idx = np.where(b["is_real"] & (~b["is_disappear"]) &
                       b["is_main_sub"] & (b["preprocess"] == -1))[0]
    mw = m[0]

    ratio = np.zeros(len(sub_idx))
    scale = np.zeros(len(sub_idx))
    mpeak = np.zeros(len(sub_idx))

    mw_mass = np.max(mw["mvir"])
    
    for j, i in enumerate(sub_idx):
        mvir_i = mvir[b["start"][i]: b["end"][i]]
        snap_i = snap[b["start"][i]: b["end"][i]]
        x_i = x[b["start"][i]: b["end"][i]]
        m_snap = merger_snap(mw, x_i, snap_i)
        
        m_snap_sub = np.searchsorted(snap_i[::-1], m_snap)
        
        mpeak[j] = np.max(mvir_i)/mw_mass
        scale[j] = a[m_snap]
        ratio[j] = mvir_i[::-1][m_snap_sub]/mw["mvir"][m_snap]

    return mpeak, scale, ratio

def mw_mass_ratio(mw, a0, mass):
    i = np.searchsorted(lib.scale_factors(), a0)
    return mass/mw["mvir"][i]
    

def main():
    palette.configure(False)
    
    for dir_i, dir_name in enumerate(DIR_NAMES):
        m_idx, m = lib.read_mergers(dir_name)
        b = lib.read_branches(dir_name)
        x, mvir, snap = lib.read_tree(dir_name, ["X", "Mvir", "Snap"])
        mw = m[0]
        
        mpeak, scale, ratio = merger_stats(b, m, x, mvir, snap)
        order = np.argsort(mpeak)[-10:]
        top10 = order[-10:]

        plt.figure()
        
        plt.plot(np.log10(scale), np.log10(ratio), ".", c="k")
        plt.plot(np.log10(scale[top10]), np.log10(ratio[top10]), "o", c=pc("r"))
        plt.plot(np.log10(scale[order[-1]]), np.log10(ratio[order[-1]]),
                 "o", c=pc("b"))
        
        zlo, zhi = 2, 1.3
        alo, ahi = 1/(1 + zlo), 1/(1 + zhi)
        ylo, yhi = -6, 0.5

        rat_lo1 = mw_mass_ratio(mw, alo, 10**10.77)
        rat_lo2 = mw_mass_ratio(mw, alo, 10**11.47)
        rat_hi1 = mw_mass_ratio(mw, ahi, 10**10.90)
        rat_hi2 = mw_mass_ratio(mw, ahi, 10**11.62)
        
        plt.ylim(ylo, yhi)
        plt.fill_between(np.log10([alo, ahi]), np.log10([rat_lo1, rat_hi1]),
                         np.log10([rat_lo2, rat_hi2]), alpha=0.2, color=pc("r"))

        zlo, zhi = 0.125, 0.0
        alo, ahi = 1/(1 + zlo), 1/(1 + zhi)

        rat_lo1 = mw_mass_ratio(mw, alo, 1.485e11)
        rat_lo2 = mw_mass_ratio(mw, alo, 2.235e11)
        rat_hi1 = mw_mass_ratio(mw, ahi, 1.485e11)
        rat_hi2 = mw_mass_ratio(mw, ahi, 2.235e11)

        plt.fill_between(np.log10([alo, ahi]), np.log10([rat_lo1, rat_hi1]),
                         np.log10([rat_lo2, rat_hi2]), alpha=0.2, color=pc("b"))
        
        plt.xlabel(r"$\log_{10}(a_{\rm merger})$")
        plt.ylabel(r"$(M_{\rm sub}/M_{\rm host})(a_{\rm merger})$")
        plt.xlim((None, 0))
        
        plt.savefig("../plots/gse_search/halo%03d_mergers.png" % HALO_NUMS[dir_i])
        

if __name__ == "__main__": main()
