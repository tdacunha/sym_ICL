import numpy as np
import matplotlib.pyplot as plt
import lib
import palette
from palette import pc

#DIR_FORMAT = "../tmp_data/Halo%03d"
#HALO_NUMS = [4]
DIR_FORMAT = "/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo%03d"
HALO_NUMS = [4, 113, 169, 170, 222, 229, 282, 327, 349, 407, 453, 476, 523, 625,
             659, 666, 719, 747, 756, 788, 858, 908, 953, 975, 983]
DIR_NAMES = [DIR_FORMAT % n for n in HALO_NUMS]
MP = 2.8e5
MVIR_CONV = MP * 300

DEVESHI_LMC_ID_MAP = {
    453: 131024491,
    222: 24112583,
    659: 30494064,
    756: 170087224,
    229: 132361106,
    4: 7910278,
    719: 18160481,
    349: 67042798,
    858: 64499712,
    523: 36478940,
    327: 9434628,
    975: 40731786,
    113: 13827755,
    983: 179606523,
    170: 45488338,
    407: 53046731,
    666: 26531868,
    747: 28384953,
    953: 54039443,
    625: 121331687,
    788: 21696396,
    282: 20237190,
    169: 39029795,
    476: -1,
    908: -1,
}

def merger_snap(h, x_sub, snap_sub):    
    dist = np.sqrt(np.sum((h["x"][snap_sub] - x_sub)**2, axis=1))
    within = dist < h["rvir"][snap_sub]
    merger = h["ok"][snap_sub] & within
    if np.sum(merger) == 0:
        return -1
    else:
        return np.min(snap_sub[merger])
    
def merger_stats(b, m, m_idx, x, mvir, snap):
    a = lib.scale_factors()
    sub_idx = np.where(b["is_real"] & (~b["is_disappear"]) &
                       b["is_main_sub"] & (b["preprocess"] == -1))[0]
    mw = m[0]

    ratio = np.zeros(len(sub_idx))
    scale = np.zeros(len(sub_idx))
    mpeak = np.zeros(len(sub_idx))
    mpeak_raw = np.zeros(len(sub_idx))

    mw_mass = np.max(mw["mvir"])
    for j, i in enumerate(sub_idx):
        mvir_i = mvir[b["start"][i]: b["end"][i]]
        snap_i = snap[b["start"][i]: b["end"][i]]
        x_i = x[b["start"][i]: b["end"][i]]
        m_snap = merger_snap(mw, x_i, snap_i)
        
        m_snap_sub = np.searchsorted(snap_i[::-1], m_snap)
        mpeak_raw[j] = np.max(mvir_i)
        mpeak[j] = np.max(mvir_i)/mw_mass
        scale[j] = a[m_snap]
        ratio[j] = mvir_i[::-1][m_snap_sub]/mw["mvir"][m_snap]

    return mpeak, scale, ratio, sub_idx

def mw_mass_ratio(mw, a0, mass):
    i = np.searchsorted(lib.scale_factors(), a0)
    return mass/mw["mvir"][i]
    
def find_deveshi_lmc(b, id, sub_idx, halo_idx):
    target_id = DEVESHI_LMC_ID_MAP[HALO_NUMS[halo_idx]]
    i = np.where(id[b["start"][sub_idx]] == target_id)[0]
    if len(i) == 0: return -1
    return i[0]

def is_GSE(mpeak, ratio, merger_scale, snap, b, sub_idx):
    a = lib.scale_factors()
    z = 1/a - 1
    disrupt_z = z[snap[b["start"][sub_idx]]]
    
    candidates = ((disrupt_z > 0.67) & (disrupt_z < 3) & (ratio > 1/5.0))

    if np.sum(candidates) == 0:
        return candidates
    return candidates & (mpeak >= np.max(mpeak[candidates]))

def is_LMC(mpeak, ratio, merger_scale, snap, b, sub_idx):
    merger_z = 1/merger_scale - 1
    candidates = ((merger_z < 0.25) & (ratio > 1/10.0))

    if np.sum(candidates) == 0:
        return candidates
    return candidates & (mpeak >= np.max(mpeak[candidates]))

def main():
    palette.configure(False)
    
    plt.figure()

    for dir_i, dir_name in enumerate(DIR_NAMES):
        print(dir_name.split("/")[-1], end=" ")
        m_idx, m = lib.read_mergers(dir_name)
        b = lib.read_branches(dir_name)
        x, mvir, snap, id, vmax = lib.read_tree(
            dir_name, ["X", "Mvir", "Snap", "ID", "Vmax"])
        mw = m[0]

        mpeak, scale, ratio, sub_idx = merger_stats(b, m, m_idx, x, mvir, snap)
        order = np.argsort(mpeak)[-10:]
        top10 = order[-10:]
        
        survive = snap[b["start"][sub_idx]] == 235
        i_lmc_deveshi = find_deveshi_lmc(b, id, sub_idx, dir_i)
        gse = is_GSE(mpeak, ratio, scale, snap, b, sub_idx)
        lmc = is_LMC(mpeak, ratio, scale, snap, b, sub_idx)
        
        idx = np.arange(len(lmc))
        gse_idx, lmc_idx = idx[gse], idx[lmc]

        print("%7d %3d %9d" % (sub_idx[i_lmc_deveshi] if i_lmc_deveshi != -1 else -1,
                               235, DEVESHI_LMC_ID_MAP[HALO_NUMS[dir_i]]),
                               end=" ")

        if len(lmc_idx) == 0:
            print("%7d %3d %9d" % (-1, -1, -1), end=" ")
        else:
            print("%7d %3d %9d" % (sub_idx[lmc_idx[0]],
                                   snap[b["start"][sub_idx[lmc_idx[0]]]],
                                   id[b["start"][sub_idx[lmc_idx[0]]]]), end=" ")

        if len(gse_idx) == 0:
            print("%7d %3d %9d" % (-1, -1, -1))
        else:
            print("%7d %3d %9d" % (sub_idx[gse_idx[0]],
                                   snap[b["start"][sub_idx[gse_idx[0]]]],
                                   id[b["start"][sub_idx[gse_idx[0]]]]))

        plt.plot(np.log10(scale[~survive]), np.log10(ratio[~survive]),
                 ".", c=pc("a"))
        plt.plot(np.log10(scale[survive]), np.log10(ratio[survive]),
                 ".", c=pc("k"))
        plt.plot(np.log10(scale[gse]), np.log10(ratio[gse]),
                 "o", c=pc("r"))
        plt.plot(np.log10(scale[lmc]), np.log10(ratio[lmc]),
                 "o", c=pc("b"))
        if i_lmc_deveshi != -1:
            plt.plot([np.log10(scale[i_lmc_deveshi])],
                     [np.log10(ratio[i_lmc_deveshi])],
                     "*", c=pc("p"))
        
        zlo, zhi = 3, 2
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
        plt.clf()
        

if __name__ == "__main__": main()
