import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import lib
import numpy.random as random

DIR_FORMAT = "../tmp_data/Halo%03d"
HALO_NUMS = [4]
DIR_NAMES = [DIR_FORMAT % n for n in HALO_NUMS]
#DIR_FORMAT = "/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo%03d"
#HALO_NUMS = [4, 113, 169, 170, 222, 229, 282, 327, 349, 407, 453, 523, 625,
#             659, 666, 719, 747, 756, 788, 858, 953, 975, 983]
MP = 2.8e5
MVIR_CONV = MP * 300
N_SNAP = 236
NORM_BY_MPEAK = True
OMEGA_M = 0.286

def merger_snap(h1, h2):
    a = lib.scale_factors()
    
    dist = np.sqrt(np.sum((h1["x"] - h2["x"])**2, axis=1))
    within = dist < np.maximum(h1["rvir"], h2["rvir"])
    snap = np.arange(len(h1))
    merger = h1["ok"] & h2["ok"] & within
    if np.sum(merger) == 0:
        return -1
    else:
        return np.min(snap[merger])
    
def merger_times(m):
    a = lib.scale_factors()
    for i in range(1, len(m)-1):
        snap = merger_snap(m[0], m[i])
        print("%2d %.3f %.3f %.3g" %
              (i, a[snap], m["mvir"][i,snap]/m["mvir"][0,snap],
               m["mvir"][1,snap]))
        

def main():
    palette.configure(False)
    
    a = lib.scale_factors()
    npre_mah, sub1_mah, sub3_mah = [], [], []
    pre_mah, sub2_mah = [], []
    for dir_name in DIR_NAMES:
        mvir, snap = lib.read_tree(dir_name, ["Mvir", "Snap"])
        m_idx, m = lib.read_mergers(dir_name)
        b = lib.read_branches(dir_name)

        merger_times(m)
        exit(1)
        
        survives = snap[b["start"]] == N_SNAP - 1
        is_mw_sub = np.where(b["is_main_sub"] & b["is_real"] &
                             survives & (~b["is_disappear"]))[0]

        for i in is_mw_sub:
            mah = np.ones(N_SNAP)*-1
            mvir_i = mvir[b["start"][i]: b["end"][i]]
            snap_i = snap[b["start"][i]: b["end"][i]]
            mpeak = np.max(mvir_i)
            if mpeak < MVIR_CONV: continue

            if NORM_BY_MPEAK:
                mah[snap_i] = mvir_i/mpeak
            else:
                mah[snap_i] = mvir_i

            if b["preprocess"][i] == -1:
                npre_mah.append(mah)
            else:
                pre_mah.append(mah)
                if b["preprocess"][i] == m_idx[1]:
                    sub1_mah.append(mah)
                elif b["preprocess"][i] == m_idx[2]:
                    sub2_mah.append(mah)
                elif b["preprocess"][i] == m_idx[3]:
                    sub3_mah.append(mah)
                    
                    
    npre_mah = np.array(npre_mah)
    pre_mah = np.array(pre_mah)
    sub1_mah = np.array(sub1_mah)
    sub2_mah = np.array(sub2_mah)
    sub3_mah = np.array(sub3_mah)

    plot_mah(pre_mah, pc("r"))
    plot_mah(npre_mah, pc("o"))
    plot_mah(sub1_mah, pc("b", 0.75))
    plot_mah(sub2_mah, pc("b", 0.5))
    plot_mah(sub3_mah, pc("b", 0.25))

    plt.xlabel(r"$\log_{10}(a)$")
    if NORM_BY_MPEAK:
        plt.ylabel(r"$M_{\rm vir}/M_{\rm peak}$")
    else:
        plt.ylabel(r"$M_{\rm vir}\,(h^{-1}M_\odot)$")

    plt.ylim(-2, 0)
        
    plt.show()

def plot_mah(mah, c):
    a = lib.scale_factors()
    med = masked_percentile(mah, -1, 0.5)
    ok = med != -1
    plt.plot(np.log10(a[ok]), np.log10(med[ok]), c=c)
    
    low, high = bootstrap_masked_percentile(mah, -1, 0.5)
    ok = (low != -1) & (high != -1)
    plt.fill_between(np.log10(a[ok]), np.log10(low[ok]), np.log10(high[ok]),
                     alpha=0.2, color=c)

    
def masked_percentile(x, mask, p):
    out = np.ones(x.shape[1])*-1
    for i in range(x.shape[1]):
        ok = x[:,i] != mask
        if np.sum(ok) == 0: continue

        out[i] = np.percentile(x[ok,i], p*100)

    return out

def bootstrap_masked_percentile(x, mask, p):
    samples = 100
    distr = np.zeros((samples, x.shape[1]))
    for i in range(samples):
        idx = random.randint(len(x), size=len(x))
        resamp_x = x[idx,:]
        distr[i,:] = masked_percentile(resamp_x, mask, p)

        
    high = masked_percentile(distr, -1, 0.50+0.68/2)
    low = masked_percentile(distr, -1, 0.50-0.68/2)

    return low, high
    
if __name__ == "__main__": main()
