import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import lib
import numpy.random as random
from colossus.cosmology import cosmology

#DIR_FORMAT = "../tmp_data/Halo%03d"
#HALO_NUMS = [4]
DIR_FORMAT = "/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo%03d"
HALO_NUMS = [4, 113, 169, 170, 222, 229, 282, 327, 349, 407, 453, 467, 523, 625,
             659, 666, 719, 747, 756, 788, 858, 953, 975, 983, 908]

DIR_NAMES = [DIR_FORMAT % n for n in HALO_NUMS]
MP = 2.8e5
MVIR_CONV = MP * 300
N_SNAP = 236
NORM_BY_MPEAK = True
OMEGA_M = 0.286
LOG_SCALING = False

params = {'flat': True, 'H0': 70.0, 'Om0': 0.286, 'Ob0': 0.049, 'sigma8': 0.82, 'ns': 1.0}
cosmo = cosmology.setCosmology("", params)

def main():
    palette.configure(False)
    
    a = lib.scale_factors()
    npre_mah, lmc_cosmo_mah, lmc_zoom_mah = [], [], []
    pre_mah, gse_mah = [], []
    for dir_name in DIR_NAMES:
        print(dir_name.split("/")[-1])
        mvir, snap = lib.read_tree(dir_name, ["Mvir", "Snap"])
        m_idx, m = lib.read_mergers(dir_name)
        b = lib.read_branches(dir_name)

        survives = snap[b["start"]] == N_SNAP - 1
        is_mw_sub = np.where(b["is_main_sub"] & b["is_real"] &
                             survives & (~b["is_disappear"]))[0]

        lmc_idx_cosmo, lmc_idx_zoom, gse_idx = lib.read_merger_idxs(dir_name)
        print(lmc_idx_cosmo, lmc_idx_zoom, gse_idx)

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
                if b["preprocess"][i] == lmc_idx_cosmo:
                    lmc_cosmo_mah.append(mah)
                if b["preprocess"][i] == lmc_idx_zoom:
                    lmc_zoom_mah.append(mah)
                if b["preprocess"][i] == gse_idx:
                    gse_mah.append(mah)
                    
                    
    npre_mah = np.array(npre_mah)
    pre_mah = np.array(pre_mah)
    lmc_cosmo_mah = np.array(lmc_cosmo_mah)
    lmc_zoom_mah = np.array(lmc_zoom_mah)
    gse_mah = np.array(gse_mah)

    plot_mah(npre_mah, pc("k"), label=r"${\rm direct\ infall}$")
    plot_mah(pre_mah, pc("r"), label=r"${\rm preprocessed}$")
    plot_mah(gse_mah, pc("o"), label=r"${\rm GSE\ subs}$")
    plot_mah(lmc_cosmo_mah, pc("b"), label=r"${\rm LMC\ subs\ (cosmo)}$")
    plot_mah(lmc_zoom_mah, pc("p"), label=r"${\rm LMC\ subs\ (zoom)}$")

    if LOG_SCALING:
        plt.xlabel(r"$\log_{10}(a)$")
        if NORM_BY_MPEAK:
            plt.ylabel(r"$\log_{10}(M_{\rm vir}/M_{\rm peak})$")
        else:
            plt.ylabel(r"$\log_{10}(M_{\rm vir}\,(h^{-1}M_\odot))$")
        plt.ylim(-2, 0)

        plt.legend(loc="lower right", fontsize=16)
    else:
        plt.xlabel(r"$1+z$")
        if NORM_BY_MPEAK:
            plt.ylabel(r"$M_{\rm vir}/M_{\rm peak}$")
        else:
            plt.ylabel(r"$(M_{\rm vir}\,(h^{-1}M_\odot)$")
        plt.xlim(20, 1)
        plt.xscale("log")


        plt.legend(loc="upper left", fontsize=16)
        
    plt.savefig("../plots/preprocess_mah.png")

def plot_mah(mah, c, label=None):
    a = lib.scale_factors()
    med = masked_percentile(mah, -1, 0.5)
    ok = med != -1
    age = cosmo.age(1/a - 1)
    lookback_t = cosmo.age(0) - age
    if LOG_SCALING:
        plt.plot(np.log10(a[ok]), np.log10(med[ok]), c=c, label=label)
    else:
        #plt.plot(lookback_t[ok], med[ok], c=c, label=label)
        plt.plot((1/a)[ok], med[ok], c=c, label=label)
    
    low, high = bootstrap_masked_percentile(mah, -1, 0.5)
    ok = (low != -1) & (high != -1)
    if LOG_SCALING:
        plt.fill_between(np.log10(a[ok]), np.log10(low[ok]), np.log10(high[ok]),
                         alpha=0.2, color=c)
    else:
        #plt.fill_between(lookback_t[ok], low[ok], high[ok],
        #                 alpha=0.2, color=c)
        plt.fill_between((1/a)[ok], low[ok], high[ok],
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
