import symlib
import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import os.path as path

def plot_mass_loss():
    palette.configure(False)

    out_dir = "../plots/core_tests/mass_loss"
    base_dir = "../tmp_data"
    suite = "SymphonyMilkyWay"
    sim_dir = symlib.get_host_directory(base_dir, suite, 0)

    param = symlib.simulation_parameters(suite)
    a = symlib.scale_factors(sim_dir)
    
    h, hist = symlib.read_subhalos(param, sim_dir)
    h = symlib.set_units_halos(h, a, param)
    c = symlib.read_cores(sim_dir)

    targets = np.arange(1, len(h), dtype=int)

    fig, ax = plt.subplots(2, figsize=(7, 14), sharex=True)
    
    for i_sub in targets:
        print(i_sub)
        ax[0].cla()
        ax[1].cla()

        ax[0].plot(a, h["mvir"][0], "--", c="k")
        
        okh, mvir = h["ok"][i_sub], h["mvir"][i_sub]
        rh = np.sqrt(np.sum(h["x"][i_sub]**2, axis=1))
    
        ax[0].plot(a[okh], mvir[okh], pc("r"), label=r"$M_{\rm vir,RS}$")
        ax[1].plot(a[okh], rh[okh], pc("r"), label=r"$r_{\rm Rockstar}$")
        okc = c["ok"][i_sub]
        m_bound, m_tidal = c["m_bound"][i_sub], c["m_tidal_bound"][i_sub]
        r_half = c["r50_bound"][i_sub]
        rc = np.sqrt(np.sum(c["x"][i_sub]**2, axis=1))
        ax[0].plot(a[okc], m_bound[okc], pc("o"), label=r"$M_{\rm bound}$")
        ax[0].plot(a[okc], m_tidal[okc], pc("b"), label=r"$M_{\rm tidal}$")
        ax[1].plot(a[okc], rc[okc], pc("b"), label=r"$r_{\rm core}$")
        ax[1].plot(a[okh], rh[okh], ":", c=pc("r"))
        ax[1].plot(a, h[0]["rvir"], "--", c="k")
        ax[1].plot(a[okc], r_half[okc], "--", c=pc("a"))

        mp = param["mp"]/param["h100"]
        ax[0].set_ylim(mp, None)
        ax[0].fill_between([0, 1], [mp]*2, [25*mp]*2, color="k", alpha=0.2)
        ax[0].set_yscale("log")
        ax[1].set_yscale("log")
        ax[1].set_xlim(0, 1)
        ax[1].set_ylim(np.min(h[0]["rvir"]), None)
        
        ax[1].set_xlabel(r"$a(z)$")
        ax[1].set_ylabel(r"$r\ ({\rm kpc})$")
        ax[0].set_ylabel(r"$M\ (M_\odot)$")

        ax[0].legend(loc="upper right", fontsize=16)
        ax[1].legend(loc="lower left", fontsize=16)

        fig.savefig(path.join(out_dir, "ml_%04d.png" % i_sub))

def basic_halo_stats():
    palette.configure(False)

    out_dir = "../plots/core_tests/mass_loss"
    base_dir = "../tmp_data"
    suite = "SymphonyMilkyWay"
    sim_dir = symlib.get_host_directory(base_dir, suite, 0)

    param = symlib.simulation_parameters(suite)
    a = symlib.scale_factors(sim_dir)
    
    h, hist = symlib.read_subhalos(param, sim_dir)
    h = symlib.set_units_halos(h, a, param)
    c = symlib.read_cores(sim_dir)

    mp = param["mp"]/param["h100"]

    m_vir = h["mvir"][1:,-1]
    m_peak = np.max(h["mvir"], axis=1)[1:]
    r_rs = np.sqrt(np.sum(h["x"][:,-1]**2, axis=1))
    okh = h["ok"][1:,-1] & (m_vir > 25*mp)
    
    r_c = np.sqrt(np.sum(c["x"][1:,-1]**2, axis=1))
    r_half = c["r50_bound"][1:,-1]
    m_bound = c["m_bound"][1:,-1]
    m_tidal = c["m_tidal"][1:,-1]
    okc = c["ok"][1:,-1] & (m_bound > 25*mp) & (r_c > r_half)

    m0 = h["mvir"][0,-1]
    
    fig, ax = plt.subplots()
    ax.plot(np.sort(m_peak[okh])/m0,
            np.arange(len(m_vir[okh]))[::-1],
            c=pc("r"))
    ax.plot(np.sort(m_peak[okc])/m0,
            np.arange(len(m_vir[okc]))[::-1],
            c=pc("b"))
    
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(r"$m=M_{\rm peak,sub}/M_{\rm host}$")
    ax.set_ylabel(r"$N(>m)$")
    
    plt.show()
        
def main():
    plot_mass_loss()
    #basic_halo_stats()
    
if __name__ == "__main__": main()
