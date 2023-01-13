import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import symlib
import os.path as path
from colossus.cosmology import cosmology
from colossus.halo import mass_so
import scipy.signal as signal

TEST_SMOOTHING = False

suite = "SymphonyMilkyWay"
base_dir = "/sdf/home/p/phil1/ZoomIns"
plot_dir_fmt = "/scratch/phil1/plots/tracking/mah/h%d"
i_host = 1

def split_plot(ax, x, y, ok, ls=None, c=None, lw=None):
    idx = np.arange(len(ok), dtype=int)[ok]
    
    mids = np.where(idx[:-1]+1 != idx[1:])[0] + 1
    edges = np.zeros(len(mids) + 2, dtype=int)
    edges[1:-1], edges[-1] = mids, len(idx)
    start, end = edges[:-1], edges[1:]
    for i in range(len(start)):
        if start[i] + 1 == end[i]:
            j = start[i]
            ax.plot(x[idx[j]], y[idx[j]], ".", c=c)
        else:
            ax.plot(x[idx[start[i]: end[i]]], y[idx[start[i]: end[i]]],
                    ls=ls, c=c, lw=lw)


def main():
    for i_host in range(1):
        plot_dir = plot_dir_fmt % i_host
        host_dir = symlib.get_host_directory(base_dir, suite, i_host)
        param = symlib.simulation_parameters(host_dir)
        cosmo = cosmology.setCosmology('', symlib.colossus_parameters(param))
        h, hist = symlib.read_subhalos(host_dir)
        c = symlib.read_cores(host_dir)
        scale = symlib.scale_factors(host_dir)
        z = 1/scale - 1
        t = cosmo.age(z)
        
        mp = param["mp"]/param["h100"]
        fig, (ax2, ax1) = plt.subplots(2, figsize=(7, 14), sharex=True)

        smooth_order = 4
        smooth_windows = [5, 11, 17]
        smooth_colors = [pc("p", 0.2), pc("p", 0.5), pc("p", 0.8)]

        #for i in range(1, len(h)):
        for i in range(1, len(h)):
            if i != 14: continue
            print("%d/%d" % (i, len(h)-1))
            ax1.cla()
            ax2.cla()
            #n_ax1 = ax1.twinx()

            R_host, m_host, ok_host = h["rvir"][0], h["mvir"][0], h["ok"][0]
            m_rs, m_c = h["mvir"][i], c["m_bound"][i]
            ok1_rs, ok2_rs, ok_c = h["ok"][i], c["ok_rs"][i], c["ok"][i]
            r_rs = np.sqrt(np.sum(h["x"][i]**2, axis=1))
            r_c = np.sqrt(np.sum(c["x"][i]**2, axis=1))
            
            ax1.plot(t[ok_c], m_c[ok_c], c=pc("b"), lw=3.5)
            ax1.plot(t[ok_host], m_host[ok_host], "--", c="k")
            ax1.plot(t[ok1_rs], m_rs[ok1_rs], c=pc("o"), lw=2.5)
            split_plot(ax1, t, m_rs, ok2_rs, ls="-", c=pc("r"), lw=2.5)

            ax2.plot(t[ok_c], r_c[ok_c], c=pc("b"), lw=3.5,
                     label=r"${\rm particle{-}tracking}$")
            ax2.plot(t[ok_host], R_host[ok_host], "--", c="k")
            ax2.plot(t[ok1_rs], r_rs[ok1_rs], c=pc("o"), lw=2.5)
            split_plot(ax2, t, r_rs, ok2_rs, ls="-", c=pc("r"), lw=2.5)
            ax2.plot([], [], c=pc("r"), label=r"$\textsc{Rockstar}$", lw=2.5)
            ax2.plot([], [], c=pc("o"), label=r"$\textsc{Rockstar}\ ({\rm error})$", lw=2.5)

            for j in range(len(smooth_windows)):
                if not TEST_SMOOTHING: break
                if np.sum(ok_c) < smooth_windows[j]: continue
                m_smooth = 10**signal.savgol_filter(
                    np.log10(m_c[ok_c]), smooth_windows[j], smooth_order)
                ax1.plot(t[ok_c], m_smooth, lw=1.5, c=smooth_colors[j])

            #ylo, yhi = ax1.get_ylim()
            #n_ax1.set_ylim(ylo/mp, yhi/mp)
            #n_ax1.set_yscale("log")
            ax1.set_yscale("log")
            ax2.set_yscale("log")
            ax2.legend(fontsize=16)
            ax1.set_ylabel(r"$m_{\rm sub}\ (M_\odot)$")
            #n_ax1.set_ylabel(r"$N_{p,{\rm sub}}$")
            ax2.set_ylabel(r"$r_{\rm sub}\ ({\rm pkpc})$")
            
            #ax1.set_xlabel(r"$a(t)$")
            ax1.set_xlabel(r"$t\ ({\rm Gyr})$")
            ok_host = h["ok"][0,:]
            ax1.set_xlim((np.min(t[ok_host]), np.max(t[ok_host])))

            fig.savefig(path.join(plot_dir, "sub_%03d.png" % i))

if __name__ == "__main__":
    palette.configure(True)
    main()
