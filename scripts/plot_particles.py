import numpy as np
import matplotlib.pyplot as plt
try:
    import palette
    from palette import pc
    palette.configure(False)
except:
    pc = lambda x: x
import array
import struct
import scipy.stats as stats
import matplotlib as mpl
import sys
import lib
sys.path.append("/home/users/phil1/code/src/github.com/phil-mansfield/read_gadget")
import read_gadget
import sys

halo = sys.argv[1]
start_snap = 0
base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo%s"%halo
p_file_fmt = "%s/particles/part_%%03d.%%d" % base_dir
id_file_fmt = "%s/particles/ids.%%d" % base_dir
mergers_file = "%s/mergers.dat" % base_dir
mp = 2.81981e5

def tracks_to_halo(tracks):
    h = np.zeros(
        len(tracks[0]), dtype=[("x", "f4", 3), ("v", "f4", 3), ("rvir", "f4")]
    )
    h["x"] = np.array([tracks[0], tracks[1], tracks[2]]).T
    h["v"] = np.array([tracks[3], tracks[4], tracks[5]]).T
    h["rvir"] = tracks[6]
    return h

def plot_circle(ax, x, y, r, lw=1.5, c="k"):
    theta = np.linspace(0, 2*np.pi, 200)
    xc = r*np.sin(theta)
    yc = r*np.cos(theta)

    ax.plot(xc+x, yc+y, lw=lw, c=c)

def read_part_file(snap, i):
    fname = p_file_fmt % (snap, i)
    f = open(fname, "rb")

    n = struct.unpack("i", f.read(4))[0]
    x, v, phi = array.array("f"), array.array("f"), array.array("f")

    x.fromfile(f, n*3)
    x = np.array(x, dtype=float)
    v.fromfile(f, n*3)
    v = np.array(v, dtype=float)
    phi.fromfile(f, n)
    phi = np.array(phi, dtype=float)
    x, v = x.reshape((n, 3)), v.reshape((n, 3))

    f.close()
    return x, v, phi

def phi_ranks(h, halo_idx):
    if halo_idx == 0:
        i_peak = len(h) - 1
    else:
        i_peak = np.argmax(h["rvir"])

    x, v, phi = read_part_file(i_peak, halo_idx)
    ok = (x[:,0] >= 0) & (phi != -1)

    ke = kinetic_energy(h["v"][i_peak], v)

    rank = np.ones(len(ok))*-1
    rank[ok] = stats.rankdata(phi[ok] + ke[ok]) / np.sum(ok)
    return rank

def bound_center(rank, x, frac=1e-4, vec=None):
    if vec is None: vec = x
    ok = (rank > 0) & (rank < frac)
    return np.median(vec[ok,:], axis=0)

def plot_halo_radius(fig, ax, i, mw_h, lmc_h, ges_h,
                     mw_rank, lmc_rank, ges_rank):
    ((ax_xy1, ax_xy2, ax_xy3), (ax_xz1, ax_xz2, ax_xz3)) = ax
    rmax = 2.5*mw_h["rvir"][-1]
    scale = 10**np.linspace(np.log10(0.05), 0, 236)

    ax_xy1.set_xlim(-rmax, rmax)
    ax_xz1.set_xlim(-rmax, rmax)
    ax_xy2.set_xlim(-rmax, rmax)
    ax_xz2.set_xlim(-rmax, rmax)
    ax_xy3.set_xlim(-rmax, rmax)
    ax_xz3.set_xlim(-rmax, rmax)

    ax_xy1.set_ylim(-rmax, rmax)
    ax_xz1.set_ylim(-rmax, rmax)
    ax_xy2.set_ylim(-rmax, rmax)
    ax_xz2.set_ylim(-rmax, rmax)
    ax_xy3.set_ylim(-rmax, rmax)
    ax_xz3.set_ylim(-rmax, rmax)

    plot_circle(ax_xy1, 0, 0, mw_h["rvir"][i], c="k")
    plot_circle(ax_xz1, 0, 0, mw_h["rvir"][i], c="k")
    plot_circle(ax_xy2, 0, 0, mw_h["rvir"][i], c="k")
    plot_circle(ax_xz2, 0, 0, mw_h["rvir"][i], c="k")
    plot_circle(ax_xy3, 0, 0, mw_h["rvir"][i], c="k")
    plot_circle(ax_xz3, 0, 0, mw_h["rvir"][i], c="k")

    x_mw, _, phi_mw = read_part_file(i, 0)
    x_lmc, _, phi_lmc = read_part_file(i, 1)
    x_ges, _, phi_ges = read_part_file(i, 2)
    ok_mw = (x_mw[:,0] >= 0) & (mw_rank >= 0)
    ok_lmc = (x_lmc[:,0] >= 0) & (lmc_rank >= 0)
    ok_ges = (x_ges[:,0] >= 0) & (ges_rank >= 0)

    for dim in range(3):
        x_mw[:,dim] -= mw_h["x"][i,dim]
        x_lmc[:,dim] -= mw_h["x"][i,dim]
        x_ges[:,dim] -= mw_h["x"][i,dim]

    x_mw, x_lmc, x_ges = x_mw[ok_mw], x_lmc[ok_lmc], x_ges[ok_ges]
    mw_rank_i = mw_rank[ok_mw]
    lmc_rank_i, ges_rank_i = lmc_rank[ok_lmc], ges_rank[ok_ges]
    mw_core = mw_rank_i < 1e-4
    lmc_core, ges_core = lmc_rank_i < 1e-4, ges_rank_i < 1e-4

    if mw_h["rvir"][i] > 0:
        plot_circle(ax_xy3, mw_h["x"][i,0], mw_h["x"][i,1],
                    lmc_h["rvir"][i], c=pc("r"))
        plot_circle(ax_xz3, mw_h["x"][i,0], mw_h["x"][i,2],
                    lmc_h["rvir"][i], c=pc("r"))
    if lmc_h["rvir"][i] > 0:
        plot_circle(ax_xy1, lmc_h["x"][i,0], lmc_h["x"][i,1],
                    lmc_h["rvir"][i], c=pc("r"))
        plot_circle(ax_xz1, lmc_h["x"][i,0], lmc_h["x"][i,2],
                    lmc_h["rvir"][i], c=pc("r"))
    if ges_h["rvir"][i] > 0:
        plot_circle(ax_xy2, ges_h["x"][i,0], ges_h["x"][i,1],
                    ges_h["rvir"][i], c=pc("b"))
        plot_circle(ax_xz2, ges_h["x"][i,0], ges_h["x"][i,2],
                    ges_h["rvir"][i], c=pc("b"))

    extent = [-rmax, rmax, -rmax, rmax]
    gridsize = 100
    vmin, vmax = 5, 1e4

    ax_xy1.hexbin(x_lmc[:,0], x_lmc[:,1], gridsize=gridsize, cmap="Greys",
                  bins="log", extent=extent, vmin=vmin, vmax=vmax)
    ax_xy1.plot(x_lmc[lmc_core,0], x_lmc[lmc_core,1],
                ".", c=pc("r"), alpha=0.15)
    ax_xz1.hexbin(x_lmc[:,0], x_lmc[:,2], gridsize=gridsize, cmap="Greys",
                  bins="log", extent=extent, vmin=vmin, vmax=vmax)
    ax_xz1.plot(x_lmc[lmc_core,0], x_lmc[lmc_core,2],
                ".", c=pc("r"), alpha=0.15)
    ax_xy2.hexbin(x_ges[:,0], x_ges[:,1], gridsize=gridsize, cmap="Greys",
                  bins="log", extent=extent, vmin=vmin, vmax=vmax)
    ax_xy2.plot(x_ges[ges_core,0], x_ges[ges_core,1],
                ".", c=pc("o"), alpha=0.15)
    ax_xz2.hexbin(x_ges[:,0], x_ges[:,2], gridsize=gridsize, cmap="Greys",
                  bins="log", extent=extent, vmin=vmin, vmax=vmax)
    ax_xz2.plot(x_ges[ges_core,0], x_ges[ges_core,2],
                ".", c=pc("o"), alpha=0.15)
    ax_xy3.hexbin(x_mw[:,0], x_mw[:,1], gridsize=gridsize, cmap="Greys",
                  bins="log", extent=extent, vmin=vmin, vmax=vmax)
    ax_xy3.plot(x_mw[mw_core,0], x_mw[mw_core,1],
                ".", c=pc("b"), alpha=0.15)
    ax_xz3.hexbin(x_mw[:,0], x_mw[:,2], gridsize=gridsize, cmap="Greys",
                  bins="log", extent=extent, vmin=vmin, vmax=vmax)
    ax_xz3.plot(x_mw[mw_core,0], x_mw[mw_core,2],
                ".", c=pc("b"), alpha=0.15)

    ax_xy1.set_ylabel(r"$Y\,(h^{-1}{\rm cMpc})$")
    ax_xz1.set_ylabel(r"$Z\,(h^{-1}{\rm cMpc})$")
    ax_xz1.set_ylabel(r"$X\,(h^{-1}{\rm cMpc})$")
    ax_xz2.set_xlabel(r"$X\,(h^{-1}{\rm cMpc})$")
    ax_xz3.set_xlabel(r"$X\,(h^{-1}{\rm cMpc})$")

    ax_xy2.set_title(r"$a(z) = %.3f$" % scale[i])

    fig.savefig("../plots/Halo%s//halo_radius.%d.png" % (halo, i))

    ax_xy1.clear()
    ax_xz1.clear()
    ax_xy2.clear()
    ax_xz2.clear()
    ax_xy3.clear()
    ax_xz3.clear()


def halo_center(i, h, rank, x, v):
    i_peak = np.argmax(h["rvir"])
    if i <= i_peak:
        return h["x"][i], h["v"][i]
    else:
        return bound_center(rank, x), bound_center(rank, x, dv, vec=v)

def kinetic_energy(v0, v):
    code_units = 149.4
    dv = np.zeros(v.shape)
    for dim in range(3):
        dv[:,dim] = v[:,dim] - v0[dim]
    
    return np.sum((dv/code_units)**2, axis=1)/2

def plot_rho(ax, c, x, rmax, scale):
    bins, log_rmin, log_rmax = 50, np.log10(rmax/1e3), np.log10(rmax)
    r = np.sqrt(np.sum((x*scale)**2, axis=1))
    ok = r > 0
    x, r = x[ok], r[ok]

    log_r = np.log10(r)
    n, edges = np.histogram(log_r, range=(log_rmin, log_rmax), bins=bins)

    mass = n*mp
    vol = 4*np.pi/3 * (10**edges)**3
    dvol = vol[1:] - vol[:-1]
    mid = (edges[1:] + edges[:-1]) / 2
    r = 10**mid/0.7
    ok = n > 5

    rho = np.log10(mass[ok] / dvol[ok])

    ax.plot(mid[ok], rho, c=c)
    ax.set_ylim(11, 18.5)

def plot_vrms(ax, c, x, v, rmax, scale):
    bins, log_rmin, log_rmax = 25, np.log10(rmax/1e3), np.log10(rmax)
    r = np.sqrt(np.sum((x*scale)**2, axis=1))
    ok = r != 0
    x, v, r = x[ok], v[ok], r[ok]

    vtot = np.sqrt(np.sum(v**2, axis=1))
    log_r = np.log10(r)
    
    vr = np.sum(x*v, axis=1)/r
    vrms, edges, _ = stats.binned_statistic(
        log_r, vr*vr, "mean", range=(log_rmin, log_rmax), bins=bins
    )
    n, _ = np.histogram(log_r, range=(log_rmin, log_rmax), bins=bins)
    vrms = np.sqrt(vrms)

    mid = (edges[1:] + edges[:-1]) / 2
    r = 10**mid/0.7
    ok = n > 5

    ax.plot(mid[ok], vrms[ok], c=c)
    ax.set_ylim(0, 400)

#def get_mw_particles(x0, r0, i):
#    for b in range(8):
#        file_name = g2_file_fmt % (i, b)
#        f = read_gadget.Gadget2Zoom(file_name, ["x", "v", "id32"])
#        
#    return None, None

def plot_profiles(fig, ax, i, mw_h, lmc_h, ges_h, lmc_rank, ges_rank):
    ((ax_rho1, ax_rho2), (ax_vrms1, ax_vrms2)) = ax
    rmax = mw_h["rvir"][-1]
    scale = 10**np.linspace(np.log10(0.05), 0, 236)

    x_mw, v_mw, phi_mw = read_part_file(i, 0)
    x_lmc, v_lmc, phi_lmc = read_part_file(i, 1)
    x_ges, v_ges, phi_ges = read_part_file(i, 2)
    ok_lmc = (x_lmc[:,0] >= 0) & (lmc_rank >= 0)
    ok_ges = (x_ges[:,0] >= 0) & (ges_rank >= 0)

    v_lmc *= np.sqrt(scale[i])
    v_ges *= np.sqrt(scale[i])

    for dim in range(3):
        x_lmc[:,dim] -= mw_h["x"][i,dim]
        x_ges[:,dim] -= mw_h["x"][i,dim]
        v_lmc[:,dim] -= mw_h["v"][i,dim]
        v_ges[:,dim] -= mw_h["v"][i,dim]

    x_lmc, x_ges = x_lmc[ok_lmc], x_ges[ok_ges]
    v_lmc, v_ges = v_lmc[ok_lmc], v_ges[ok_ges]
    lmc_rank, ges_rank = lmc_rank[ok_lmc], ges_rank[ok_ges]

    x0_lmc, v0_lmc = halo_center(i, lmc_h, lmc_rank, x_lmc, v_lmc)
    x0_ges, v0_ges = halo_center(i, ges_h, ges_rank, x_ges, v_ges)

    dx_lmc, dv_lmc = np.zeros(x_lmc.shape), np.zeros(v_lmc.shape)
    dx_ges, dv_ges = np.zeros(x_ges.shape), np.zeros(v_ges.shape)
    for dim in range(3):
        dx_lmc[:,dim] = x_lmc[:,dim] - x0_lmc[dim]
        dx_ges[:,dim] = x_ges[:,dim] - x0_ges[dim]
        dv_lmc[:,dim] = v_lmc[:,dim] - v0_lmc[dim]
        dv_ges[:,dim] = v_ges[:,dim] - v0_ges[dim]

    r_search = min(mw_h["rvir"][i]*3, rmax/scale[i])
    x_mw, v_mw = get_mw_particles(mw_h["x"][i], r_search, i)

    plot_rho(ax_rho1, pc("r"), x_lmc, rmax, scale[i])
    plot_rho(ax_rho1, pc("o"), x_ges, rmax, scale[i])
    plot_rho(ax_rho2, pc("r"), dx_lmc, rmax, scale[i])
    plot_rho(ax_rho2, pc("o"), dx_ges, rmax, scale[i])
    plot_vrms(ax_vrms1, pc("r"), x_lmc, v_lmc, rmax, scale[i])
    plot_vrms(ax_vrms1, pc("o"), x_ges, v_ges, rmax, scale[i])
    plot_vrms(ax_vrms2, pc("r"), dx_lmc, dv_lmc, rmax, scale[i])
    plot_vrms(ax_vrms2, pc("o"), dx_ges, dv_ges, rmax, scale[i])
    
    ax_rho1.set_title(r"$a(z) = %.3f$" % scale[i])
    ax_rho1.set_ylabel(r"$\log_{10}(\rho)\,(M_\odot/{\rm pMpc})$")
    ax_vrms1.set_ylabel(r"$v_{\rm rms}\,({\rm km\,s^{-1}})$")
    ax_vrms1.set_xlabel(r"$\log_{10}(r_{\rm MW})\,({\rm pMpc})$")
    ax_vrms2.set_xlabel(r"$\log_{10}(r_{\rm sub})\,({\rm pMpc})$")

    fig.savefig("../plots/Halo%s/profiles.%d.png" % (halo, i))
    ax_rho1.clear()
    ax_rho2.clear()
    ax_vrms1.clear()
    ax_vrms2.clear()

def main():
    idx, mergers = lib.read_mergers(base_dir)
    print(idx)
    mw_h, lmc_h, ges_h = mergers[0], mergers[1], mergers[2]
    #tracks = np.loadtxt(tracks_file).T
    #mw_h = tracks_to_halo(tracks[0:7])
    #lmc_h = tracks_to_halo(tracks[7:14])
    #ges_h = tracks_to_halo(tracks[14:21])

    mw_rank = phi_ranks(mw_h, 0)
    lmc_rank = phi_ranks(lmc_h, 1)
    ges_rank = phi_ranks(ges_h, 2)

    pix = 1024
    dpi = pix/16
    fig_part, ax_part = plt.subplots(2, 3, sharex=True, sharey=False,
                                     figsize=(3*pix/2/dpi, pix/dpi), dpi=dpi)
    fig_prof, ax_prof= plt.subplots(2, 2, sharex=True, sharey=False,
                                    figsize=(pix/dpi, pix/dpi), dpi=dpi)

    ges_h["x"] -= mw_h["x"]
    lmc_h["x"] -= mw_h["x"]
    ges_h["v"] -= mw_h["v"]
    lmc_h["v"] -= mw_h["v"]

    for i in range(start_snap, 236):
        print("snap", i)
        plot_halo_radius(
            fig_part, ax_part, i, mw_h, lmc_h, ges_h,
            mw_rank, lmc_rank, ges_rank
        )
        #plot_profiles(
        #    fig_part, ax_part, i, mw_h, lmc_h, ges_h, lmc_rank, ges_rank
        #)

if __name__ == "__main__": main()
