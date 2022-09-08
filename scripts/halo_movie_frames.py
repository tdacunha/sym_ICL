import numpy as np
import matplotlib.pyplot as plt
import lib
import star_tagging
import os
import sys
import os.path as path
import matplotlib.colors as mpl_colors
import symlib
import scipy.spatial as spatial
import scipy.special as special
from colossus.halo import mass_so

try:
    import palette
    palette.configure(False)
    from palette import pc
except:
    def pc(c): return c

SUITE = "SymphonyMilkyWayLR"
SUB_INDEX = 33
HOST_INDEX = 4

BASE_OUT_DIR = "/home/users/phil1/code/src/github.com/phil-mansfield/symphony_pipeline/plots/movie_frames"

def spline(x):
    r1 = x <= 1
    r2 = (x > 1) & (x <= 2)
    
    out = np.zeros(len(x))
    out[r1] = 1 - (3/2)*x[r1]**2 + 0.75*x[r1]**3
    out[r2] = 0.25*(2 - x[r2])**3

    return out
    
def rho_sph(ri, mi, n_dim):
    r = np.max(ri)/2
    Vr = np.pi**(n_dim/2)*r**n_dim / special.gamma(n_dim/2 + 1)
    # No need to add the central particle's mass: It's already in the sample
    # if it's a real particle and its mass is zero if it's a test point.
    return np.sum(mi*spline(ri/r))/Vr

def sph_density(x, m, test, k=64, return_tree=False):
    k = min(k, len(x))
    tree = spatial.cKDTree(x)
    
    rho = np.zeros(len(test))
    for i in range(len(test)):
        ri, idx = tree.query(test[i], k)
        rho[i] = rho_sph(ri, m[idx], len(x[0]))

    if return_tree:
        return rho, tree
    else:
        return rho

def eval_grid(r, center, pts):
    x = np.linspace(center[0] - r, center[0] + r, pts)
    y = np.linspace(center[1] - r, center[1] + r, pts)
    xy, yx = np.meshgrid(x, y)
    xy, yx = xy.flatten(), yx.flatten()
    return np.vstack((xy, yx)).T

def density_2d_proj(x, m, test, dim_x=0, dim_y=1, k=64, return_tree=False):
    xx = np.zeros((len(x), 2))
    xx[:,0], xx[:,1] = x[:,dim_x], x[:,dim_y]
    if test.shape[1] == 2:
        yy = test
    else:
        yy = np.zeros((len(test), 2))
        yy[:,0], yy[:,1] = test[:,dim_x], test[:,dim_y]

    return sph_density(xx, m, yy, k=k)

def is_bound(param, dx, dv, ok=None, order=None):
    rmax, vmax, pe, order = symlib.profile_info(param, dx, ok, order)
    pe *= vmax**2
    ke = 0.5*np.sum(dv**2, axis=1)

    if ok is None:
        return pe + ke < 0, order
    else:
        return (pe + ke < 0) & ok, order

def is_bound_iter(n_iter, param, dx, dv, ok=None, order=None):
    ok, order = is_bound(param, dx, dv, ok, order)
    n_bound = np.sum(ok)

    for i in range(n_iter - 1):
        ok, order = is_bound(param, dx, dv, ok, order)
        n_bound_i = np.sum(ok)
        if n_bound_i == n_bound: break

    return ok, order

def main():
    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
        
    # Work out subdirectory names and create as needed.
    print(HOST_INDEX)
    out_dir = path.join(BASE_OUT_DIR, "host_%d" % HOST_INDEX,
                        "sub_%d" % SUB_INDEX)
    print(out_dir)
    os.makedirs(out_dir, exist_ok=True)    
    sim_dir = symlib.get_host_directory(base_dir, SUITE, HOST_INDEX)
    i_sub = SUB_INDEX

    # Load simulation information.
    param = symlib.simulation_parameters(sim_dir)
    cosmo = symlib.colossus_parameters(param)
    h100 = param["h100"]
    mp, eps = param["mp"]/h100, param["eps"]/h100
    h, hist = symlib.read_subhalos(sim_dir)
    h_cmov, hist_cmov = symlib.read_subhalos(sim_dir, comoving=True)
    c = symlib.read_cores(sim_dir)
    scale = symlib.scale_factors(sim_dir)
    snaps = np.arange(len(scale), dtype=int)    

    r_c = np.sqrt(np.sum(c["x"]**2, axis=2))
    r_half = c["r50_bound"]
    c["ok"] = c["ok"] & (r_c > c["r50_bound"]) & (c["m_bound"] > 32*mp)

    frame_idx = 0
    # fig, ax = plt.subplots(1, 2, figsize=(16, 8))
    fig, axs = plt.subplots(1, 2, figsize=(16, 8))
    ax, ax_b = axs[0], axs[1]
    already_found = False
    c_range = None
    
    info = symlib.ParticleInfo(sim_dir)

    cores = symlib.read_particles(info, sim_dir, None, "infall_core")

    for snap in range(len(snaps)):
        #if snap < len(snaps) - 10: continue
        if snap < hist["first_infall_snap"][i_sub]: continue
        print(snap)

        x = symlib.read_particles(info, sim_dir, snap, "x", owner=i_sub)
        v = symlib.read_particles(info, sim_dir, snap, "v", owner=i_sub)
        x_core = x[cores[i_sub]]
        ok = symlib.read_particles(info, sim_dir, snap, "valid", owner=i_sub)

        x = symlib.set_units_x(x[ok], h_cmov[0,snap], scale[snap], param)
        v = symlib.set_units_v(v[ok], h_cmov[0,snap], scale[snap], param)
        x_core = symlib.set_units_x(x_core, h_cmov[0,snap], scale[snap], param)
        mp_dm = np.ones(len(x)) * mp

        if c[i_sub,snap]["ok"]:
            dx = x - c[i_sub,snap]["x"]
            dv = v - c[i_sub,snap]["v"]
            ok_b, order = is_bound_iter(10, param, dx, dv, ok=ok)

        # Everything from here on is just plotting nonsense.
        r_max, grid_pts = np.max(h[0,:]["rvir"]), 201
        r_max_b = np.max(c[i_sub,:]["r50_bound"])
        grid = eval_grid(r_max, [0, 0], grid_pts)
        grid_b = eval_grid(r_max_b, [0, 0], grid_pts)

        vmin, vmax = 2.5, 6.5

        rho_proj = density_2d_proj(x[ok], mp_dm[ok], grid, k=64)
        rho_proj = np.log10(np.reshape(rho_proj, (grid_pts, grid_pts)))
        if c[i_sub,snap]["ok"] and np.sum(ok_b) > 1:
            rho_proj_b = density_2d_proj(dx[ok_b], mp_dm[ok_b], grid_b, k=64)
            rho_proj_b = np.log10(np.reshape(rho_proj_b, (grid_pts, grid_pts)))
        else:
            rho_proj_b = vmin * np.ones((grid_pts, grid_pts))

        ax.cla()
        ax_b.cla()

        ax.imshow(rho_proj, vmin=vmin, vmax=vmax, origin="lower",
                  cmap="bone", extent=[-r_max, r_max, -r_max, r_max])
        ax.set_xlim((-r_max, r_max))
        ax.set_ylim((-r_max, r_max))
        ax_b.imshow(rho_proj_b, vmin=vmin, vmax=vmax, origin="lower",
                    cmap="bone", extent=[-r_max_b, r_max_b, -r_max_b, r_max_b])
        ax_b.set_xlim((-r_max_b, r_max_b))
        ax_b.set_ylim((-r_max_b, r_max_b))

        if c["ok"][i_sub,snap]:
            symlib.plot_circle(ax, c["x"][i_sub,snap,0], c["x"][i_sub,snap,1],
                               c["r50_bound"][i_sub,snap], lw=1.5, c=pc("r"))
            symlib.plot_circle(ax_b, 0, 0, c["r50_bound"][i_sub,snap],
                               lw=1.5, c=pc("r"))
            dx_core = x_core - c["x"][i_sub,snap]
            ax_b.plot(dx_core[:10,0], dx_core[:10,1], ".", alpha=0.5, c=pc("r"))
            if h["ok"][i_sub,snap]:
                dx_h = h["x"][i_sub,snap] - c["x"][i_sub,snap]
                symlib.plot_circle(ax_b, dx_h[0], dx_h[1], h["rvir"][i_sub,snap],
                                   lw=1.5, c=pc("b"))
        if h["ok"][i_sub,snap]:
            symlib.plot_circle(ax, h["x"][i_sub,snap,0], h["x"][i_sub,snap,1],
                               h["rvir"][i_sub,snap], lw=1.5, c=pc("b"))

        ax.plot(x_core[:10,0], x_core[:10,1], ".", alpha=0.5, c=pc("r"))
        symlib.plot_circle(ax, 0, 0, h["rvir"][0,snap], lw=3, c="w")

        fig.savefig(path.join(out_dir, "frame_%03d.png" % frame_idx))
        frame_idx += 1

def calc_r_half(dx, mp):
    r = np.sqrt(np.sum(dx**2, axis=1))
    r_min, r_max = max(np.min(r), 0.1), np.max(r)
    bins = np.linspace(np.log10(r_min), np.log10(r_max), 101)
    m, _ = np.histogram(np.log10(r), bins=bins, weights=mp)

    m_enc = np.cumsum(m)
    return 10**bins[np.searchsorted(m_enc, m_enc[-1]/2)+1]

def get_c_range(x, frac):
    x = np.sort(x.flatten())
    x_sum = np.cumsum(x)
    x_tot = x_sum[-1]
    i_min = np.searchsorted(x_sum, x_tot*(1-frac)/2)
    i_max = np.searchsorted(x_sum, x_tot*(1+frac)/2)
    return x[i_min], x[-1]
    
if __name__ == "__main__": main()
