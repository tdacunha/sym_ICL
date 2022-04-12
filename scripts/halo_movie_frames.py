import numpy as np
import matplotlib.pyplot as plt
import lib
import star_tagging
import os
import sys
import os.path as path
import matplotlib.colors as mpl_colors

try:
    import palette
    palette.configure(False)
except:
    pass
    
H100 = 0.7
MP_DM = 2.8e5/0.7 # In Msun
EPS = 170e-3/H100 # In ckpc

def main():
    # Parse command line arguements
    sim_dir = sys.argv[1] # directory containing the HaloXXX directories
    host_id = int(sys.argv[2]) # ID of the host halo
    h_idx = int(sys.argv[3]) # index of the subhalo in mergers.dat
    out_base_dir = sys.argv[4] # base dir where the frames are stored
    # You can choose what snapshot to tag at, if you want.
    if len(sys.argv) > 5:
        tag_snap = int(sys.argv[5])
    else:
        tag_snap = None
        
    snaps, scales = np.arange(236), lib.scale_factors()
        
    # Work out subdirectory names and create as needed.
    base_dir = path.join(sim_dir, "Halo%03d" % host_id)
    out_dir = path.join(path.join(
        out_base_dir, "Halo%03d" % host_id), "h%03d" % h_idx)
    os.makedirs(out_dir, exist_ok=True)

    # *****
    # Read in mergers.dat, which contains main branch information for the
    # biggest mergers in the host's history.
    _, mergers = lib.read_mergers(base_dir)

    # *****
    # If needed, compute the snapshot where the halo first falls into the host.
    if tag_snap is None:
        ok = mergers[h_idx]["mvir"] > 0
        tag_snap = lib.merger_snap(
            mergers[0], mergers[h_idx,ok]["x"], snaps[ok]) - 1
    
    # *****
    # Create the galaxy-halo model. This is done by loading your favorite
    # Mstar-Mhalo relation, Rgalaxy-Rhalo relation, and stellar mass profile
    # model.
    galaxy_halo_model = star_tagging.GalaxyHaloModel(
        star_tagging.UniverseMachineMStar(),
        star_tagging.Jiang2019RHalf(),
        star_tagging.PlummerProfile()
    )

    # *****
    # Find the stellar masses. This function handles all the boilerplate code
    # for you.
    mp_star, ranks = star_tagging.default_tag(
        base_dir, MP_DM, galaxy_halo_model, mergers, h_idx, tag_snap)
    
    # Loop over the snapshots to create frames
    frame_idx = 0
    fig, ax = plt.subplots(1, 2, figsize=(16, 8))
    already_found = False
    for snap in range(len(snaps)):
        # Don't start plotting until the halo actually exists.
        if not already_found and (mergers[h_idx,snap]["mvir"] == -1 or 
                                  snap < tag_snap or
                                  mergers[0,snap]["mvir"] == -1):
            continue
        already_found=True

        # When testing locally, you might not have all the snapshots. but in
        # general, don't just try-and-except things.
        try:
            # *****
            # Read in particles. You can also read in velocities and potential
            # energies by adding "v" and "phi" to the list argument. They'll
            # be the scond and third return values.
            x, _, _ = lib.read_part_file(base_dir, snap, h_idx, ["x"])
        except:
            continue

        print(snap)
        
        # *****
        # Remove as-yet unaccreted particles, convert to physical, center
        # around the host. You could also give velocities as the second argument
        # and get them cleaned as the second return value.
        # idx is the indices of accreted particles into the full particle array
        # (e.g. something with the same length as the mp array)
        x, _, idx = star_tagging.clean_particles(
            x, None, mergers[0,snap], H100, scales[snap])

        # *****
        # Load particles into the ranks so you can calculate things like the
        # core position.
        ranks.load_particles(x, None, idx)

        # Everything from here on is just plotting nonsense.
        r_plot = 4*rvir_max(mergers[h_idx,:], scales)
        plot_frame(fig, ax, mergers[0,snap], mergers[h_idx,snap], scales[snap],
                   r_plot, x, idx, ranks.xc, mp_star)
        
        plt.savefig(path.join(out_dir, "frame_%03d.png" % frame_idx))
        frame_idx += 1

def rvir_max(halo, scale):
    i = np.argmax(halo["rvir"])
    return halo["rvir"][i]*scale[i]/H100*1e3

initial_r_half = None
dm_c_range = None
star_c_range = None

def plot_frame(fig, ax, host, sub, scale, r_max, x, idx, x_core, mp_star):
    global initial_r_half
    global dm_c_range
    global star_c_range
    
    dx = np.zeros(x.shape)
    for dim in range(3): dx[:,dim] = x[:,dim] - x_core[dim]

    mp_dm = np.ones(len(idx))*MP_DM
    mp_star = mp_star[idx]


    r_half = calc_r_half(dx, mp_star)
    if initial_r_half is None:
        initial_r_half = r_half
        
    ax[0].clear()
    ax[1].clear()
    
    extent = [-r_max, r_max, -r_max, r_max]
    r_max_star = initial_r_half*10
    
    if dm_c_range is None:
        dm_c_range = get_c_range(
            dx[:,0], dx[:,1], mp_dm, r_max, 100, 0.99)
    if star_c_range is None:
        star_c_range = get_c_range(
            dx[:,0], dx[:,1], mp_star, r_max_star, 50, 0.99)
    
    H, _, _, im = ax[0].hist2d(
        dx[:,0], dx[:,1], bins=100, weights=mp_dm,
        range=((-r_max, r_max), (-r_max, r_max)),
        norm=mpl_colors.LogNorm(dm_c_range[0], dm_c_range[1]),
        cmap="Greys"
    )

    H, _, _, im = ax[1].hist2d(
        dx[:,0], dx[:,1], bins=100, weights=mp_star,
        range=((-r_max_star, r_max_star), (-r_max_star, r_max_star)),
        norm=mpl_colors.LogNorm(star_c_range[0], star_c_range[1]),
        cmap="Greys",
    )

    rvir_sub = sub["rvir"]*1e3/H100*scale
    rvir_host = host["rvir"]*1e3/H100*scale
    x_sub = sub["x"]*1e3/H100*scale
    x_host = host["x"]*1e3/H100*scale

    if host["mvir"] > 0:
        plot_circle(ax[0], -x_core[0], -x_core[1], rvir_host, "tab:orange", 4)
        plot_circle(ax[1], -x_core[0], -x_core[1], rvir_host, "tab:orange", 5)

    if sub["mvir"] > 0:
        plot_circle(ax[0], (x_sub[0]-x_host[0])-x_core[0],
                    (x_sub[1]-x_host[1])-x_core[1], rvir_sub,
                    "tab:blue", 3)
        ax[1].plot([(x_sub[0]-x_host[0])-x_core[0]],
                   [(x_sub[1]-x_host[1])-x_core[1]], "x", color="tab:blue")
        plot_circle(ax[1], (x_sub[0]-x_host[0])-x_core[0],
                    (x_sub[1]-x_host[1])-x_core[1], rvir_sub,
                    "tab:blue", 4)

    plot_circle(ax[0], 0, 0, r_half, "tab:red", 1)
    plot_circle(ax[1], 0, 0, r_half, "tab:red", 3)
    
    ax[0].set_xlim(-r_max, r_max)
    ax[1].set_xlim(-12*initial_r_half, 12*initial_r_half)
    ax[0].set_ylim(-r_max, r_max)
    ax[1].set_ylim(-12*initial_r_half, 12*initial_r_half)
    ax[0].set_xlabel(r"$X\ ({\rm kpc})$")
    ax[1].set_xlabel(r"$X\ ({\rm kpc})$")
    ax[0].set_ylabel(r"$Y\ ({\rm kpc})$")

def plot_circle(ax, x0, y0, r, color, lw):
    th = np.linspace(0, 2*np.pi, 100)
    y = np.sin(th)*r + y0
    x = np.cos(th)*r + x0
    ax.plot(x, y, c=color, lw=lw)

def calc_r_half(dx, mp):
    r = np.sqrt(np.sum(dx**2, axis=1))
    r_min, r_max = max(np.min(r), 0.1), np.max(r)
    bins = np.linspace(np.log10(r_min), np.log10(r_max), 101)
    m, _ = np.histogram(np.log10(r), bins=bins, weights=mp)

    m_enc = np.cumsum(m)
    return 10**bins[np.searchsorted(m_enc, m_enc[-1]/2)+1]

def get_c_range(x, y, mp, r_max, bins, frac):
    range = ((-r_max, r_max), (-r_max, r_max))
    H, _, _ = np.histogram2d(x, y, weights=mp, range=range, bins=bins)
    H = np.sort(H.flatten())

    H_sum = np.cumsum(H)
    H_tot = H_sum[-1]
    i_min = np.searchsorted(H_sum, H_tot*(1 - frac)/2)
    i_max = np.searchsorted(H_sum, H_tot*(1 + frac)/2)

    return H[i_min], H[i_max]
    
    
if __name__ == "__main__": main()
