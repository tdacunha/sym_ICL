import star_tagging
import lib
import numpy as np
import matplotlib.pyplot as plt
import os.path as path
import sys

import palette
from palette import pc
palette.configure(False)

BASE_DIR_FMT ="/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo%03d/"
N_TRIALS = 100

gal_halo = star_tagging.GalaxyHaloModel(
    star_tagging.UniverseMachineMStar,
    star_tagging.Nadler2020RHalf,
    star_tagging.PlummerProfile,
)

def main():
    halo_id = int(sys.argv[1])
    base_dir = BASE_DIR_FMT % halo_id

    # Tree headers
    b_idx, m = lib.read_mergers(base_dir)
    b = lib.read_branches(base_dir)
    # Tree data
    x, mvir, snap = lib.read_tree(base_dir, ["X", "Mvir", "Snap"])
    scale = lib.scale_factors()

    # Merger information
    mpeak, merger_snap, merger_ratio = lib.merger_stats(b, m, x, mvir, snap)
    p_sub_idx = lib.pristine_merger_indices(b)

    # Plotting
    colors =[pc("r", 0.7), pc("o", 0.7), pc("g", 0.7),
             pc("b", 0.7), pc("p", 0.7), pc("r", 0.45),
             pc("o", 0.45), pc("g", 0.45), pc("b", 0.45), 
             pc("p", 0.45)]
    for i in range(len(m)):
        # index of each major merger in the list of main branches
        j = b_idx[i]
        # start and end of the branch
        start, end = b["start"][j], b["end"][j]
        # snapshtos of that branch
        snap_i = snap[start:end]
        # mvir of that branch
        mvir_i = mvir[start:end]
        # snapshots where the halo was pre-merger
        pre_merger = snap_i < merger_snap[np.searchsorted(p_sub_idx, j)]
        if i == 0: pre_merger = np.ones(len(snap_i), dtype=bool)

        c = colors[(i-1)%len(colors)]
        if i == 0: c = pc("k")
        # plot the whole MAH
        plt.plot(scale[snap_i], mvir_i, "--", c=c, lw=1.5)
        # plot the MAH
        plt.plot(scale[snap_i][pre_merger], mvir_i[pre_merger], lw=2.5, c=c)

    plt.yscale("log")
    plt.xscale("log")

    plt.ylabel(r"$M_{\rm vir}\,(h^{-1}\,M_\odot)$")
    plt.xlabel(r"$z+1$")

    plt.savefig("../plots/sub_mah_%03d.png" % halo_id)

if __name__ == "__main__": main()
