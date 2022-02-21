import numpy as np
import lib

DIR_FORMAT = "/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo%03d"
HALO_NUMS = [4, 113, 169, 170, 222, 229, 282, 327, 349, 407, 453, 476, 523, 625,
             659, 666, 719, 747, 756, 788, 858, 953, 975, 983, 908]
DIR_NAMES = [DIR_FORMAT % n for n in HALO_NUMS]
MP = 2.8e5
MVIR_CONV = MP * 300
N_SNAP = 236

def main():
    f = open("/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/sub_association.txt", "w")
    print("# 0 - Halo name", file=f)
    print("# 1 - Flag (0: direct infall, 1: preprocessed by small halo,",file=f)
    print("#           2: preprocessed by GSE, 3: preprocessed by LMC)", file=f)
    print("# 2 - Snapshot", file=f)
    print("# 3 - ID", file=f)

    for dir_name in DIR_NAMES:
        halo_name = dir_name.split("/")[-1]
        print(halo_name)

        snap, id, mvir = lib.read_tree(dir_name, ["Snap", "ID", "Mvir"])
        m_idx, m = lib.read_mergers(dir_name)
        b = lib.read_branches(dir_name)

        survives = snap[b["start"]] == N_SNAP - 1
        is_mw_sub = np.where(b["is_main_sub"] & b["is_real"] &
                             survives & (~b["is_disappear"]))[0]

        mw_idx = m_idx[0]
        lmc_idx, _, gse_idx = lib.read_merger_idxs(dir_name)

        id0 = id[b["start"][is_mw_sub]]
        snap0 = snap[b["start"][is_mw_sub]]
        flag = np.zeros(len(id0), dtype=int)
        flag[b["preprocess"][is_mw_sub] == -1] = 1
        flag[b["preprocess"][is_mw_sub] == gse_idx] = 2
        flag[b["preprocess"][is_mw_sub] == lmc_idx] = 3

        for i in range(len(is_mw_sub)):
            j = is_mw_sub[i]
            mpeak = np.max(mvir[b["start"][j]: b["end"][j]])
            if mpeak > MVIR_CONV:
                print(halo_name, flag[i], snap0[i], id0[i], file=f)

    f.close()

if __name__ == "__main__": main()
