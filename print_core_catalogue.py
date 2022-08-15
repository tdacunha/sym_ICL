import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import subhalo_tracking as sh
import os.path as path
import symlib
import matplotlib.colors as mpl_colors
from colossus.cosmology import cosmology
from colossus.halo import mass_so
import subfind
import os
import sys

def get_sim_dirs(config_name):
    with open(config_name, "r") as fp: text = fp.read()
    lines = [line for line in text.split("\n") if len(line) > 0]
    for i in range(len(lines)):
        lines[i] = [tok for tok in lines[i].split(" ") if len(tok) > 0]
    return [line[7] for line in lines]

def parse_sim_dir(sim_dir):
    suite_dir, halo = path.split(sim_dir)
    base_dir, suite = path.split(suite_dir)
    print(sim_dir, base_dir, suite, halo)
    return base_dir, suite, halo

def get_matching_file_names(file_fmt):
    file_names = []
    while True:
        file_name = file_fmt % len(file_names)
        if not path.exists(file_name): return file_names
        file_names.append(file_name)

def calculate_min_snap(prev_file_names):
    snaps = []
    for name in prev_file_names:
        snaps.append(np.loadtxt(name, usecols=(0,), dtype=int))

    if len(snaps) == 0: return 0
    snaps = np.hstack(snaps)
    if len(snaps) == 0: return 0
    # We don't know how far the program was in calculating this snapshot, so
    # redo the last one we found.
    return int(np.max(snaps))

def main():
    palette.configure(False)

    config_name, idx_str = sys.argv[1], sys.argv[2]
    target_idx = int(idx_str)
    sim_dirs = get_sim_dirs(config_name)

    if len(sys.argv) == 4:
        if sys.argv[3] == "--reset":
            reset_files = True
        else:
            raise ValueError("Unrecognized flag, '%s'" % sys.argv[3])
    else:
        reset_files = False

    if target_idx == -1:
        raise ValueError("Must target a single halo.")

    sim_dir = sim_dirs[target_idx]
    if sim_dir[-1] == "/": sim_dir = sim_dir[:-1]

    out_file_fmt = path.join(sim_dir, "halos", "core.%d.txt")
    prev_file_names = get_matching_file_names(out_file_fmt)

    if reset_files:
        for name in prev_file_names: os.remove(name)
        out_file = out_file_fmt % 0
    else:
        out_file = out_file_fmt % len(prev_file_names)
    
    base_dir, suite_name, halo_name = parse_sim_dir(sim_dir)
    
    param = symlib.parameter_table[suite_name]
    h, hist = symlib.read_subhalos(sim_dir, include_false_selections=True,
                                   comoving=True)
    h_cmov = np.copy(h)
    info = symlib.ParticleInfo(sim_dir)

    scale = symlib.scale_factors(sim_dir)
    cosmo = cosmology.setCosmology("", symlib.colossus_parameters(param))
    h = symlib.set_units_halos(h_cmov, scale, param)

    targets = np.arange(1, len(h), dtype=int)
    
    n_core = 32
    tracks = [None]*len(h)

    if reset_files:
        min_snap = 0
    else:
        min_snap = calculate_min_snap(prev_file_names)

    max_snap = len(scale) - 1
    starting_snap = np.maximum(hist["merger_snap"], min_snap)
    infall_cores = [None]*len(h)

    with open(out_file, "a") as fp:
        print("""# 0 - snap
# 1 - subhalo index
# 2-4 - x (pkpc)
# 5-7 - v (km/s)
# 8 - R_tidal (pkpc)
# 9 - R_50,bound (pkpc)
# 10 - R_95,bound (pkpc)
# 11 - M_tidal (Msun)
# 12 - M_tidal,bound (Msun)
# 13 - M_bound (msun)""", file=fp)

    for snap in range(np.min(starting_snap[targets]), max_snap + 1):

        print(snap)

        sd = sh.SnapshotData(info, sim_dir, snap, scale[snap], h_cmov, param,
                             include_false_selections=True)
        prof = sh.MassProfile(sd.param, snap, h, sd.x, sd.owner, sd.valid)
        
        for j in range(len(targets)):
            i_sub = targets[j]
            if snap < starting_snap[i_sub]: continue
            if hist["false_selection"][i_sub]: continue
            if -1 in sd.infall_cores[i_sub]: continue
            print("   ", i_sub)

            if tracks[i_sub] is None:
                print(tracks)
                tracks[i_sub] = sh.SubhaloTrack(
                    i_sub, sd, sd.infall_cores[i_sub], param)
            else:
                tracks[i_sub].next_snap(sd, sd.infall_cores[i_sub])
                
            xp, vp = sd.x[i_sub], sd.v[i_sub]
            mp = np.ones(len(xp))*sd.mp
        
            xc = tracks[i_sub].x[snap]
            vc = tracks[i_sub].v[snap]

            dxp, dvp = sh.delta(xp, xc), sh.delta(vp, vc)
            ok, valid = sd.ok[i_sub], sd.valid[i_sub]

            bound_only=True
            r_tidal, m_tidal = prof.tidal_radius(
                xc, xp[valid], vc, vp[valid], bound_only=False)            

            is_bound,_ = sh.is_bound_iter(10, sd.param, dxp, dvp, ok=valid)
            m_bound = np.sum(mp[is_bound])
            
            r = np.sqrt(np.sum(dxp**2, axis=1))
            if np.sum(is_bound) > 2:
                r_50_bound = np.quantile(r[is_bound], 0.5)
                r_95_bound = np.quantile(r[is_bound], 0.95)
            else:
                r_50_bound, r_95_bound = 0, 0

            m_tidal_bound = np.sum(mp[(r < r_tidal) & is_bound])
            
            with open(out_file, "a") as fp:
                print(("%d %d "+"%.4f "*3+"%.4f "*3+"%.4f "*3+"%.4g "*3) %
                      (snap, i_sub, xc[0], xc[1], xc[2], vc[0], vc[1], vc[2],
                       r_tidal, r_50_bound, r_95_bound,
                       m_tidal, m_tidal_bound, m_bound), file=fp)

        
if __name__ == "__main__": main()
