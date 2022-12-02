import numpy as np
import matplotlib.pyplot as plt
import subhalo_tracking as sh
import os.path as path
import symlib
import matplotlib.colors as mpl_colors
from colossus.cosmology import cosmology
from colossus.halo import mass_so
import subfind
import os
import sys
import gravitree

SUBFIND_K = 16
VOTING_NP = 32

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

def parse_flags(flags):
    out = { }
    for flag in flags:
        if len(flag) < 3 or flag[:2] != "--": continue
        tok = flag[2:].split("=")
        if len(tok) == 1:
            out[tok[0]] = ""
        else:
            out[tok[0]] = tok[1]
    return out

def main():
    config_name, idx_str = sys.argv[1], sys.argv[2]
    target_idx = int(idx_str)

    flags = parse_flags(sys.argv[3:])
    n_halo = len(np.loadtxt(config_name, dtype=str))

    for i in range(n_halo):
        if target_idx == -1 or i == target_idx:
            print_halo(config_name, i, flags)

def print_halo(config_name, target_idx, flags):
    if "reset" in flags:
        reset_files = True
    else:
        reset_files = False

    if target_idx == -1:
        raise ValueError("Must target a single halo.")

    sim_dirs = get_sim_dirs(config_name)
    sim_dir = sim_dirs[target_idx]
    if sim_dir[-1] == "/": sim_dir = sim_dir[:-1]

    if "suffix" in flags:
        out_file_fmt = path.join(sim_dir, "halos", "core_%s.%%d.txt" %
                                 flags["suffix"])
    else:
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
# 10 - R_50,bound,Rockstar (pkpc)
# 11 - M_tidal (Msun)
# 12 - M_tidal,bound (Msun)
# 13 - M_bound (msun)
# 14 - Vmax (km/s)
# 15 - f_core,32
# 16 - f_core,32,rs
# 17 - d_core,mbp (pkpc)""", file=fp)

    start_snap = np.min(starting_snap[targets])
    end_snap = max_snap
    if "snap_range" in flags:
        low, high = flags["snap_range"].split(":")
        start_snap = max(start_snap, int(low))
        end_snap = min(end_snap, int(high))

    for snap in range(start_snap, end_snap + 1):
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
                tracks[i_sub] = sh.SubhaloTrack(
                    i_sub, sd, sd.infall_cores[i_sub][:VOTING_NP],
                    param, SUBFIND_K)
            else:
                tracks[i_sub].next_snap(sd, sd.infall_cores[i_sub][:VOTING_NP])
                
            xp, vp = sd.x[i_sub], sd.v[i_sub]
            mp = np.ones(len(xp))*sd.mp
        
            xc = tracks[i_sub].x[snap]
            vc = tracks[i_sub].v[snap]

            dxp, dvp = sh.delta(xp, xc), sh.delta(vp, vc)
            ok, valid = sd.ok[i_sub], sd.valid[i_sub]

            r_tidal, m_tidal = prof.tidal_radius(
                xc, xp[valid], vc, vp[valid], bound_only=True)

            is_bound = gravitree.binding_energy(
                dxp[valid], dvp[valid], sd.mp, sd.eps, n_iter=10) < 0
            
            m_bound = np.sum(is_bound)*sd.mp
            
            r = np.sqrt(np.sum(dxp[valid]**2, axis=1))
            if np.sum(is_bound) > 2:
                r_50_bound = np.quantile(r[is_bound], 0.5)
                _, vmax = symlib.rmax_vmax(param, dxp[valid], is_bound)
            else:
                r_50_bound = 0.0
                vmax = 0.0

            m_tidal_bound = np.sum((r < r_tidal) & is_bound)*sd.mp

            dxh = dxp + xc - h[i_sub,snap]["x"]
            if h["ok"][i_sub,snap]:
                dvh = dvp + vc - h[i_sub,snap]["v"]
                is_bound_h = gravitree.binding_energy(
                    dxh[valid], dvh[valid], sd.mp, sd.eps, n_iter=10) < 0
                rh = np.sqrt(np.sum(dxh[valid]**2, axis=1))
                if np.sum(is_bound_h) > 2:
                    r_50_bound_h = np.quantile(rh[is_bound_h], 0.5)
                else:
                    r_50_bound_h = 0.0
            else:
                r_50_bound_h = 0.0
            
            cores = sd.infall_cores[i_sub]
            r_core = np.sqrt(np.sum(dxp[cores]**2, axis=1))
            r_h = np.sqrt(np.sum(dxh[cores]**2, axis=1))

            f_core_32 = np.sum(r_core < r_50_bound)/32
            f_core_32_rs = np.sum(r_h < r_50_bound_h)/32

            with open(out_file, "a") as fp:
                print(("%d %d "+"%.6f "*3+"%.6f "*3+"%.6f "*3+"%.6g "*3+"%.6g "+"%.6f "+"%.6f "+"%6.f") %
                      (snap, i_sub, xc[0], xc[1], xc[2], vc[0], vc[1], vc[2],
                       r_tidal, r_50_bound, r_50_bound_h,
                       m_tidal, m_tidal_bound, m_bound,
                       vmax, f_core_32, f_core_32_rs, r_core[0]), file=fp)
        
if __name__ == "__main__": main()
