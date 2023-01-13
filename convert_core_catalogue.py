import numpy as np
import struct
import symlib
import sys
import os
import os.path as path

def get_sim_dirs(config_name):
    with open(config_name, "r") as fp: text = fp.read()
    lines = [line for line in text.split("\n") if len(line) > 0]
    for i in range(len(lines)):
        lines[i] = [tok for tok in lines[i].split(" ") if len(tok) > 0]
    return [line[7] for line in lines]

def parse_sim_dir(sim_dir):
    suite_dir, halo = path.split(sim_dir)
    base_dir, suite = path.split(suite_dir)
    return base_dir, suite, halo

def get_matching_file_names(file_fmt):
    file_names = []
    while True:
        file_name = file_fmt % len(file_names)
        if not path.exists(file_name): return file_names
        file_names.append(file_name)

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
    flags = parse_flags(sys.argv[3:])

    target_idx = int(idx_str)
    sim_dirs = get_sim_dirs(config_name)

    if target_idx == -1:
        print("Must target a single halo.")
        exit(1)

    sim_dir = sim_dirs[target_idx]
    if sim_dir[-1] == "/": sim_dir = sim_dir[:-1]
    
    if "suffix" not in flags:
        file_fmt = path.join(sim_dir, "halos", "core.%d.txt")
    else:
        file_fmt = path.join(sim_dir, "halos", "core_%s.%%d.txt" %
                             flags["suffix"])
    file_names = get_matching_file_names(file_fmt)

    base_dir, suite_name, halo_name = parse_sim_dir(sim_dir)
    param = symlib.parameter_table[suite_name]
    h, hist = symlib.read_subhalos(sim_dir, include_false_selections=True)

    out = np.ones(h.shape, dtype=symlib.CORE_DTYPE)
    out["x"], out["v"] = -1, -1
    out["r_tidal"], out["r50_bound"], out["r50_bound_rs"] = -1, -1, -1
    out["m_tidal"], out["m_tidal_bound"], out["m_bound"] = -1, -1, -1

    for file_name in file_names:
        cols = np.loadtxt(file_name, dtype=int, usecols=(0,1)).T
        if cols.shape == (0,): continue

        snaps, subs = cols
        cols = np.loadtxt(file_name).T
        x, v = cols[2:5].T, cols[5:8].T
        r_tidal, r50_bound, r50_bound_rs = cols[8:11]
        m_tidal, m_tidal_bound, m_bound = cols[11:14]
        vmax, f_core, f_core_rs, d_core_mbp = cols[14:18]

        for i in range(len(subs)):
            sub, snap = subs[i], snaps[i]
            out[sub,snap]["x"], out[sub,snap]["v"] = x[i], v[i]
            out[sub,snap]["r_tidal"] = r_tidal[i]
            out[sub,snap]["m_tidal"] = m_tidal[i]
            out[sub,snap]["r50_bound"] = r50_bound[i]
            out[sub,snap]["r50_bound_rs"] = r50_bound_rs[i]
            out[sub,snap]["m_tidal"] = m_tidal[i]
            out[sub,snap]["m_tidal_bound"] = m_tidal_bound[i]
            out[sub,snap]["m_bound"] = m_bound[i]
            out[sub,snap]["vmax"] = vmax[i]
            out[sub,snap]["f_core"] = f_core[i]
            out[sub,snap]["f_core_rs"] = f_core_rs[i]
            out[sub,snap]["d_core_mbp"] = d_core_mbp[i]

    if "suffix" not in flags:
        file_root = "cores.dat"
    else:
        file_root = "cores_%s.dat"
        
    file_root = ("cores.dat" if "suffix" not in flags else
                 "cores_%s.dat" % flags["suffix"])
    out_file_name = path.join(sim_dir, "halos", file_root)
    out.reshape(out.shape[0]*out.shape[1])
    with open(out_file_name, "wb") as fp:
        fp.write(struct.pack("qq", out.shape[0], out.shape[1]))
        out["x"].tofile(fp)
        out["v"].tofile(fp)
        out["r_tidal"].tofile(fp)
        out["r50_bound"].tofile(fp)
        out["r50_bound_rs"].tofile(fp)
        out["m_tidal"].tofile(fp)
        out["m_tidal_bound"].tofile(fp)
        out["m_bound"].tofile(fp)
        out["vmax"].tofile(fp)
        out["f_core"].tofile(fp)
        out["f_core_rs"].tofile(fp)
        out["d_core_mbp"].tofile(fp)

        
    if "suffix" in flags:
        c = symlib.read_cores(sim_dir, suffix=flags["suffix"])
    else:
        c = symlib.read_cores(sim_dir)

if __name__ == "__main__": main()
