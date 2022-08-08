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

def main():
    config_name, idx_str = sys.argv[1], sys.argv[2]
    target_idx = int(idx_str)
    sim_dirs = get_sim_dirs(config_name)

    if target_idx == -1:
        print("Must target a single halo.")
        exit(1)

    sim_dir = sim_dirs[target_idx]
    if sim_dir[-1] == "/": sim_dir = sim_dir[:-1]
    
    file_fmt = path.join(sim_dir, "halos", "core.%d.txt")
    file_names = get_matching_file_names(file_fmt)

    base_dir, suite_name, halo_name = parse_sim_dir(sim_dir)
    param = symlib.parameter_table[suite_name]
    h, hist = symlib.read_subhalos(param, sim_dir)

    out = np.ones(h.shape, dtype=symlib.CORE_DTYPE)
    out["x"], out["v"] = -1, -1
    out["r_tidal"], out["r50_bound"], out["r95_bound"] = -1, -1, -1
    out["m_tidal"], out["m_tidal_bound"], out["m_bound"] = -1, -1, -1

    for file_name in file_names:
        snaps, subs = np.loadtxt(file_name, dtype=int, usecols=(0,1)).T
        cols = np.loadtxt(file_name).T
        x, v = cols[2:5].T, cols[5:8].T
        r_tidal, r50_bound, r95_bound = cols[8:11]
        m_tidal, m_tidal_bound, m_bound = cols[11:14]
        
        for i in range(len(subs)):
            sub, snap = subs[i], snaps[i]
            out[sub,snap]["x"], out[sub,snap]["v"] = x[i], v[i]
            out[sub,snap]["r_tidal"] = r_tidal[i]
            out[sub,snap]["m_tidal"] = m_tidal[i]
            out[sub,snap]["r50_bound"] = r50_bound[i]
            out[sub,snap]["r95_bound"] = r95_bound[i]
            out[sub,snap]["m_tidal"] = m_tidal[i]
            out[sub,snap]["m_tidal_bound"] = m_tidal_bound[i]
            out[sub,snap]["m_bound"] = m_bound[i]

    out_file_name = path.join(sim_dir, "halos", "cores.dat")
    out.reshape(out.shape[0]*out.shape[1])
    with open(out_file_name, "wb") as fp:
        fp.write(struct.pack("qq", out.shape[0], out.shape[1]))
        out["x"].tofile(fp)
        out["v"].tofile(fp)
        out["r_tidal"].tofile(fp)
        out["r50_bound"].tofile(fp)
        out["r95_bound"].tofile(fp)
        out["m_tidal"].tofile(fp)
        out["m_tidal_bound"].tofile(fp)
        out["m_bound"].tofile(fp)

    c = symlib.read_cores(sim_dir)

if __name__ == "__main__": main()
