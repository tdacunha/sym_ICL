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
import struct

N_CORE = 32

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

def write_infall_cores(sim_dir, idxs):
    file_name = path.join(sim_dir, "halos", "infall_cores.dat")

    with open(file_name, "wb") as fp:
        fp.write(struct.pack("qq", idxs.shape[0], idxs.shape[1]))
        idxs = idxs.flatten()
        idxs.tofile(fp)

def main():
    palette.configure(False)

    config_name, idx_str = sys.argv[1], sys.argv[2]
    target_idx = int(idx_str)
    sim_dirs = get_sim_dirs(config_name)

    for host_i in range(len(sim_dirs)):
        if host_i != target_idx and target_idx != -1: continue
        sim_dir = sim_dirs[host_i]
        if sim_dir[-1] == "/": sim_dir = sim_dir[:-1]
        
        base_dir, suite_name, halo_name = parse_sim_dir(sim_dir)
        print(suite_name, halo_name)

        param = symlib.simulation_parameters(sim_dir)
        h, hist = symlib.read_subhalos(sim_dir, comoving=True,
                                       include_false_selections=True)
        h_cmov = np.copy(h)
        info = symlib.ParticleInfo(sim_dir)

        scale = symlib.scale_factors(sim_dir)
        cosmo = cosmology.setCosmology("", symlib.colossus_parameters(param))
        h = symlib.set_units_halos(h_cmov, scale, param)
    
        infall_cores = np.zeros((len(h), N_CORE), dtype=np.int32)

        for snap in range(len(scale)):
            if snap not in hist["merger_snap"][1:]: continue
            print("  snap %3d" % snap)

            sd = sh.SnapshotData(info, sim_dir, snap, scale[snap], h_cmov,
                                 param, include_false_selections=True)
            prof = sh.MassProfile(sd.param, snap, h, sd.x, sd.owner, sd.valid)
        
            sub_idxs = np.where(hist["merger_snap"] == snap)[0]

            for i_sub in sub_idxs:
                if i_sub == 0 or hist["false_selection"][i_sub]: continue
                print("   ", i_sub)

                merger_snap = hist["merger_snap"][i_sub]
                m_merger = h["mvir"][i_sub,hist["merger_snap"][i_sub]]
                mp = param["mp"]
                peak_snap = np.argmax(h["mvir"][i_sub,:])
                m_peak = hist["mpeak"][i_sub]
                
                infall_cores[i_sub] = sh.n_most_bound(
                    h["x"][i_sub,snap], h["v"][i_sub,snap],
                    sd.x[i_sub], sd.v[i_sub], sd.ok[i_sub],
                    N_CORE, sd.param
                )

        write_infall_cores(sim_dir, infall_cores)

        
if __name__ == "__main__": main()
