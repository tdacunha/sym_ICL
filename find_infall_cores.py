import sys
import os.path as path
import symlib
import numpy as np
import struct

N_CORE_MAX = 100

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

def write_cores(sim_dir, idxs):
    file_name = path.join(sim_dir, "halos", "cores.dat")

    with open(file_name, "wb") as fp:
        fp.write(struct.pack("qq", idxs.shape[0], idxs.shape[1]))
        idxs = idxs.flatten()
        idxs.tofile(fp)

def main():
    config_name, idx_str = sys.argv[1], sys.argv[2]
    target_idx = int(idx_str)

    sim_dirs = get_sim_dirs(config_name)
    
    for host_i in range(len(sim_dirs)):
        if host_i != target_idx and target_idx != -1: continue
        sim_dir = sim_dirs[host_i]
        if sim_dir[-1] == "/": sim_dir = sim_dir[:-1]

        base_dir, suite, halo_name = parse_sim_dir(sim_dir)

        param = symlib.parameter_table[suite]
        scale = symlib.scale_factors(sim_dir)
        n_snap = len(scale)

        h, hist = symlib.read_subhalos(param, sim_dir)

        info = symlib.ParticleInfo(sim_dir)

        core_idxs = np.ones((len(h), N_CORE_MAX), dtype=np.int32) * -1

        owner_all = symlib.read_particles(info, sim_dir, 0, "ownership")

        for snap in range(n_snap):
            sub_idxs = np.where(snap == hist["merger_snap"])[0]
            if len(sub_idxs) == 0: continue
            if len(sub_idxs) == 1 and sub_idxs[0] == 0: continue

            x_all = symlib.read_particles(info, sim_dir, snap, "x")
            v_all = symlib.read_particles(info, sim_dir, snap, "v")

            valid_all = symlib.read_particles(info, sim_dir, snap, "valid")
            
            for sub_i in sub_idxs:
                if sub_i == 0: continue
                print("snap: %d, sub: %d" % (snap, sub_i))
                valid, owner = valid_all[sub_i], owner_all[sub_i]
                x = x_all[sub_i][valid]
                v = v_all[sub_i][valid]
                p_idx = np.arange(len(x_all[sub_i]), dtype=np.int32)[valid]

                v *= np.sqrt(scale[snap])

                dx, dv = np.zeros(x.shape), np.zeros(v.shape)
                for dim in range(3):
                    dx[:,dim] = x[:,dim] - h[sub_i,snap]["x"][dim]
                    dv[:,dim] = v[:,dim] - h[sub_i,snap]["v"][dim]

                dx *= scale[snap]/param["h100"]

                rmax, vmax, PE, order = symlib.profile_info(param, x)
                KE = np.sum(dv**2, axis=1)/2 / vmax**2
                E = PE + KE

                E, p_idx = E[owner[valid] == 0], p_idx[owner[valid] == 0]
                order = np.argsort(E)
                E, p_idx = E[order], p_idx[order]

                if len(p_idx) > N_CORE_MAX:
                    core_idxs[sub_i,:] = p_idx[:N_CORE_MAX]
                else:
                    core_idxs[sub_i,:len(p_idx)] = p_idx
                    
        write_cores(sim_dir, core_idxs)
        print(core_idxs[1])
        
        core_idxs_rd = symlib.read_particles(info, sim_dir, 0, "core")
        print(core_idxs_rd[1])

if __name__ == "__main__": main()
