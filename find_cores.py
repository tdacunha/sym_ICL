import numpy as np
import symlib
import find_infall_cores
import scipy.spatial.distance as distance
import struct
import sys
import scipy.special as special
import os.path as path

N_CORE_MAX = 32

def spline(x):
    r1 = x <= 1
    r2 = (x > 1) & (x <= 2)
    
    out = np.zeros(len(x))
    out[r1] = 1 - (3/2)*x[r1]**2 + 0.75*x[r1]**3
    out[r2] = 0.25*(2 - x[r2])**3

    return out

def rho_sph(ri, n_dim):
    r = np.max(ri, axis=0)/2
    Vr = np.pi**(n_dim/2)*r**n_dim / special.gamma(n_dim/2 + 1)
    rho = np.zeros(len(r))
    
    # ri.shape[0] is small so this is okay
    for i in range(ri.shape[0]): rho += spline(ri[i,:]/r)
    return rho/Vr

def phase_space_density(x_core, v_core, x, v):
    # Use the Rockstar metric to get phase space densities.
    dx, dv = np.std(x_core), np.std(v_core)
    x_core, x = x_core/dx, x/dx
    v_core, v = v_core/dv, v/dv

    # Combine into one big vector
    z, z_core = np.zeros((len(x), 6)), np.zeros((len(x_core), 6))
    z[:,:3], z[:,3:] = x, v
    z_core[:,:3], z_core[:,3:] = x_core, v_core

    # len(x_core) is so small and all the points are so close together that
    # a KD tree is orders of magnitude slower than just doing it with brute
    # force.
    r = distance.cdist(z_core, z, "euclidean")
    return rho_sph(r, 6)
    
def write_cores(sim_dir, idxs):
    file_name = path.join(sim_dir, "halos", "cores.dat")

    with open(file_name, "wb") as fp:
        fp.write(struct.pack("qqq", idxs.shape[0],
                             idxs.shape[1], idxs.shape[2]))
        idxs = idxs.flatten()
        idxs = np.asarray(idxs, dtype=np.int32)
        idxs.tofile(fp)

def write_core_positions(sim_dir, x, v):
    file_name = path.join(sim_dir, "halos", "core_pos.dat")
    with open(file_name, "wb") as fp:
        fp.write(struct.pack("qq", x.shape[0], x.shape[1]))
        x, v = x.flatten(), v.flatten()
        x = np.asarray(x, dtype=np.float32)
        v = np.asarray(v, dtype=np.float32)
        x.tofile(fp)
        v.tofile(fp)
        
def core_position(idxs, x):
    return x[idxs[0]]

def core_velocity(idxs, v, scale):
    return np.sqrt(scale)*v[idxs[0]]

def track_cores_snap(info, sim_dir, snap, h, hist, infall_core_idxs,
                     core_idxs, core_x, core_v):
    print("snap", snap)

    base_dir, suite, halo_name = find_infall_cores.parse_sim_dir(sim_dir)    
    param = symlib.parameter_table[suite]
    scale = symlib.scale_factors(sim_dir)

    owner_all = symlib.read_particles(info, sim_dir, snap, "ownership")
    valid_all = symlib.read_particles(info, sim_dir, snap, "valid")
    x_all = symlib.read_particles(info, sim_dir, snap, "x")
    v_all = symlib.read_particles(info, sim_dir, snap, "v")

    for sub_i in range(1, len(h)):
        if snap < hist["merger_snap"][sub_i]:
            continue
        elif snap == hist["merger_snap"][sub_i]:
            n_core = np.min(N_CORE_MAX)
            x, v = x_all[sub_i], v_all[sub_i]

            idx = infall_core_idxs[sub_i,:n_core]
            core_idxs[sub_i,snap,:n_core] = idx
            core_x[sub_i,snap,:] = core_position(idx, x)
            core_v[sub_i,snap,:] = core_velocity(idx, v, scale[snap])
        else:
            owner, valid = owner_all[sub_i], valid_all[sub_i]
            x, v = x_all[sub_i], v_all[sub_i]
            
            n_core = np.sum(core_idxs[sub_i,snap-1,:] >= 0)
            x_core = x[core_idxs[sub_i,snap-1,:n_core]]
            v_core = v[core_idxs[sub_i,snap-1,:n_core]]

            ok = valid & (owner == 0)
            
            idx = np.arange(len(x), dtype=int)[ok]
            rho = phase_space_density(x_core, v_core, x[ok], v[ok])
            
            top_n_idx = np.argpartition(rho, -n_core)[-n_core:]
            top_n_rho = rho[top_n_idx]
            top_n_order = np.argsort(top_n_rho)[::-1]
            idx = idx[top_n_idx[top_n_order]]

            core_idxs[sub_i,snap,:n_core] = idx
            core_x[sub_i,snap,:] = core_position(idx, x)
            core_v[sub_i,snap,:] = core_velocity(idx, v, scale[snap])

            
def track_cores(sim_dir):
    base_dir, suite, halo_name = find_infall_cores.parse_sim_dir(sim_dir)
    
    param = symlib.parameter_table[suite]
    scale = symlib.scale_factors(sim_dir)
    n_snap = len(scale)
        
    h, hist = symlib.read_subhalos(param, sim_dir)
    
    info = symlib.ParticleInfo(sim_dir)

    core_idxs = np.ones((len(h), n_snap, N_CORE_MAX), dtype=int)*-1
    core_x = np.zeros((len(h), n_snap, 3))
    core_v = np.zeros((len(h), n_snap, 3))
    
    infall_core_idxs = symlib.read_particles(info, sim_dir, 0, "infall_core")
    for snap in range(np.min(hist["merger_snap"][1:]), n_snap):
        track_cores_snap(info, sim_dir, snap, h, hist, infall_core_idxs,
                         core_idxs, core_x, core_v)
        
    write_cores(sim_dir, core_idxs)
    write_core_positions(sim_dir, core_x, core_v)

def main():
    config_name, idx_str = sys.argv[1], sys.argv[2]
    target_idx = int(idx_str)

    sim_dirs = find_infall_cores.get_sim_dirs(config_name)
    
    for host_i in range(len(sim_dirs)):
        if host_i != target_idx and target_idx != -1: continue
        sim_dir = sim_dirs[host_i]
        if sim_dir[-1] == "/": sim_dir = sim_dir[:-1]

        track_cores(sim_dir)

if __name__ == "__main__": main()
