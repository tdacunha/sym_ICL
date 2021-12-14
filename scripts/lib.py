import array
import struct
import numpy as np
import os
import os.path as path

""" MERGER_DTYPE is the numpy datatype used by the main return value of
read_mergers(). Positions and distances are in comvoing Mpc/h, velocities are
physical peculiar velocities, and masses are in Msun/h.
"""
MERGER_DTYPE = [("id", "i4"), ("mvir", "f4"), ("vmax", "f4"),
                ("x", "f4", (3,)), ("v", "f4", (3,))]

""" BRANCHES_DTYPE is the numpy datatype used by the main return value of 
read_branches(). 
 - The main branch of a given halo within the depth-first merger tree can be 
   found with x[start:end].
 - is_real: false if a halo is definitely a numerical artefact (e.g. if it's
   already a subhalo in its first snapshot).
 - is_disappear: true if the subhalo disappears without merging. This is also
   bad and probably means the thing is a numerical artefact.
 - is_main_sub: true is the halo was ever a subhalo of the zoom-in box's main 
   halo (includes splashback subhaloes).
 - preprocess: the index of the branch of the largest halo that hosted this 
   halo before it became a subhalo of the main halo. If the halo was never a
   subhalo of the main halo, this is just the index of largest halo to have
   every hosted this halo. If no other halo every hosted this halo before
   it entered Rvir of the main halo, this is -1. Includes splashback subhaloes
   and doesn't include collisions between subhaloes once they've already been
   accreted.
"""
BRANCHES_DTYPE = [("start", "i4"), ("end", "i4"), ("is_real", "?"),
                  ("is_disappear", "?"), ("is_main_sub", "?"),
                  ("preprocess", "i4")]

""" TREE_COL_NAMES is the mapping of variable names to columns in the
consistent-trees file. These are the variable names you need to pass to
read_tree().
"""
TREE_COL_NAMES = {
    "DFID": 28,
    "ID": 1,
    "DescID": 3,
    "UPID": 6,
    "Phantom": 8,
    "Snap": 31,
    "NextProg": 32,
    "Mvir": 10,
    "Rs": 12,
    "Vmax": 16,
    "M200b": 39,
    "M200c": 40,
    "M500c": 41,
    "Xoff": 43,
    "SpinBullock": 45,
    "BToA": 46,
    "CToA": 47,
    "VirialRatio": 56,
    "X": 17,
    "V": 20,
    "J": 23,
    "A": 48,
}

def scale_factors(n_snap=236, a_start=1/20.0, a_end=1.0):
    """ scale_factors returns the scale factors used by the simulation. The
    defaults are set to correspond to the MWest suite's parameters.
    """
    return 10**np.linspace(np.log10(a_start), np.log10(a_end), n_snap)

def mvir_to_rvir(mvir, a, omega_M):
    """ mvir_to_rvir converts a Bryan & Norman virial mass in Msun/h to a virial
    radius in comoving Mpc/h at a given scale factor a, and omega_M.
    """
    omega_L = 1 - omega_M
    Ez = np.sqrt(omega_M/a**3 + omega_L)
    rho_crit = 2.77519737e11*Ez**2
    omega_Mz = (omega_M/a**3)/Ez**2

    rho_m = omega_Mz * rho_crit

    x = omega_Mz - 1
    delta_vir = 18*np.pi**2 + 82*x - 39.0*x**2
    rho_vir = rho_crit*delta_vir

    r_phys = (mvir/(rho_vir * (4*np.pi / 3)))**(1.0/3)
    r_cmov = r_phys/a

    return r_cmov

def flatten(arrays):
    """ flatten takes a list of numpy arrays and flattens it into a single
    array. arrays[i].shape[0] can be any value, but all other components of the
    shape vectors must be the same. THis is needed because hstack doesn't work
    on tensor-arrays.
    """

    N = sum(arr.shape[0] for arr in arrays)

    shape, dtype = arrays[0].shape, arrays[0].dtype
    if len(shape) == 1:
        out = np.zeros(N, dtype=dtype)
    else:
        out = np.zeros((N,) + shape[1:], dtype=dtype)

    start, end = 0, 0
    for i in range(len(arrays)):
        end += arrays[i].shape[0]
        out[start: end] = arrays[i]
        start = end

    return out

def read_mergers(dir_name):
    """ read_mergers reads major merger data from the halo directory dir_name.
    It returns two arrays. The first, m_idx, is the indices of the major
    mergers within the branches arrays. Index 0 is the main halo, index 1 is the
    biggest merger (by Mpeak), index 2 is the second biggest, etc. As a
    convenience, data from the main branches of those haloes are returned as
    the second argument, m. m[halo_idx, snap] gives the  properties of the
    halo at the index halo_idx in m_idx and the snapshot snap. The fields are
    given by MERGER_DTYPE (see the header of this file). If a halo doesn't
    exist at a given snapshot, values are set -1 ok will be set to false.
    """
    fname = path.join(dir_name, "mergers.dat")
    f = open(fname, "rb")

    n_snap = struct.unpack("i", f.read(4))[0]
    n_merger = struct.unpack("i", f.read(4))[0]

    idx = np.fromfile(f, np.int32, n_merger)    
    out = np.zeros((n_merger+1, n_snap), dtype=MERGER_DTYPE)

    for i in range(n_merger):
        out["mvir"][i,:] = np.fromfile(f, np.float32, n_snap)
    for i in range(n_merger):
        out["vmax"][i,:] = np.fromfile(f, np.float32, n_snap)
    for i in range(n_merger):
        out["id"][i,:] = np.fromfile(f, np.int32, n_snap)
    for i in range(n_merger):
        out["x"][i,:,:] = np.fromfile(f, (np.float32, (3,)), n_snap)
    for i in range(n_merger):
        out["v"][i,:,:] = np.fromfile(f, (np.float32, (3,)), n_snap)
        
    f.close()

    return idx, out

def read_branches(dir_name):
    """ read_branches reads main branch data from the halo directory dir_name.
    It returns an array with length n_branches where each element has type
    BRANCHES_DTYPE.
    """
    fname = path.join(dir_name, "branches.dat")
    f = open(fname, "rb")
    
    n = struct.unpack("i", f.read(4))[0]
    central_idx = struct.unpack("i", f.read(4))[0]
    out = np.zeros(n, dtype=BRANCHES_DTYPE)

    edges = np.fromfile(f, np.int32, n+1)
    out["start"] = edges[:-1]
    out["end"] = edges[1:]
    out["is_real"] = np.fromfile(f, np.bool, n)
    out["is_disappear"] = np.fromfile(f, np.bool, n)
    out["is_main_sub"] = np.fromfile(f, np.bool, n)
    out["preprocess"] = np.fromfile(f, np.int32, n)
    
    return out

def read_tree(dir_name, var_names):
    """ read_tree reads variables from the halo directory halo_dir in
    depth-first order. var_names is a list of variables to be read (see
    TREE_COL_NAMES). A list of arrays is returned. Use the branches and merger
    files to identify main branches and important haloes, respectively.
    """
    paths = [path.join(dir_name, fname) for fname in os.listdir(dir_name)]
    tree_files = [p for p in paths if path.isfile(p) and
                  len(p) > 6 and p[-6:] == "df.bin"]

    out = []
    for i in range(len(var_names)):
        var = []
        for j in range(len(tree_files)):
            hd = read_tree_header(tree_files[j])
            offset = tree_var_offset(hd, var_names[i])
            f = open(tree_files[j], "rb")
            col = tree_var_col(hd, var_names[i])
            f.seek(offset)
        
            if is_int(col, hd):
                var.append(np.fromfile(f, np.int32, hd.n))
            elif is_float(col, hd):
                var.append(np.fromfile(f, np.float32, hd.n))
            else:
                var.append(np.fromfile(f, (np.float32, (3,)), hd.n))

        out.append(flatten(var))
    return out

# Everything else in this file is an internal helper function

def is_int(col, hd): return col < hd.n_int
def is_float(col, hd): return col >= hd.n_int and col < hd.n_float
     
def read_tree_header(fname):
    f = open(fname, "rb")
    class Header(object): pass
    hd = Header()
    hd.n, hd.n_int, hd.n_float, hd.n_vec = struct.unpack("iiii", f.read(16))
    cols = array.array("i")
    cols.fromfile(f, hd.n)
    hd.cols = np.asarray(cols, dtype=np.int32)
    return hd

def tree_var_col(hd, var_name):
    if var_name not in TREE_COL_NAMES:
        raise ValueError("Variable name '%s' unrecognized" % var_name)
    
    target_col = TREE_COL_NAMES[var_name]
    for i in range(len(hd.cols)):
        if hd.cols[i] == target_col: break

    return i

def tree_var_offset(hd, var_name):
    i = tree_var_col(hd, var_name)
    hd_size = 4*4 + 4*(hd.n_int + hd.n_float + hd.n_vec)
    return hd.n*4*(i + 2*max(0, i - hd.n_int - hd.n_float)) + hd_size    

def main():
    import sys
    dir_name = sys.argv[1]
    m_idx, m = read_mergers(path.join(dir_name, "mergers.dat"))
    ci, b = read_branches(path.join(dir_name, "branches.dat"))
    
    mvir, x, snap = read_tree(dir_name, ["Mvir", "X", "Snap"])

    for i in range(len(m_idx)):
        print(i)
        s, e = b["start"][m_idx[i]], b["end"][m_idx[i]]
        print(mvir[s:s+5])
        print(x[s:s+5])
        print(snap[s:s+5])
        print()
    

if __name__ == "__main__": main()
