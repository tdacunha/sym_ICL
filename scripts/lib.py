import array
import struct
import numpy as np
import os
import os.path as path

MERGER_DTYPE = [("id", "i4"), ("mvir", "f4"), ("vmax", "f4"),
                ("x", "f4", (3,)), ("v", "f4", (3,))]

BRANCHES_DTYPE = [("start", "i4"), ("end", "i4"), ("is_real", "?"),
                  ("is_disappear", "?"), ("is_main_sub", "?"),
                  ("preprocess", "i4")]

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
    return 10**np.linspace(np.log10(a_start), np.log10(a_end), n_snap)

def mvir_to_rvir(mvir, a, omega_M):
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

def read_mergers(fname):
    f = open(fname, "rb")

    n_snap = struct.unpack("i", f.read(4))[0]
    n_merger = struct.unpack("i", f.read(4))[0]
    idx = array.array("i")
    
    idx.fromfile(f, n_merger)
    idx = np.array(idx, dtype=int)
    
    out = np.zeros((n_merger+1, n_snap), dtype=MERGER_DTYPE)

    for i in range(n_merger):
        mvir = array.array("f")
        mvir.fromfile(f, n_snap)
        out["mvir"][i,:] = np.array(mvir, dtype=float)
        
    for i in range(n_merger):
        vmax = array.array("f")
        vmax.fromfile(f, n_snap)
        out["vmax"][i,:] = np.array(vmax, dtype=float)

    for i in range(n_merger):
        id = array.array("i")
        id.fromfile(f, n_snap)
        out["id"][i,:] = np.array(id, dtype=int)

    for i in range(n_merger):
        x = array.array("f")
        x.fromfile(f, n_snap*3)
        out["x"][i,:,:] = np.array(x, dtype=float).reshape((n_snap, 3))

    for i in range(n_merger):
        v = array.array("f")
        v.fromfile(f, n_snap*3)
        out["v"][i,:,:] = np.array(x, dtype=float).reshape((n_snap, 3))
        
    f.close()

    return idx, out

def read_branches(fname):
    f = open(fname, "rb")
    
    n = struct.unpack("i", f.read(4))[0]
    central_idx = struct.unpack("i", f.read(4))[0]
    out = np.zeros(n, dtype=BRANCHES_DTYPE)
    
    edges = array.array("i")
    edges.fromfile(f, n+1)
    
    edges = np.array(edges, dtype=np.int32)
    out["start"] = edges[:-1]
    out["end"] = edges[1:]
    
    is_real = array.array("b")
    is_real.fromfile(f, n)
    out["is_real"] = np.asarray(is_real, dtype=np.bool)
    
    is_disappear = array.array("b")
    is_disappear.fromfile(f, n)
    out["is_disappear"] = np.asarray(is_disappear, dtype=np.bool)
    
    is_main_sub = array.array("b")
    is_main_sub.fromfile(f, n)
    out["is_main_sub"] = np.asarray(is_main_sub, dtype=np.bool)
    
    preprocess = array.array("i")
    preprocess.fromfile(f, n)
    out["preprocess"] = np.asarray(preprocess, dtype=np.int32)
    
    return central_idx, out

def read_tree(dir_name, var_names):
    paths = [path.join(dir_name, fname) for fname in os.listdir(dir_name)]
    tree_files = [p for p in paths if path.isfile(p) and
                  len(p) > 6 and p[-6:] == "df.bin"]
    hd = read_tree_header(tree_files[0])

    out = []
    for i in range(len(var_names)):
        offset = tree_var_offset(hd, var_names[i])

        var = []
        for j in range(len(tree_files)):
            f = open(tree_files[j], "rb")
            col = tree_var_col(hd, var_names[i])
            f.seek(offset)
        
            if is_int(col, hd):
                x = array.array("i")
                x.fromfile(f, hd.n)
                var.append(np.array(x, dtype=np.int32))
            elif is_float(col, hd):
                x = array.array("f")
                x.fromfile(f, hd.n)
                var.append(np.array(x, dtype=np.float))
            else:
                x = array.array("f")
                x.fromfile(f, hd.n*3)
                var.append(np.array(x, dtype=np.float).reshape((hd.n, 3)))

        out.append(np.hstack(var))
    return out
             
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
