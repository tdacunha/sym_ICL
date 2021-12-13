import array
import struct
import numpy as np
import sys

MERGER_DTYPE = [("id", "i4"), ("mvir", "f4"), ("vmax", "f4"),
                ("x", "f4", (3,)), ("v", "f4", (3,))]

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

    return out

def main():
    h = read_mergers(sys.argv[1])
    print(h["mvir"][0,:])
    print(h["mvir"][2,:])

if __name__ == "__main__": main()
