import array
import struct
import numpy as np
import sys

MERGER_DTYPE = [("id", "i4"), ("mvir", "f4"), ("vmax", "f4"),
                ("x", "f4", (3,)), ("v", "f4", (3,))]

def read_mergers(fname):
    f = open(fname, "rb")
    
    n_merger = struct.unpack("q", f.read(8))[0]
    n_snap = struct.unpack("q", f.read(8))[0]
    
    out = np.zeros((n_merger+1, n_snap), dtype=MERGER_DTYPE)
    
    for i in range(n_merger + 1):
        id = array.array("i")
        x, v = array.array("f"), array.array("f")
        mvir, vmax = array.array("f"), array.array("f")

        id.fromfile(f, n_snap)
        mvir.fromfile(f, n_snap)
        vmax.fromfile(f, n_snap)
        x.fromfile(f, n_snap*3)
        v.fromfile(f, n_snap*3)

        out["id"][i,:] = np.array(id, dtype=int)
        out["mvir"][i,:] = np.array(mvir, dtype=float)
        out["vmax"][i,:] = np.array(vmax, dtype=float)
        out["x"][i,:,:] = np.array(x, dtype=float).reshape((n_snap, 3))
        out["v"][i,:,:] = np.array(v, dtype=float).reshape((n_snap, 3))

    f.close()

    return out

def main():
    h = read_mergers(sys.argv[1])
    print(h["mvir"][0,:])
    print(h["mvir"][2,:])

if __name__ == "__main__": main()
