import numpy as np
import symlib
import os
import os.path as path

#suites = ["LMC", "MilkyWay", "Group", "LCluster", "Cluster"]
suites = ["MilkyWay"]

snap_flag = [0, 0, 0, 1, 1]
snap_sets = [[235, 203, 181, 149, 109],
             [199, 168, 146, 115, 75]]

for ii in range(len(suites)):
    config_name = "../configs/%s/config.txt" % suites[ii]
    n_file = np.loadtxt(config_name, dtype=int, usecols=(4,))
    snap_fmt, out_dirs = np.loadtxt(config_name, dtype=str, usecols=(5, 7)).T

    for i in range(len(snap_fmt)):
        out_dir = path.join(out_dirs[i], "full_snapshots")
        cmd = "mkdir %s" % out_dir
        print(cmd)
        os.system(cmd)
        for k in range(5):
            snap = snap_sets[snap_flag[ii]][k]
            for j in range(n_file[i]):
                in_name = snap_fmt[i] % (snap, j)
                cmd = "cp %s %s/" % (in_name, out_dir)
                #print(cmd)
                os.system(cmd)
