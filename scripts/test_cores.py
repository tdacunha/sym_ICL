import numpy as np
import symlib
import matplotlib.pyplot as plt
import os.path as path

def main():
    suite_name, halo_name = "SymphonyMilkyWay", "Halo023"
    base = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    sim_dir = path.join(base, suite_name, halo_name)

    param = symlib.parameter_table[suite_name]
    
    h, hist = symlib.read_subhalos(param, sim_dir)

    print(h["x"][20,-10:])
    print(h["x_core"][20,-10:])
    
if __name__ == "__main__": main()
