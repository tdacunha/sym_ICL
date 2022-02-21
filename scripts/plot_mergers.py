import read_mergers
import matplotlib.pyplot as plt
import numpy as np
import palette
from palette import pc

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

def plot_circle(ax, x, y, r, lw=1.5, c="k"):
    theta = np.linspace(0, 2*np.pi, 200)
    xc = r*np.sin(theta)
    yc = r*np.cos(theta)

    ax.plot(xc+x, yc+y, lw=lw, c=c)

def main():
    palette.configure(False)

    fig, ax = plt.subplots(figsize=(8.04, 8))
    
    m = read_mergers.read_mergers("../tmp_data/Halo983_merger_info.dat")
    a_start, a_end = 1/(20), 1
    a = 10**np.linspace(a_start, a_end, 236)

    for i in range(0, 236):
        plot_haloes(m[:,i], a[i], i, ax)
    #plt.show()
        
def plot_haloes(m, a, snap, ax):
    ok = m["x"][:,0] > 0
    rvir = mvir_to_rvir(m["mvir"], a, 0.27)

    ax.clear()
    
    ax.set_xlim(-2.1, 2.1)
    ax.set_ylim(-2.1, 2.1)
    x0, y0 = m["x"][0, 0], m["x"][0,1]
    
    for i in range(len(rvir)):
        if not ok[i]: continue
        c = pc("b")
        if i == 0:
            c = pc("r")
        plot_circle(ax, m["x"][i,0]-x0, m["x"][i,1]-y0, rvir[i], c=c)

    plt.savefig("../plots/merger_frames/merger_%d.png" % snap)
        
if __name__ == "__main__": main()
