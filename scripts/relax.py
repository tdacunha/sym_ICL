import numpy as np
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology
from colossus.halo import mass_so
import palette
from palette import pc
import m_star_profile_test_3 as m_star_profile
import lib

MP = 2.8e5
EPS = 170e-6
QUANTILES = 10**np.linspace(-4, 0, 12)
h100 = 0.7

params = { 'flat': True, 'H0': 70.0, 'Om0': 0.286,
           'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.95 }
cosmo = cosmology.setCosmology('myCosmo', params)

def t_relax_t_orbit(Nr, r, eps):
    return Nr/4 * (np.log(r**2/eps**2 + 1) +
                   eps**2-2*r**2/(3*(eps**2+r**2)) -
                   np.log(3/2))**-1

def t_orbit(Mr, r):
    return 7.5 * (r/40)**1.5 * (1e10/Mr)**0.5

def t_relax_profile(mp, r, eps):
    # mp, r, and eps are in Msun, pkpc, and kpc, repsectively. Answer is in Gyr.
    order = np.argsort(r)
    Nr = np.zeros(len(r))
    Nr_sort = np.arange(len(Nr)) + 1
    Nr[order] = Nr_sort
    
    t_orb = t_orbit(Nr*mp, r)
    t_relax = t_relax_t_orbit(Nr, r, eps)*t_orb
    
    return t_relax
    
def main():
    STARTING_T_RELAX = True
    
    palette.configure(False)

    host_id = 169
    base_dir = "../tmp_data/Halo%03d" % host_id

    _, m = lib.read_mergers(base_dir)
    
    scale = lib.scale_factors()
    z = 1/scale - 1
    age = cosmo.age(z)
    t_orbit = mass_so.dynamicalTime(z, "vir", "orbit")
    
    h_idx, scale_start, scale_infall = {
        169: (4, 0.25, 0.9),
        222: (6, 0.3, 0.95),
        523: (10, 0.35, 0.9),
        719: (1, 0.4, 0.85)
    }[host_id]
    scale = lib.scale_factors()
    snap_start, snap_infall = 127, 227
    
    dt = (age - age[snap_start])
    dt_orbit = dt/t_orbit[snap_start]

    comp_snaps = np.searchsorted(dt_orbit, [0.1, 0.25, 0.5, 1.0, 2.0, 4.0])
    
    ranks = m_star_profile.rank_by_radial_energy(base_dir, snap_start, h_idx)

    snap0 = comp_snaps[0]
    _rr0, x0, idx0 = m_star_profile.rank_radii(
        m, base_dir, snap0, h_idx, [ranks], [0.5])    
    rr0 = _rr0[0,:,0]
    ok = rr0 >= 0

    if STARTING_T_RELAX:
        r = np.sqrt(np.sum(x0**2, axis=1))
        t_relax = t_relax_profile(
            MP/h100, r*1e3/h100, EPS*scale[snap_start]*1e3/h100)
        med_t_relax = np.zeros(len(rr0))
        for i in range(len(med_t_relax)):
            if not ok[i]:
                med_t_relax[i] = -1
            else:
                med_t_relax[i] = np.mean(t_relax[ranks[idx0] == i])

    print(med_t_relax)
                
    colors = [pc("r"), pc("o"), pc("g"), pc("b"), pc("p"), pc("k")]
    
    plt.figure(0)
    for k in range(len(comp_snaps)):
        snap = comp_snaps[k]
        
        _rr, x, idx = m_star_profile.rank_radii(
            m, base_dir, snap, h_idx, [ranks], [0.5])
        rr = _rr[0,:,0]

        if not STARTING_T_RELAX:
            r = np.sqrt(np.sum(x**2, axis=1))
            t_relax = t_relax_profile(
                MP/h100, r/h100*1e3, 1e3*EPS*scale[snap]/h100)
            med_t_relax = np.zeros(len(rr))
            for i in range(len(med_t_relax)):
                if not ok[i]:
                    med_t_relax[i] = -1
                else:
                    med_t_relax[i] = np.median(t_relax[ranks[idx] == i])

        plt.plot(dt[snap]/med_t_relax[ok], rr[ok]/rr0[ok], c=colors[k],
                 label=r"$\Delta t/t_{\rm orbit,vir}$ = %.2f" % dt_orbit[snap])
        plt.plot(dt[snap]/med_t_relax[ok], rr[ok]/rr0[ok], "o", c=colors[k])

        
    ylo, yhi = plt.ylim()
    plt.ylim(ylo, yhi)
    plt.plot([0.177, 0.177], [ylo, yhi], "--", c="k")
        
    plt.legend(loc="upper left", fontsize=16)
    plt.xscale("log")
    plt.ylabel(r"$R_{1/2}/R_{1/2,0}$")
    plt.xlabel(r"$\Delta t/t_{\rm relax}$")
    
    plt.savefig(("../plots/m_star_profs/profile_evolution/" +
                 "relaxation.%d.%d.png") % (host_id, h_idx))

if __name__ == "__main__": main()
