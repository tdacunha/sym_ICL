import orbit_model
import struct
import array
import numpy as np
import galpy.potential as potential
import scipy.stats as stats
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import palette
import numpy.random as random
import lib
from palette import pc
from colossus.cosmology import cosmology
import m_star_profile_test_3 as profile_test
import star_tagging
import scipy.optimize as optimize
from colossus.halo import mass_so

params = { 'flat': True, 'H0': 70.0, 'Om0': 0.286,
           'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.95 }
f_baryon = params["Ob0"]/params["Om0"]
cosmo = cosmology.setCosmology('myCosmo', params)

QUANTILES = 10**np.linspace(-4, 0, 20)

P_FILE_FMT = "%s/particles/part_%03d.%d"
MP = 2.8e5
EPS = 170e-6

MIN_TRACERS = 5

def ranked_np_profile_matrix(ranks, idx, r, bins):
    n_rank = 1 + np.max(ranks)
    n_bins = len(bins)

    M = np.zeros((n_bins-1, n_rank))

    for i in range(n_rank):
        ri = r[ranks[idx] == i]
        N, _ = np.histogram(ri, bins=bins)
        #M[:,i] =  np.cumsum(N)
        M[:,i] = N
        
    return M

def cvir_nfw(vmax, mvir, rvir):
    vvir = 655.8 * (mvir/1e14)**0.5 * (rvir/1.0)**-0.5
    cv = vmax/vvir
    return orbit_model.cv_to_c_nfw(cv)

def fit_mp_star(M, m_enc):
    n = M.shape[1]
    res = optimize.lsq_linear(
        M, m_enc, bounds=(np.zeros(n), np.inf*np.ones(n))
    )
    return res.x
    

def m_half_from_prof(dm_prof, r):
    m_enc_prof = np.cumsum(dm_prof)
    return r[np.searchsorted(m_enc_prof, m_enc_prof[-1]/2)]

def main():
    random.seed(1)
    
    palette.configure(False) 
    #base_dir="/oak/stanford/orgs/kipac/users/phil1/simulations/MWest/Halo169/"
    base_dir = "../tmp_data/Halo169"
    
    #p_snaps = np.arange(127, 227 + 1)
    #p_snaps = np.arange(177, 227 + 1)
    p_snaps = [127, 177, 227]
    h_idx = 4
    host_id = 169
    
    _, m = lib.read_mergers(base_dir)

    scale = lib.scale_factors()
    z = 1/scale - 1
    age = cosmo.age(z)
    t_orbit = mass_so.dynamicalTime(z, "vir", "orbit")

    snap = p_snaps[0]
    ranks = [
        profile_test.rank_by_radius(base_dir, snap, h_idx),
        profile_test.rank_by_radial_energy(base_dir, snap, h_idx),
        profile_test.rank_by_nfw_energy(base_dir, snap, h_idx),
        profile_test.rank_by_infall(base_dir, snap, h_idx)
    ]
    
    r_bins = np.zeros(102)
    r_bins[1:] = 10**np.linspace(-2.5, 0, 101)

    comp_snaps = [127, 177, 227]

    rvir = None
    Nr_orig = None
    for snap in comp_snaps:
        rank0 = ranks[2]
        dt = (age[snap] - age[comp_snaps[0]])
        dt_orbit = dt/t_orbit[comp_snaps[0]]
        dt_string = ((r"$\Delta t = %.1f\ {\rm Gyr} = " +
                      r"%.1f\cdot t_{\rm orbit,vir}$") %
                     (dt, dt_orbit))
        plt.figure()
        if rvir is None:
            m_snap = snap
            
            m_star_enc = np.zeros((5, len(r_bins)-1))
            mp_star = np.zeros((5, len(QUANTILES)))
            
            profile_model = star_tagging.PlummerProfile()
            r_half_model = star_tagging.Jiang2019RHalf()
            m_star_model = star_tagging.UniverseMachineMStar()

            z = 1/scale[m_snap] - 1
            rvir = scale[m_snap]*m[h_idx,m_snap]["rvir"]
            cvir = cvir_nfw(m[h_idx,m_snap]["vmax"],
                            m[h_idx,m_snap]["mvir"], rvir)

            x, _, _ = lib.read_part_file(base_dir, snap, h_idx, ["x"])
            x, _, _, idx = profile_test.clean_particles(
                x,None,None, m[h_idx,snap], scale[snap])
            r = np.sqrt(np.sum(x**2, axis=1))
            
            M = ranked_np_profile_matrix(rank0, idx, r/rvir, r_bins)
            
            for i in range(5):
                r_half = r_half_model.r_half(rvir, cvir, z)
            
                mpeak = np.max(m[h_idx,:m_snap+1]["mvir"])
                m_star = m_star_model.m_star(mpeak, z)

                m_star_enc[i,:] = profile_model.m_enc(
                    m_star, r_half, r_bins[1:]*rvir)
                dm_star_enc = np.zeros(len(m_star_enc[i,:]))
                dm_star_enc[1:] = m_star_enc[i,1:] - m_star_enc[i,:-1]
                dm_star_enc[0] = m_star_enc[i,0]
                
                mp_star[i,:] = fit_mp_star(M, dm_star_enc)
        else:
            x, _, _ = lib.read_part_file(base_dir, snap, h_idx, ["x"])
            x, _, _, idx = profile_test.clean_particles(
                x,None,None, m[h_idx,snap], scale[snap])
            r = np.sqrt(np.sum(x**2, axis=1))
        
            M = ranked_np_profile_matrix(rank0, idx, r/rvir, r_bins)

                
        plt.plot(r_bins[1:], np.cumsum(np.sum(M*MP, axis=1)), c="k")

        if Nr_orig is None:
            Nr_orig = np.sum(M, axis=1)
        else:
            plt.plot(r_bins[1:], np.cumsum(Nr_orig*MP), "--", c="k")

        for i in range(M.shape[1]):
            brightness = 0.2 + 0.6*i/(M.shape[1]-1)
            plt.plot(r_bins[1:], np.cumsum(M[:,i]*MP),
                     pc("r", brightness), lw=2)

        plt.title(dt_string)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel(r"$r/R_{\rm vir,0}$")
        plt.ylabel(r"$M(<r)$")
        
        ylo, yhi = plt.ylim(MP/2, np.sum(Nr_orig)*MP*2)
        xlo, xhi = plt.xlim()
        plt.xlim(xlo, xhi)
        eps_ratio = 4*EPS*scale[snap]/rvir
        plt.fill_between([xlo, eps_ratio], [ylo]*2, [yhi]*2,
                         alpha=0.2, color="k")

        plt.figure()

        m_star_max = 0
        for i in range(m_star_enc.shape[0]):
            plt.plot(r_bins[1:], m_star_enc[i,:], pc("r"))
            mass_profile = M @ mp_star[i,:]
            plt.plot(r_bins[1:], np.cumsum(mass_profile),
                     "--", c=pc("b"), lw=2)
            if m_star_enc[i,-1] > m_star_max: m_star_max = m_star_enc[i,-1]

            r_half_true = m_half_from_prof(m_star_enc[i,:], r_bins[1:])
            r_half_tag = m_half_from_prof(np.cumsum(mass_profile), r_bins[1:])
            print("r_half,true/r_half,tag = %.3f" % (r_half_true/r_half_tag))

        print()

        plt.title(dt_string)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel(r"$r/R_{\rm vir,0}$")
        plt.ylabel(r"$M(<r)$")
        
        ylo, yhi = plt.ylim(MP/2, 2*m_star_max)
        plt.ylim(ylo, yhi)
        xlo, xhi = plt.xlim()
        plt.xlim(xlo, xhi)
        plt.fill_between([xlo, eps_ratio], [ylo]*2, [yhi]*2,
                         alpha=0.2, color="k")
        
    plt.show()

if __name__ == "__main__": main()
