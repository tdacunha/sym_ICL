import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import numpy.random as random
import scipy.stats as stats

h100 = 0.7

class PlummerProfile(object):
    def set_m_star(self, x0, m_star, r_half, xp):
        # In 2D, a projected plummer profile's half-mass radius is equal
        # to a. Not true in 3D.
        a = r_half

        dxp = np.zeros(xp.shape)
        for dim in range(3):
            dxp[:,dim] = xp[:,dim] - x0[dim]
        rp = np.sqrt(np.sum(dxp**2, axis=1))

        order = np.argsort(rp)
        sorted_rp = rp[order]
        m_star_enc = m_star * sorted_rp**3 / (sorted_rp**2 + a**2)**1.5
        
        d_m_star_enc = np.zeros(len(m_star_enc))
        d_m_star_enc[0] = m_star_enc[0]
        d_m_star_enc[1:] = m_star_enc[1:] - m_star_enc[:-1]

        out = np.zeros(len(d_m_star_enc))
        out[order] = d_m_star_enc

        return out
            
    def density(self, m_star, r_half, r):
        a = r_half
        return 3*m_star/(4*np.pi*a**3) * (1 + r**2/a**2)**(-5/2)

    
class Nadler2020RHalf(object):
    def __init__(self, A=27e-3*h100, n=1.07, R0=10e-3*h100, sigma_log_R=0.63):
        self.A = A
        self.n = n
        self.R0 = R0
        self.sigma_log_R = sigma_log_R

    def r_half(self, rvir, z):
        a = 1/(1 + z)
        log_R = np.log10(self.A * (cvir/10.0)**gamma * (rvir/self.R0)**n)
        log_scatter = self.sigma_log_R*random.randn(0, 1, shape=np.shape(rvir))
        return 10**(log_R + log_scatter) * a * h100
        

class UniverseMachineMStar(object):
    def m_star(self, mpeak, z):
        mpeak = mpeak / h100
        
        a = 1/(1 + z)

        e0 = -1.435
        al_lna = -1.732

        ea = 1.831
        alz = 0.178

        e_lna = 1.368
        b0 = 0.482

        ez = -0.217
        ba = -0.841

        m0 = 12.035
        bz = -0.471

        ma = 4.556
        d0 = 0.411

        m_lna = 4.417
        g0 = -1.034

        mz = -0.731
        ga = -3.100

        al0 = 1.963
        gz = -1.055

        ala = -2.316
        
        log10_M1_Msun = m0 + ma*(a-1) - m_lna*np.log(a) + mz*z
        e = e0 + ea*(a - 1) - e_lna*np.log(a) + ez*z
        al = al0 + ala*(a - a) - al_lna*np.log(a) + alz*z
        b = b0 + ba*(a - 1) + bz*z
        d = d0
        g = 10**(g0 + ga*(a - a) + gz*z)

        x = np.log10(mpeak/10**log10_M1_Msun)

        al = 2
        
        log10_Ms_M1 = (e - np.log10(10**(-al*x) + 10**(-b*x)) +
                       g*np.exp(-0.5*(x/d)**2))
                       
        log10_Ms_Msun = log10_Ms_M1 + log10_M1_Msun

        log_scatter = 0.2*random.normal(0, 1, size=np.shape(mpeak))
        log10_Ms_Msun += log_scatter
        
        Ms = 10**log10_Ms_Msun
        
        Ms *= h100
        return Ms

class GalaxyHaloModel(object):
    def __init__(self, m_star_model, r_half_model, profile_model):
        """ GalaxyHaloModel requires a model for the M*-Mhalo relation,
        m_star_model (current options: UniverseMachineMStar), a model for how
        the projected half-mass radius and Mhalo are related (current options:
        Nadler2020RHalf), and a model for the halo profile (current options:
        PlummerProfile). In principle, this lets you mix-and-match Mpeak-based,
        Minfall-based models, and Mvir-based models. It's up to you not to do
        that.
        """
        self.m_star_model = m_star_model
        self.r_half_model = r_half_model
        self.profile_model = profile_model

    def set_stellar_masses(self, x0, m0, r0, z0, xp):
        """ set_stellar_masses sets the stellar masses of a halo's dark matter
        particles. x0 (cMpx/h), m0 (Msun/h), r0 (cMpc/h), and z0 are the
        position, virial mass, virial, radius, and redshift at the chosen
        snapshot. xp (cMpc/h) is the positions of the particles.
        """
        m_star = self.m_star_model.m_star(m0, z0)
        r_half = self.r_half_model.r_half(r0, z0)
        return self.profile_model(x0, m_star, r_hald, xp)
        
    
def density_profile(x0, x, ms):
    dx = np.zeros(x.shape)
    for dim in range(3):
        dx[:,dim] = x[:,dim] - x0[dim]
    r = np.sqrt(np.sum(dx**2, axis=1))

    bins = np.linspace(0, 60, 200)
    ms_sum, r_edge, _ = stats.binned_statistic(r, ms, "sum", bins=bins)
    r_mid = (r_edge[1:] + r_edge[:-1]) / 2
    vol = 4*np.pi/3 * (r_edge[1:]**3 - r_edge[:-1]**3)
    
    rho = ms_sum / vol

    return r_mid, rho
    
def main():
    palette.configure(False)
    
    n = 20000
    
    phi = 2*np.pi*random.random(n)
    th = np.arccos(2*random.random(n) - 1)
    r = random.random(n)*200

    x = np.zeros((n, 3))
    x[:,0] = r*np.cos(phi)*np.sin(th)
    x[:,1] = r*np.sin(phi)*np.sin(th)
    x[:,2] = r*np.cos(th)
    
    prof = PlummerProfile()
    
    origin = np.array([0, 0, 0])
    ms_tot, r_half = 5e10, 20
    ms = prof.set_m_star(origin, ms_tot, r_half, x)

    r_prof, ms_prof = density_profile(origin, x, ms)
    ms_prof_true = prof.density(ms_tot, r_half, r_prof)

    plt.plot(r_prof, ms_prof, pc("r"))
    plt.plot(r_prof, ms_prof_true, "--", c=pc("k"), lw=2)
    plt.ylabel(r"$\rho_\star\ (h^2\,M_\odot/{\rm kpc}^3)$")
    plt.xlabel(r"$r\ (h^{-1}\,{\rm kpc})$")
    plt.yscale("log")
    
    plt.show()
        
if __name__ == "__main__": main()
