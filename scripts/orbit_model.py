import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import galpy
import galpy.potential as potential
import galpy.orbit as orbit
import scipy.interpolate as interpolate
import scipy.signal as signal
import scipy.optimize as optimize
import lib

SHOW_T_R_CURVES = False

def _two_way_interpolate(x, y):
    """ two_way_interpolate returns interpolation functions that map from x -> y
    and y -> x. Both input arrays must monotonic.
    """
    if x[0] < x[1]:
        x_to_y = interpolate.interp1d(x, y)
    else:
        x_to_y = interpolate.interp1d(x[::-1], y[::-1])

    if y[0] < y[1]:
        y_to_x = interpolate.interp1d(y, x)
    else:
        y_to_x = interpolate.interp1d(y[::-1], x[::-1])

    return x_to_y, y_to_x

def _f(x): return np.log(1+x) - x/(1+x)
def _x_max_nfw(): return 2.1626
def _v_vmax_nfw(x): return 2.1506 * np.sqrt(_f(x) / x)
def _m_enc_nfw(x): return _f(x)
def _alpha_nfw(x): return -1 - 2*x/(1 + x)

def _half_mass_nfw(x, mass_fraction):
    def f_mass_ratio(xx):
        return _m_enc_nfw(xx) / _m_enc_nfw(x) - mass_fraction
    sol = optimize.root_scalar(f_mass_ratio, bracket=[1e-4*x, x])
    return sol.root

_c = 10**np.linspace(0, 3, 1000)
_cv = np.sqrt(_f(_x_max_nfw()) / _x_max_nfw() * _c/_f(_c))
# Make sure we're solving the correct part of the cV - cvir relation
_cv[_c < _x_max_nfw()] = 1.0 
c_to_cv_nfw, cv_to_c_nfw = _two_way_interpolate(_c, _cv)

def radius_hist(n_bins, r):
    i_max = signal.argrelextrema(r, np.greater)[0][0]
    i_min = signal.argrelextrema(r, np.less)[0][0]
    hist_range = (r[i_min], r[i_max])
    new_n_bins = int(round(n_bins / (r[i_max] - r[i_min])))
    if i_max < i_min: i_min, i_max = i_max, i_min
    
    n, r_edge = np.histogram(r[i_min: i_max+1], bins=n_bins, range=hist_range,
                             density=True)
    r_mid = (r_edge[1:] + r_edge[:-1]) / 2
    r_mean = np.mean(r[i_min: i_max+1])
    return r_mid, n, r_mean

def r_mean_dalal(o, pot):
    r_apo = o.rap(pot=pot, analytic=True)
    r_peri = o.rperi(pot=pot, analytic=True)
    return r_peri + (r_apo - r_peri)/np.sqrt(2)

def poly3(x, p0, p1, p2):
    return p0 + p1*x + p2*x**2

def poly2(x, p0, p1):
    return p0*x + p1*x*x

def r_avg_fit(r_apo, r_peri, r_half_apo_mass):
    r0 = r_peri/r_apo
    a_half = -np.log10(0.5)/np.log10(r_half_apo_mass/r_apo)
    C = -0.0425*a_half + 0.227
    D = -0.135*a_half + 2.108
    return 0.42 + 0.55*r0 + C*np.exp(-D*r0)

def rmax_nfw(vmax, mvir, rvir):
    # vmax must be in km/s, mvir in Msun and Rvir in pMpc. mvir and rivr must
    # either both have little-h or both not have it.
    vvir = 655.8 * (mvir/1e14)**0.5 * (rvir/1.0)**-0.5
    cv = vmax/vvir
    cvir = cv_to_c_nfw(cv)

    return rvir/cvir*_x_max_nfw()

def project(x, y, z, px_hat, py_hat, pz_hat):
    """ project projects the vectors given by (x, y, z) along the orthognoal
    unit vectors px_hat, py_hat, and pz_hat.
    """
    vec = np.array([x, y, z]).T
    return np.dot(vec, px_hat), np.dot(vec, py_hat), np.dot(vec, pz_hat)

def decompose_velocity(x, v):
    r = np.sqrt(np.sum(x**2))
    x_hat = x / r
    vr = np.sum(x_hat*v)

    vTot = np.sqrt(np.sum(v**2))
    vT = np.sqrt(vTot**2 - vr**2)

    return r, vr, vT

def peri_apo(rmax, vmax, x, v, vesc_func):
    param = np.zeros(6)
    r, vR, vT = decompose_velocity(x/rmax, v/vmax)
    
    v = np.sqrt(vR**2 + vT**2)
    vesc_r = vesc_func(r)
    if v >= vesc_r:
        return -1, -1, False

    L_m = r*vT
    e_m = (v**2 - vesc_r**2)/2
    
    def u_m_eff(rr):
        u_m = -vesc_func(rr)**2 / 2
        return L_m*L_m/(2*rr*rr) + u_m

    r_min, r_max = 1e-5, 100
    if u_m_eff(r_max) < e_m:
        return -1, -1, False    

    r_u_eff_min = optimize.minimize_scalar(u_m_eff, bracket=(r_min, r_max)).x
    f = lambda rr: u_m_eff(rr) - e_m
    
    r_apo = optimize.root_scalar(f, bracket=(r_u_eff_min, r_max)).root
    
    if u_m_eff(r_min) < e_m:
        r_peri = r_min
    else:
        r_peri = optimize.root_scalar(f, bracket=(r_min, r_u_eff_min)).root

    return r_peri*rmax, r_apo*rmax, True

def energy(rmax, vmax, x, v,
           pot=potential.NFWPotential(a=1/_x_max_nfw(), normalize=True)):
    r = np.sqrt(np.sum((x/rmax)**2, axis=1))
    v = np.sqrt(np.sum((v/vmax)**2, axis=1))    
    return 0.5*(v)**2, 0.5*pot.vesc(r)**2

def r_avg_nfw(rmax, r_peri, r_apo):
    r_half = _half_mass_nfw(r_apo/rmax*_x_max_nfw(), 0.5)/_x_max_nfw()
    r_avg = r_apo*r_avg_fit(r_apo/rmax, r_peri/rmax, r_half)
    return r_avg

def main():
    palette.configure(False)

    t = np.linspace(0, 10, 20001)

    vT = np.linspace(0, 1, 22)[1:-1]
    n_bins = 100
    #colors = [pc("r", c_min + dc*j) for j in range(len(vT))]

    r_rs = [1/1000, 1/100, 1/10, 1, 10, 100, 1000]
    #r_rs = [1]
    
    r_avg = np.zeros((len(r_rs), len(vT)))
    r_peri_r_apo = np.zeros((len(r_rs), len(vT)))
    
    for j in range(len(r_rs)):
        pot = potential.NFWPotential(normalize=1.0, a=r_rs[j])

        if SHOW_T_R_CURVES: plt.figure()
        
        for i in range(len(vT)):
            r_bins = np.linspace(0, 1, 101)

            o = orbit.Orbit([1, 0.00, vT[i]])
            o.integrate(t, pot, method='dopr54_c')

            r_mids, n, r_avg[j,i] = radius_hist(
                n_bins, o.r(t)/o.rap(pot=pot, analytic=True))

            r_peri_r_apo[j,i] = (o.rperi(pot=pot, analytic=True) /
                                 o.rap(pot=pot, analytic=True))
            
            ok = n>0
            if SHOW_T_R_CURVES:
                plt.plot(r_mids[ok], n[ok],  colors[i],
                         label=r"$$")
                plt.plot([r_avg[j,i]]*2, [0.1, 0.5], lw=1.5, c=colors[i])
        
        if SHOW_T_R_CURVES:
            #plt.xscale("log")
            #plt.yscale("log")
            plt.ylim(0.1, 40)
        
            plt.xlabel(r"$r/r_{\rm peri}$")
            plt.ylabel(r"$Pr(r/r_{\rm peri})\ ({\rm time-averaged})$")


    colors = [pc("k"), pc("r"), pc("o"), pc("g"), pc("b"), pc("p"), pc("a")] 
    def fit_func(x, p0, p1):
        return 0.42 + 0.55*x + p0*np.exp(-p1*x)

    #def fit_func()

    slopes = np.zeros(len(r_rs))
    p_opt = np.zeros((2, len(r_rs)))
    
    for j in range(len(r_rs)):
        r0 = r_peri_r_apo[j,:]
        r1 = r_avg[j,:]

        alpha_a = _alpha_nfw(r_rs[j])
        alpha_p = _alpha_nfw(r0*r_rs[j])
        mass_slope = np.log10(_m_enc_nfw(r_rs[j])/_m_enc_nfw(r_rs[j]*r0)) / np.log10(r0)
        x_half = _half_mass_nfw(r_rs[j], 0.5)
        mass_slope_2 = -np.log10(0.5) / np.log10(x_half/r_rs[j])
        slopes[j] = mass_slope_2
        plt.figure(0)
        plt.plot(r0, r1, c=colors[j])

        p_opt[:,j], _ = optimize.curve_fit(fit_func, r0, r1)
        #plt.figure(0)
        #plt.plot(r0, fit_func(r0, p_opt[0,j], p_opt[1,j]) - (0.42 + 0.55*r0),
         #        "--", c="k", lw=1.5)
        
        plt.figure(1)
        plt.plot(r0, mass_slope, "--", c=colors[j])
        plt.plot(r0, mass_slope_2*np.ones(len(r0)), c=colors[j])
        
        print("%.3f %.3f %.3f %.3f %.3f" % (np.mean(mass_slope), alpha_a,
                                            np.mean(alpha_p),
                                            p_opt[0,j], p_opt[1,j]))
        """
        plt.figure(1)
        plt.plot(alpha, [p_opt[0]], "o", c=colors[j])
        
        plt.figure(2)
        plt.plot(alpha, [p_opt[0]], "o", c=colors[j])
        """
        
    def linear(x, a, b): return a*x + b

    fit1, _ = optimize.curve_fit(linear, slopes, p_opt[0,:])
    fit2, _ = optimize.curve_fit(linear, slopes, p_opt[1,:])

    print(fit1, fit2)

    for j in range(len(slopes)):
        plt.figure(0)
        r0 = r_peri_r_apo[j,:]
        r1 = r_avg[j,:]
        x_half = _half_mass_nfw(r_rs[j], 0.5)

        fit = r_avg_fit(r_rs[j], r0*r_rs[j], x_half)
        
        plt.plot(r0, fit, "--", c="k", lw=1.5)
        print(np.mean(fit - r1),
              np.sqrt(np.mean((fit - r1)**2)))

    
    #plt.yscale("log")
    #plt.xscale("log")
        
    plt.show()

def main2():
    params = {'flat': True, 'H0': 70.0, 'Om0': 0.286, 'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.95}
    cosmo = cosmology.setCosmology('myCosmo', params)

    base_dir = "../tmp_data/Halo169/"
    p_snaps = [127, 177, 227]
    h_idx = 4
    p_files = [
        "../tmp_data/Halo169/particles/part_127.4",
        "../tmp_data/Halo169/particles/part_177.4",
        "../tmp_data/Halo169/particles/part_227.4"
    ]

    _, m = lib.read_mergers(base_dir)

    scale = lib.scale_factors()
    z = 1/scale - 1
    
    for i_snap in range(len(p_snaps)):
        snap = p_snaps[i_snap]
        rvir, mvir = m["rvir"][h_idx, snap], m["mvir"][h_idx, snap]
        rvir *= scale[snap]
        
        rvir_colossus = mass_so.M_to_R(mvir, z[snap], "vir")
        
        vmax = m["vmax"][h_idx, snap]

        rmax = rmax_nfw(vmax, mvir, rvir)
        
        print("%g %g %g %g" % (rvir, mvir, vmax, rmax))
    
if __name__ == "__main__": main()
