import numpy as np
import matplotlib.pyplot as plt
import symlib
import subfind
import matplotlib.colors as mpl_colors
from colossus.cosmology import cosmology
from colossus.halo import mass_so
import os.path as path
import scipy.signal as signal
import scipy.interpolate as interpolate
import gravitree

def delta(x, x0):
    dx = np.zeros(x.shape)
    for dim in range(3):
        dx[:,dim] = x[:,dim] - x0[dim]
    return dx

def distance(x, x0):
    return np.sqrt(np.sum(delta(x, x0)**2, axis=1))

def n_most_bound(xc, vc, x, v, ok, n_core, param):
    if np.sum(ok) < n_core:
        return np.ones(n_core)*-1

    dx, dv = delta(x, xc), delta(v, vc)
    
    mp, eps = param["mp"]/param["h100"], param["eps"]/param["h100"]
    # We need to do exactly one unbinding pass here because we're only using
    # smoothly accreted particles.
    E_ok = gravitree.binding_energy(dx[ok], dv[ok], mp, eps, n_iter=1)
    E = np.inf*np.ones(len(ok))
    E[ok] = E_ok

    part = np.argpartition(E, n_core)[:n_core]
    part_order = np.argsort(E[part])
    return part[part_order]

"""
def is_bound(param, dx, dv, ok=None, order=None):
    rmax, vmax, pe, order = symlib.profile_info(param, dx, ok, order)
    pe *= vmax**2
    ke = 0.5*np.sum(dv**2, axis=1)

    if ok is None:
        return pe + ke < 0, order
    else:
        return (pe + ke < 0) & ok, order

def is_bound_iter(n_iter, param, dx, dv, ok=None, order=None):
    ok, order = is_bound(param, dx, dv, ok, order)
    n_bound = np.sum(ok)

    for i in range(n_iter - 1):
        ok, order = is_bound(param, dx, dv, ok, order)
        n_bound_i = np.sum(ok)
        if n_bound_i == n_bound: break
        
    return ok, order
"""

def rockstar_cores(snap_info, h, sub_idxs, n_core):
    core_idxs = np.ones((len(h), n_core), dtype=int)*-1
    param = snap_info.param

    snap = snap_info.snap
    
    for i_sub in sub_idxs:
        x, v = snap_info.x[i_sub], snap_info.v[i_sub]
        valid, owner =  snap_info.valid[i_sub], snap_info.owner[i_sub]
        # Note: we aren't actually unbinding the whole halo, just findng the
        # most bound of the smoothly accreted particles.
        ok = valid & (owner == 0)
        
        xc, vc = h["x"][i_sub,snap], h["v"][i_sub,snap]
        core_idxs[i_sub,:] = n_most_bound(xc, vc, x, v, ok, n_core, param)
        
    return core_idxs

class SnapshotData(object):
    def __init__(self, info, sim_dir, snap, a, h_cmov, param,
                 include_false_selections=False):
        ifs = include_false_selections
        self.x = symlib.read_particles(info, sim_dir, snap, "x",
                                       include_false_selections=ifs)
        self.v = symlib.read_particles(info, sim_dir, snap, "v",
                                       include_false_selections=ifs)
        self.valid = symlib.read_particles(info, sim_dir, snap, "valid",
                                           include_false_selections=ifs)
        self.owner = symlib.read_particles(info, sim_dir, snap, "ownership",
                                           include_false_selections=ifs)
        try:
            self.infall_cores = symlib.read_particles(
                info, sim_dir, snap, "infall_core",
                include_false_selections=ifs)
        except:
            self.infall_cores = None

        self.ok = [None]*len(self.x)
        
        # Center on the host halo, not the subhalo itself.
        for i in range(len(self.x)):
            self.x[i] = symlib.set_units_x(self.x[i], h_cmov[0,snap], a, param)
            self.v[i] = symlib.set_units_v(self.v[i], h_cmov[0,snap], a, param)
            self.ok[i] = self.valid[i] & (self.owner[i] == 0)
            
        self.mp, self.eps = symlib.set_units_parameters(a, param)
        self.snap = snap
        self.param = param

class SubhaloTrack(object):
    def __init__(self, i_sub, snap_data_init, core, param, k):
        x, v = snap_data_init.x[i_sub], snap_data_init.v[i_sub]
        
        self.i_sub = i_sub

        snap = snap_data_init.snap
        n_snap = param["n_snap"]
        n_core = len(core)
        self.n_core = n_core
        self.snaps = np.arange(n_snap, dtype=int)
        self.cores = np.ones((n_snap, n_core), dtype=int)*-1
        self.x = np.ones((n_snap, 3))*np.nan
        self.v = np.ones((n_snap, 3))*np.nan
        self.param = param
        self.k = k

        self._set_xv(snap_data_init, core)
        
    def next_snap(self, snap_data, core):    
        snap = snap_data.snap
        self._set_xv(snap_data, core)

    def _set_xv(self, snap_data, prev_core):
        snap = snap_data.snap
    
        x, v = snap_data.x[self.i_sub], snap_data.v[self.i_sub]
        valid, owner = snap_data.valid[self.i_sub], snap_data.owner[self.i_sub]
        mp, eps = snap_data.mp, snap_data.eps
        
        ok = valid & (owner == 0)

        x_sf, v_sf, rho_sf, owner_sf, peak_n = subfind.subfind(
            x[ok], v[ok], mp, k=self.k)
        owner_all = np.ones(len(x), dtype=int)*-1
        owner_all[ok] = owner_sf

        owner_votes = owner_all[prev_core]

        assert(-1 not in owner_votes)
        i_owner = np.argmax(np.bincount(owner_votes))
        xc_sf, vc_sf = x_sf[i_owner], v_sf[i_owner]

        # This isn't used any more and will just slow us down.
        #core_sf = n_most_bound(xc_sf, vc_sf, x, v, ok, self.n_core, self.param)
        
        self.cores[snap] = prev_core
        self.x[snap] = xc_sf
        self.v[snap] = vc_sf   
        
class MassProfile(object):
    # Positions should already be in physical kpc, no h100.
    def __init__(self, param, snap, h, x, own, valid):
        r_low = param["eps"] / 10
        r_high = h[0,snap]["rvir"]*10
        r_edges = 10**np.linspace(np.log10(r_low), np.log10(r_high), 100)

        h0 = h[0,snap]["x"]
        dm_hist = np.zeros(len(r_edges) - 1)
        
        for i in range(len(x)):
            dx = np.copy(x[i])
            for dim in range(3): dx[:,dim] -= h0[dim]
            r = np.sqrt(np.sum(dx**2, axis=1))

            ok = (own[i] == 0) & (valid[i])
            n, _, = np.histogram(r[ok], bins=r_edges)
            dm_hist += n*param["mp"]
            dm_hist[0] += np.sum(r < r_low)

        dlog_r = np.log10(r_edges[1]) - np.log10(r_edges[0])
        mr = np.cumsum(dm_hist)
        
        ok = (mr > 0) & (r_edges[1:] > param["eps"])
        dlog_m_dlog_r = signal.savgol_filter(
            np.log10(mr[ok]), 25, 4, deriv=1,
            delta=dlog_r, mode="interp"
        )
        dlog_m_dlog_r[dlog_m_dlog_r < 0] = 0

        log_r = np.log10(r_edges[1:][ok])
        
        self._deriv_interp = interpolate.interp1d(
            log_r, dlog_m_dlog_r, kind="linear",
            fill_value=(dlog_m_dlog_r[0], dlog_m_dlog_r[-1]),
            bounds_error=False
        )

        self._m_interp = interpolate.interp1d(
            log_r, np.log10(mr[ok]), kind="linear",
            fill_value=(np.log10(mr[ok])[0], np.log10(mr[ok])[-1]),
            bounds_error=False
        )

        self.param = param
        self.mp = param["mp"]/param["h100"]
        self.eps = param["eps"]/param["h100"]

    def m(self, r):
        return 10**self._m_interp(np.log10(r))

    def dln_m_dln_r(self, r):
        return self._deriv_interp(np.log10(r))

    def _tidal_radius_iter(self, m, dx, dv=None, method="centrifugal"):
        """ Helper function for tidal_radius

        methods are:
        "jacobi" -  Eq. 2 from vdB & Ogiya 2018
        "radial" - Eq. 3 from vdB & Ogiya 2018
        "centrifugal" - Eq. 4 from vdB & Ogia 2018
        "circular" - Eq. 5 from vdB & Ogiya 2018
        "klypin99" - Eq. 6 from vdB & Ogiya 2018

        This is not vectorized
        """

        if dx.shape != (3,):
            raise ValueError("dx must be single 3-vector.")

        R = np.sqrt(np.sum(dx**2))
        dln_M_dln_R = self.dln_m_dln_r(R)
        M = self.m(R)
        
        if method == "jacobi":
            factor = 3
        elif method == "radial":
            factor = 2 - dln_M_dln_R
        elif method == "centrifugal":
            if dv is None:
                raise ValueError("dv must be set if method = 'centrifugal'")
            # In units of  (km/s * kpc)/kpc^2 = km/kpc/s:
            Omega = np.sqrt(np.sum(np.cross(dv, dx)**2))/R**2
            # In units of 1/Gyr (funny how closw those two numbers are...)
            Omega *= 1.022
            Omega_circ = 0.750 * (M/1e12)**0.5 * (R/200)**-1.5
            factor = 2 + (Omega/Omega_circ)**2 - dln_M_dln_R
        elif method == "circular":
            factor = 3 - dln_M_dln_R
        elif method == "klypin99":
            factor = 1
        else:
            raise ValueError("Unrecognized tidal radius method, '%s'" % method)

        if factor <= 0:
            return np.inf
        elif factor < 1:
            factor = 1.0
        
        return R*(m/(factor*M))**(1.0/3)

    def tidal_radius(self, dx_core, dx, dv_core=None, dv=None,
                     method="centrifugal", bound_only=True):
        if dx_core[0] == 0 and dx_core[1] == 0 and dx_core[2] == 0:
            return 0.0

        dx = delta(dx, dx_core)
        dv = delta(dv, dv_core)
        
        dr = np.sqrt(np.sum(dx**2, axis=1))
        order = np.argsort(dr)

        m = len(dr)*self.mp
        conv_limit = 0.001
        
        mp, eps = self.mp, self.eps

        while True:
            r_tidal = self._tidal_radius_iter(m, dx_core, dv_core,
                                              method=method)
            ok = dr < r_tidal
            n_tidal = np.sum(ok)
            if bound_only:
                is_bound = np.zeros(len(ok), dtype=bool)
                is_bound_ok = gravitree.binding_energy(
                    dx[ok], dv[ok], mp, eps, n_iter=2) < 0
                is_bound[ok] = is_bound_ok
            else:
                is_bound = ok
                
            m_tidal = np.sum(is_bound)*mp

            if m_tidal == 0:
                return 0.0, 0.0
            if np.abs(m/m_tidal - 1) < conv_limit:
                return r_tidal, m_tidal

            
            m = m_tidal

        assert(0)
    
        
def main():
    import palette
    from palette import pc

    palette.configure(False)
    
    base_dir = "/sdf/home/p/phil1/ZoomIns"
    suite_name = "SymphonyMilkyWay"
    sim_dir = path.join(base_dir, suite_name, "Halo023")
    
    param = symlib.parameter_table[suite_name]
    h, hist = symlib.read_subhalos(sim_dir)
    c = symlib.read_cores(sim_dir)
    h_cmov, hist_cmov = symlib.read_subhalos(sim_dir, comoving=True)

    info = symlib.ParticleInfo(sim_dir)

    scale = symlib.scale_factors(sim_dir)
    cosmo = cosmology.setCosmology("", symlib.colossus_parameters(param))
    
    targets = np.where(c["ok"][:,235])[0]
    targets = targets[targets >= 3] # I don't want to deal with the big ones
    
    snap = 235
    sd = SnapshotData(info, sim_dir, snap, scale[snap], h_cmov, param)
    prof = MassProfile(sd.param, snap, h, sd.x, sd.owner, sd.valid)
    
    print("""# 0 - subhalo index
# 1 - log10 M_bound (no tides)
# 2 - log10 M_tidal (tides computed relative to bound mass)
# 3 - log10 M_tidal (tides computed relative to all mass)
# 4 - log10 R_tidal (tides computed relative to bound mass)
# 5 - log10 R_tidal (tides computed relative to all mass)
    """)

    for i in range(len(targets)):
        i_sub = targets[i]
        ok = sd.ok[i_sub]
        r_tidal, m_tidal = prof.tidal_radius(
            c["x"][i_sub,snap], sd.x[i_sub][ok],
            c["v"][i_sub,snap], sd.v[i_sub][ok],
            bound_only=True
        )
        r_tidal_ub, m_tidal_ub = prof.tidal_radius(
            c["x"][i_sub,snap], sd.x[i_sub][ok],
            c["v"][i_sub,snap], sd.v[i_sub][ok],
            bound_only=False
        )
        
        if m_tidal > 0:
            print("%4d %.3f %.3f %.3f %.3f %.3f" % 
                  (i_sub, np.log10(c["m_bound"][i_sub,235]),
                   np.log10(m_tidal), np.log10(m_tidal_ub),
                   np.log10(r_tidal), np.log10(r_tidal_ub)))
        else:
            print("%4d %.3f   nil %.3f   nil %.3f" % 
                  (i_sub, np.log10(c["m_bound"][i_sub,235]),
                   np.log10(m_tidal_ub), np.log10(r_tidal_ub)))
        
if __name__ == "__main__": main()
