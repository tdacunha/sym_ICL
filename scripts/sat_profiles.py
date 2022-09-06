import numpy as np
import observations
import matplotlib.pyplot as plt
import palette
from palette import pc
import numpy.random as random
import symlib
import core_tests

SAGA_LIMIT = -12.3
ELVES_LIMIT = -9
R_MAX = 270

DISK_DISRUPTION=True

class Projector(object):
    def __init__(self, samples):
        theta = np.arccos(2*random.random(samples) - 1)
        r_proj = np.sin(theta)
        
        n_bins = int(np.ceil(2*samples**(1/3.0)))
        N, edges = np.histogram(r_proj, range=(0, 1), bins=n_bins)
        self.f = N / np.sum(N)
        self.r_proj_r = (edges[1:] + edges[:-1]) / 2

def project(x, y, z, px_hat, py_hat, pz_hat):
    """ project projects the vectors given by (x, y, z) along the orthogonal
    unit vectors px_hat, py_hat, and pz_hat.
    """
    vec = np.array([x, y, z]).T
    return np.dot(vec, px_hat), np.dot(vec, py_hat), np.dot(vec, pz_hat)

def axes(pz, incline):
    """ axes returns three orthogonal unit vectors, pz_hat, py_hat, and pz_hat.
    pz_hat will point in the direction of pz, px_hat will be orthogonal to 
    both pz and incline, and py_hat will be orthogonal to both of them.
    In the case where pz and incline or orthogonal, py_hat will be in the same
    direciton as incline.
    """
    pz_hat = pz / np.sqrt(pz[0]**2 + pz[1]**2 + pz[2]**2)
    px_hat = np.cross(incline, pz_hat)
    px_hat /= np.sqrt(px_hat[0]**2 + px_hat[1]**2 + px_hat[2]**2)
    py_hat = np.cross(px_hat, pz_hat)
    return px_hat, py_hat, pz_hat

def random_vector():
    """ random_vec returns a unit vector pointing in a random direction.
    """
    phi = 2 * np.pi * random.rand(1)[0]
    theta = np.arccos(2*random.rand(1)[0] - 1)

    x = np.sin(theta)*np.cos(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(theta)
    
def random_view(x, v):
    px_hat, py_hat, pz_hat = axes(random_vector(), random_vector())
    x1, x2, x3 = project(x[:,0], x[:,1], x[:,2], px_hat, py_hat, pz_hat)
    v1, v2, v3 = project(v[:,0], v[:,1], v[:,2], px_hat, py_hat, pz_hat)

    v_los = v2
    r_2d = np.sqrt(x1*x1 + x2*x2)
    return r_2d, v_los

def p_disk_disrupt(r, ok):
    # The differential correction Mpeak > 8e8 Msun from Samuel et al. 2020
    a, d0, d1 = 0.8, 8, 78

    p = np.zeros(len(r))
    p_ok = r > d0
    p[ok & p_ok] = a*(1 - np.exp(-(r[ok & p_ok] - d0)/d1))

    return p

def jackknife_histogram(x, bins, weights, flags):
    flag_max = int(np.max(flags))

    N_mean, _ = np.histogram(x, bins=bins, weights=weights)
    N_mean = np.cumsum(N_mean)

    N_all = []
    for flag in range(flag_max + 1):
        ok = flags != flag
        N, _ = np.histogram(x[ok], bins=bins, weights=weights[ok])
        N_all.append(np.cumsum(N))

    N_err = np.std(N_all, axis=0) * np.sqrt(flag_max + 1)
    N_log_err = np.std(np.log10(N_all), axis=0) * np.sqrt(flag_max + 1)
    return N_mean, N_err, N_log_err

def main():
    palette.configure(False)

    elves_host = observations.read_ELVES_hosts()
    elves = observations.read_ELVES_sats()
    saga = observations.read_SAGA()
    saga = saga[saga["M_r"] < SAGA_LIMIT]
    
    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    n_hosts = symlib.n_hosts("SymphonyMilkyWay")
    param = symlib.simulation_parameters("SymphonyMilkyWay")
    invalid_hosts = [6, 9, 10, 16, 17, 31, 36, 37, 40, 42, 43]
    n_hosts_sim = 0

    proj = Projector(1000000)

    weights_h, weights_c = [], []
    flags_h, flags_c = [], []
    r_proj_h, r_proj_c = [], []
    mpeak_h, mpeak_c = [], []
    for i_host in range(n_hosts):
        if i_host in invalid_hosts: continue
        sim_dir = symlib.get_host_directory(
            base_dir, "SymphonyMilkyWay", i_host)
        h, hist = symlib.read_subhalos(sim_dir)
        c = symlib.read_cores(sim_dir)

        c_ok = core_tests.core_ok(c, param)
        h_ok = core_tests.rockstar_ok(h, c, param)
        h_ok[0] = False

        r_c = np.sqrt(np.sum(c["x"]**2, axis=2))[:,-1]
        r_h = np.sqrt(np.sum(h["x"]**2, axis=2))[:,-1]

        if DISK_DISRUPTION:
            p_c = p_disk_disrupt(r_c, c_ok[:,-1])
            p_h = p_disk_disrupt(r_h, h_ok[:,-1])
        else:
            p_c, p_h = 1.0, 1.0

        for i in range(len(h)):
            if c_ok[i,-1]:
                weights_c.append(proj.f*p_c[i])
                r_proj_c.append(proj.r_proj_r*r_c[i])
                mpeak_c.append(np.ones(len(proj.f)) * 
                               hist["mpeak"][i]/h["mvir"][0,-1])
                flags_c.append(np.ones(len(proj.f))*i_host)
            if h_ok[i,-1]:
                weights_h.append(proj.f*p_h[i])
                r_proj_h.append(proj.r_proj_r*r_h[i])
                mpeak_h.append(np.ones(len(proj.f)) *
                               hist["mpeak"][i]/h["mvir"][0,-1])
                flags_h.append(np.ones(len(proj.f))*i_host)

        n_hosts_sim += 1

    weights_c = np.hstack(weights_c) / n_hosts_sim
    weights_h = np.hstack(weights_h) / n_hosts_sim
    r_proj_c = np.hstack(r_proj_c)
    r_proj_h = np.hstack(r_proj_h)
    flags_c = np.hstack(flags_c)
    flags_h = np.hstack(flags_h)


    # Make them psuedo-magnitude gaps
    dM_c = -2.5*np.log10(np.hstack(mpeak_c))
    dM_h = -2.5*np.log10(np.hstack(mpeak_h))

    fig, ax = plt.subplots()

    # -20.96 is the average of the reported Mr's in SAGA 1.
    MV_elves_host = elves_host["M_V"][elves["host_idx"]]
    dM_r_saga = saga["M_r"] - saga["host_M_r"]
    dM_V_elves = elves["M_V"] - MV_elves_host

    elves_ok = ((elves_host["r_cover"][elves["host_idx"]] >= R_MAX) &
                (elves["r_proj"] < R_MAX))
    saga_ok = (saga["M_r"] < SAGA_LIMIT) & (saga["r_proj"] < R_MAX)
    c_ok = r_proj_c < R_MAX
    h_ok = r_proj_h < R_MAX

    elves_n_host = np.sum(elves_host["r_cover"] >= R_MAX)
    saga_n_host = 1 + np.max(saga["host_idx"])
    saga_dM_r_lim = np.min(SAGA_LIMIT - saga["host_M_r"])
    elves_dM_V_lim = np.min(ELVES_LIMIT - MV_elves_host)
    bins_saga = np.linspace(0, saga_dM_r_lim, 400)
    bins_elves = np.linspace(0, elves_dM_V_lim, 400)

    print("jackknifing")
    weights_saga = 100/saga["spec_coverage"][saga_ok]/saga_n_host
    N_saga, N_err_saga, N_log_err_saga = jackknife_histogram(
        dM_r_saga[saga_ok], bins_saga, weights_saga,
        saga["host_idx"][saga_ok]
    )
    weights_elves = elves["p_sat"][elves_ok]/elves_n_host
    N_elves, N_err_elves, N_log_err_elves = jackknife_histogram(
        dM_V_elves[elves_ok], bins_elves, weights_elves,
        elves["host_idx"][elves_ok]
    )

    print(N_elves)
    print(N_err_elves)
    print(N_log_err_elves)

    edges_c = np.linspace(np.min(dM_c), np.max(dM_c), 400)

    N_c, N_err_c, N_log_err_c = jackknife_histogram(
        dM_c[c_ok], edges_c, weights_c[c_ok], flags_c[c_ok])
    edges_h = np.linspace(np.min(dM_h), np.max(dM_h), 400)
    N_h, N_err_h, N_log_err_h = jackknife_histogram(
        dM_h[h_ok], edges_h, weights_h[h_ok], flags_h[h_ok])
    print("done jackknifing")
    
    ax.plot(bins_saga[1:], N_saga, color=pc("r"), label=r"${\rm SAGA}$")
    high_saga = 10**(np.log10(N_saga) + N_log_err_saga)
    low_saga = 10**(np.log10(N_saga) - N_log_err_saga)
    ax.fill_between(bins_saga[1:], low_saga, high_saga,
                    color=pc("r"), alpha=0.2)
    ax.plot(bins_elves[1:], N_elves, color=pc("b"), label=r"${\rm ELVES}$")
    high_elves = 10**(np.log10(N_elves) + N_log_err_elves)
    low_elves = 10**(np.log10(N_elves) - N_log_err_elves)
    ax.fill_between(bins_elves[1:], low_elves, high_elves,
                    color=pc("b"), alpha=0.2)
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Delta m = m_{\rm sat} - m_{\rm host}$")
    ax.set_ylabel(r"$N(<\Delta m)$")
    ax.legend(loc="lower right", fontsize=16)
    ylo, yhi = ax.get_ylim()
    ax.set_ylim(ylo, yhi)
    ax.set_xlim(0, elves_dM_V_lim)

    fig.savefig("../plots/core_plots/stellar_mass_func.png")
    
    fig, ax = plt.subplots()
    
    bins = np.linspace(15, 310, 200)

    N_cut_saga = np.array([0.3, 1, 3])
    N_cut_elves = np.array([0.3, 1, 3, 8])
    i_cut_saga = np.searchsorted(N_saga, N_cut_saga)
    i_cut_elves = np.searchsorted(N_elves, N_cut_elves)

    i_cut_c = np.searchsorted(N_c, N_cut_elves)
    i_cut_h = np.searchsorted(N_h, N_cut_elves)
    dM_r_cut_saga = bins_saga[i_cut_saga + 1]
    dM_V_cut_elves = bins_elves[i_cut_elves + 1]
    dM_cut_c = edges_c[i_cut_c + 1]
    dM_cut_h = edges_h[i_cut_h + 1]

    print("dM_cut_c", dM_cut_c)
    print("dM_cut_h", dM_cut_h)
    print("dM_r_cut_saga", dM_r_cut_saga)
    print("dM_V_cut_elves", dM_V_cut_elves)

    saga_colors = [pc("k", 0.5), pc("k", 0.6), pc("k", 0.7)]
    elves_colors = [pc("k", 0.5), pc("k", 0.6), pc("k", 0.7), pc("k", 0.8)]
    h_colors = [pc("b", 0.2), pc("b", 0.4), pc("b", 0.6), pc("b", 0.8)]
    c_colors = [pc("r", 0.2), pc("r", 0.4), pc("r", 0.6), pc("r", 0.8)]

    r_weights = np.zeros(len(elves))
    for i in range(len(elves_host)):
        r_weights[elves["r_proj"] < elves_host["r_cover"][i]] += 1
    r_weights[r_weights == 0] = 1
    r_weights = np.max(r_weights)/r_weights

    """
    print(dM_r_cut_saga)
    for i in range(len(dM_r_cut_saga)):
        ok = (dM_r_saga < dM_r_cut_saga[i]) & (saga["r_proj"] < R_MAX)
        ax.hist(saga["r_proj"][ok], cumulative=True, color=saga_colors[i],
                histtype="step", ls="--",
                weights=100/saga["spec_coverage"][ok]/saga_n_host,
                bins=bins, lw=3)
    """

    for i in range(len(dM_V_cut_elves)):
        ok = ((dM_V_elves < dM_V_cut_elves[i]) & 
              (elves_host["r_cover"][elves["host_idx"]] >= R_MAX) &
              (elves["r_proj"] < R_MAX))
        ax.hist(elves["r_proj"][ok], cumulative=True, color=elves_colors[i],
                weights=elves["p_sat"][ok]/elves_n_host,
                histtype="step", bins=bins, lw=3)

        ok_c = (r_proj_c < R_MAX) & (dM_c < dM_cut_c[i])
        ax.hist(r_proj_c[ok_c], cumulative=True, color=c_colors[i],
                weights=weights_c[ok_c], histtype="step", bins=bins, lw=3)

        ok_h = (r_proj_h < R_MAX) & (dM_h < dM_cut_h[i])
        ax.hist(r_proj_h[ok_h], cumulative=True, color=h_colors[i],
                weights=weights_h[ok_h], histtype="step", bins=bins, lw=3)


    ax.plot([], [], pc("k"), label=r"${\rm ELVES}$")
    #ax.plot([], [], "--", c=pc("k"), label=r"${\rm SAGA}$")
    ax.plot([], [], pc("r"), label=r"${\rm core-tracking}$")
    ax.plot([], [], pc("b"), label=r"${\rm Rockstar}$")
        
    ax.legend(loc="upper left", fontsize=16)
    
    ax.set_xlim(15, R_MAX)
    ax.set_yscale("log")
    ax.set_xlabel(r"$r\ ({\rm kpc})$")
    ax.set_ylabel(r"$N(<r)$")
    ax.set_ylim(0.05, 30)

    if DISK_DISRUPTION:
        fig.savefig("../plots/core_plots/sat_radial_dist_disk.png")
    else:
        fig.savefig("../plots/core_plots/sat_radial_dist_no_disk.png")
    
    
if __name__ == "__main__": main()
