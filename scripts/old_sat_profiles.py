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

def p_disk_disrupt(r, ok):
    # The differential correction Mpeak > 8e8 Msun from Samuel et al. 2020
    a, d0, d1 = 0.8, 8, 78

    p = np.zeros(len(r))
    p_ok = r > d0
    p[ok & p_ok] = a*(1 - np.exp(-(r[ok & p_ok] - d0)/d1))

    return p

def main():
    palette.configure(False)

    proj = Projector(1000000)
    
    fig, ax = plt.subplots()

    ax.plot(proj.r_proj_r, proj.f, "o", c="k")
    ax.set_ylabel(r"$N/N_{\rm tot}$")
    ax.set_xlabel(r"$r_{\rm proj}/r_{\rm 3D}$")

    fig.savefig("../plots/core_plots/projection_kernel.png")

    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    n_hosts = symlib.n_hosts("SymphonyMilkyWay")
    param = symlib.simulation_parameters("SymphonyMilkyWay")
    invalid_hosts = [6, 9, 10, 16, 17, 31, 36, 37, 40, 42, 43]
    n_hosts_sim = 0

    weights_h, weights_c = [], []
    r_proj_h, r_proj_c = [], []
    mpeak_h, mpeak_c = [], []
    for i in range(n_hosts):
        if i in invalid_hosts: continue
        sim_dir = symlib.get_host_directory(base_dir, "SymphonyMilkyWay", i)
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
            if h_ok[i,-1]:
                weights_h.append(proj.f*p_h[i])
                r_proj_h.append(proj.r_proj_r*r_h[i])
                mpeak_h.append(np.ones(len(proj.f)) *
                               hist["mpeak"][i]/h["mvir"][0,-1])

        n_hosts_sim += 1

    weights_c = np.hstack(weights_c) / n_hosts_sim
    weights_h = np.hstack(weights_h) / n_hosts_sim
    r_proj_c = np.hstack(r_proj_c)
    r_proj_h = np.hstack(r_proj_h)

    # Make them psuedo-magnitude gaps
    dM_c = -2.5*np.log10(np.hstack(mpeak_c))
    dM_h = -2.5*np.log10(np.hstack(mpeak_h))

    elves_host = observations.read_ELVES_hosts()
    elves = observations.read_ELVES_sats()
    saga = observations.read_SAGA()
    saga = saga[saga["M_r"] < SAGA_LIMIT]

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

    N_saga, edges_saga = np.histogram(
        dM_r_saga[saga_ok], bins=bins_saga,
        weights=100/saga["spec_coverage"][saga_ok]/saga_n_host
    )
    N_elves, edges_elves = np.histogram(
        dM_V_elves[elves_ok], bins=bins_elves,
        weights=elves["p_sat"][elves_ok]/elves_n_host
    )
    N_c, edges_c = np.histogram(dM_c[c_ok], bins=400, weights=weights_c[c_ok])
    N_h, edges_h = np.histogram(dM_h[h_ok], bins=400, weights=weights_h[h_ok])

    N_saga, N_elves = np.cumsum(N_saga), np.cumsum(N_elves)
    N_c, N_h = np.cumsum(N_c), np.cumsum(N_h)

    ax.plot(edges_saga[1:], N_saga, color=pc("r"), label=r"${\rm SAGA}$")
    ax.plot(edges_elves[1:], N_elves, color=pc("b"), label=r"${\rm ELVES}$")
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
    N_cut_elves = np.array([0.3, 1, 3, 10])
    i_cut_saga = np.searchsorted(N_saga, N_cut_saga)
    i_cut_elves = np.searchsorted(N_elves, N_cut_elves)
    print(N_saga[i_cut_saga - 1])
    print(N_elves[i_cut_elves - 1])
    i_cut_c = np.searchsorted(N_c, N_cut_elves)
    i_cut_h = np.searchsorted(N_h, N_cut_elves)
    dM_r_cut_saga = edges_saga[i_cut_saga + 1]
    dM_V_cut_elves = edges_elves[i_cut_elves + 1]
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

    print(dM_r_cut_saga)
    for i in range(len(dM_r_cut_saga)):
        ok = (dM_r_saga < dM_r_cut_saga[i]) & (saga["r_proj"] < R_MAX)
        ax.hist(saga["r_proj"][ok], cumulative=True, color=saga_colors[i],
                histtype="step", ls="--",
                weights=100/saga["spec_coverage"][ok]/saga_n_host,
                bins=bins, lw=3)

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
    ax.plot([], [], "--", c=pc("k"), label=r"${\rm SAGA}$")
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
