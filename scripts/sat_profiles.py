import numpy as np
import observations
import matplotlib.pyplot as plt
import palette
from palette import pc
import numpy.random as random

SAGA_LIMIT = -12.3
ELVES_LIMIT = -9
R_MAX = 300

class Projector(object):
    def __init__(self, samples):
        theta = np.arccos(2*random.random(samples) - 1)
        r_proj = np.sin(theta)
        
        n_bins = int(np.ceil(2*samples**(1/3.0)))
        N, edges = np.histogram(r_proj, range=(0, 1), bins=n_bins)
        self.f = N / np.sum(N)
        self.r_proj_r = (edges[1:] + edges[:-1]) / 2
        
def main():
    palette.configure(False)
    
    elves_host = observations.read_ELVES_hosts()
    elves = observations.read_ELVES_sats()
    saga = observations.read_SAGA()
    saga = saga[saga["M_r"] < SAGA_LIMIT]

    fig, ax = plt.subplots()

    # -20.96 is the average of the reported Mr's in SAGA 1.
    Mr_saga_host = -20.96*np.ones(len(saga))
    MV_elves_host = elves_host["M_V"][elves["host_idx"]]
    dM_r_saga = saga["M_r"] - Mr_saga_host
    dM_V_elves = elves["M_V"] - MV_elves_host

    elves_ok = ((elves_host["r_cover"][elves["host_idx"]] >= R_MAX) &
                (elves["r_proj"] < R_MAX))
    saga_ok = (saga["M_r"] < SAGA_LIMIT) & (saga["r_proj"] < R_MAX)

    elves_n_host = np.sum(elves_host["r_cover"] >= R_MAX)
    saga_n_host = 1 + np.max(saga["host_idx"])
    saga_dM_r_lim = np.min(SAGA_LIMIT - Mr_saga_host)
    elves_dM_V_lim = np.min(ELVES_LIMIT - MV_elves_host)
    saga_dM_r_lim_2 = np.min(-15 - Mr_saga_host)
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
    N_saga, N_elves = np.cumsum(N_saga), np.cumsum(N_elves)
    
    ax.plot(edges_saga[1:], N_saga, color=pc("r"), label=r"${\rm SAGA}$")
    ax.plot(edges_elves[1:], N_elves, color=pc("b"), label=r"${\rm ELVES}$")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Delta m = m_{\rm sat} - m_{\rm host}$")
    ax.set_ylabel(r"$N(<\Delta m)$")
    ax.legend(loc="lower right", fontsize=16)
    ylo, yhi = ax.get_ylim()
    ax.set_ylim(ylo, yhi)
    ax.plot([saga_dM_r_lim_2]*2, [ylo, yhi], "--", c=pc("r", 0.7))
    ax.set_xlim(0, elves_dM_V_lim)

    fig.savefig("../plots/core_tests/stellar_mass_func.png")
    
    fig, ax = plt.subplots()
    
    bins = np.linspace(15, 310, 200)

    N_cut_saga = [0.3, 1, 3]
    N_cut_elves = [0.3, 1, 3, 10]
    i_cut_saga = np.searchsorted(N_saga, N_cut_saga)
    i_cut_elves = np.searchsorted(N_elves, N_cut_elves)
    dM_r_cut_saga = edges_saga[i_cut_saga + 1]
    dM_V_cut_elves = edges_elves[i_cut_elves + 1]

    saga_colors = [pc("r", 0.3), pc("r", 0.5), pc("r", 0.7)]
    elves_colors = [pc("b", 0.2), pc("b", 0.4), pc("b", 0.6), pc("b", 0.8)]

    labels = [r"$N(<\Delta m) = 0.3$",
              r"$N(<\Delta m) = 1$",
              r"$N(<\Delta m) = 3$",
              r"${\rm ELVES};\ N(<\Delta m) = 10$"]
    
    r_weights = np.zeros(len(elves))
    for i in range(len(elves_host)):
        r_weights[elves["r_proj"] < elves_host["r_cover"][i]] += 1
    r_weights[r_weights == 0] = 1
    r_weights = np.max(r_weights)/r_weights

    for i in range(len(dM_r_cut_saga)):
        ok = (dM_r_saga < dM_r_cut_saga[i]) & (saga["r_proj"] < R_MAX)
        ax.hist(saga["r_proj"][ok], cumulative=True, color=saga_colors[i],
                histtype="step",
                weights=100/saga["spec_coverage"][ok]/saga_n_host,
                bins=bins, lw=3)

    for i in range(len(dM_V_cut_elves)-1, -1, -1):
        ok = ((dM_V_elves < dM_V_cut_elves[i]) & 
              (elves_host["r_cover"][elves["host_idx"]] >= R_MAX) &
              (elves["r_proj"] < R_MAX))
        ax.hist(elves["r_proj"][ok], cumulative=True, color=elves_colors[i],
                weights=elves["p_sat"][ok]/elves_n_host,
                histtype="step", bins=bins, lw=3, label=labels[i])

    ax.plot([], [], pc("r"), label=r"${\rm SAGA}$")
        
    ax.legend(loc="upper left", fontsize=16)
    
    ax.set_xlim(15, R_MAX)
    ax.set_yscale("log")
    ax.set_xlabel(r"$r\ ({\rm kpc})$")
    ax.set_ylabel(r"$N(<r)$")
    ax.set_ylim(None, 30)

    fig.savefig("../plots/core_tests/stellar_radial_dist.png")
    
    
if __name__ == "__main__": main()
