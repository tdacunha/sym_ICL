import numpy as np
import matplotlib.pyplot as plt
import symlib

try:
    import palette
    palette.configure(False)
except:
    pass

def main():
    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    suite = "SymphonyMilkyWay"

    param = symlib.simulation_parameters(suite)
    n_hosts = symlib.n_hosts(suite)

    macc_host_mpeak = []
    macc_1st_mpeak = []
    scale_1st = []
    scale_host = []
    scale_peak = []

    for i in range(n_hosts):
        sim_dir = symlib.get_host_directory(base_dir, suite, i)
        h, hist = symlib.read_subhalos(sim_dir)
        scale = symlib.scale_factors(sim_dir)

        host_idx = np.arange(len(h))
        macc_host_mpeak.append(h["mvir"][host_idx,hist["merger_snap"]]/
                               hist["mpeak"])
        macc_1st_mpeak.append(h["mvir"][host_idx,hist["first_infall_snap"]]/
                               hist["mpeak"])
        scale_1st.append(scale[hist["first_infall_snap"]])
        scale_host.append(scale[hist["merger_snap"]])
        scale_peak.append(scale[np.argmax(h["mvir"], axis=1)])

    macc_host_mpeak = np.hstack(macc_host_mpeak)
    macc_1st_mpeak = np.hstack(macc_1st_mpeak)
    scale_1st = np.hstack(scale_1st)
    scale_host = np.hstack(scale_host)
    scale_peak = np.hstack(scale_peak)

    print(np.mean(macc_host_mpeak))
    print(np.mean(macc_1st_mpeak))
    print(np.mean(scale_1st))
    print(np.mean(scale_host))
    print(np.mean(scale_peak))
    
    bins = np.logspace(-2, 0, 40)

    plt.hist(macc_host_mpeak, bins=bins, histtype="step", color="tab:red",
             label=r"${\rm host}$", density=True, lw=2)
    plt.hist(macc_1st_mpeak, bins=bins, histtype="step", color="tab:blue",
             label=r"${\rm first\ infall}$", density=True, lw=2)
    plt.xlabel(r"$M_{\rm acc}/M_{\rm peak}$")
    plt.ylabel(r"$P(M_{\rm acc}/M_{\rm peak})$")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc="upper left")

    plt.savefig("../plots/macc_mpeak.png")

if __name__ == "__main__": main()
