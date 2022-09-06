import symlib
import numpy as np

def main():
    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    suites = ["SymphonyLMC", "SymphonyMilkyWay", "SymphonyGroup", "SymphonyLCluster", "SymphonyCluster"]

    cuts = [1.0, 1.1, 2.0, 10.0]
    n_cut_merger = np.zeros(len(cuts), dtype=int)
    n_cut_infall = np.zeros(len(cuts), dtype=int)
    n_total = 0
    n_core_err = 0
    n_false_selection = 0
    n_false_selection_2 = 0

    for suite in suites:
        print(suite)

        param = symlib.simulation_parameters(suite)
        mp = param["mp"] / param["h100"]
        
        for i_host in range(symlib.n_hosts(suite)):
            print("    host", i_host)

            sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
            h, hist = symlib.read_subhalos(sim_dir)

            mpeak_infall = np.zeros(len(h))
            mpeak_merger = np.zeros(len(h))

            for i_sub in range(1, len(h)):
                infall_snap = hist["first_infall_snap"][i_sub]
                merger_snap = hist["merger_snap"][i_sub]
                mpeak_infall[i_sub] = np.max(h["mvir"][i_sub,:infall_snap+1])
                mpeak_merger[i_sub] = np.max(h["mvir"][i_sub,:merger_snap+1])

            mpeak_infall[0] = hist["mpeak"][0]
            mpeak_merger[0] = hist["mpeak"][0]

            n_total += len(h)
            n_core_err +=np.sum(mpeak_infall < 32*mp)
            n_false_selection += np.sum(mpeak_infall < 300*mp)
            n_false_selection_2 += np.sum(hist["false_selection"])
            for i_cut in range(len(cuts)):
                n_cut_merger[i_cut] += np.sum(
                    hist["mpeak"]/mpeak_merger > cuts[i_cut])
                n_cut_infall[i_cut] += np.sum(
                    hist["mpeak"]/mpeak_infall > cuts[i_cut])

            break
        break

    print()
    print("Counts")
    print("Total: %d, Core Error: %d, False Selection: %d False Selection 2: %d" %
          (n_total, n_core_err, n_false_selection, n_false_selection_2))
    print("Merger Cuts:", n_cut_merger)
    print("Infall Cuts:", n_cut_infall)
    print()
    print("Fractions")
    print("Core Error: %.3f, False Selection: %.3f False Selection 2: %.3f" %
          (n_core_err/n_total, n_false_selection/n_total,
           n_false_selection_2/n_total))
    print("Merger Cuts:", n_cut_merger/n_total)
    print("Infall Cuts:", n_cut_infall/n_total)


if __name__ == "__main__": main()
