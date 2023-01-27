import os.path as path

# Change to the z=0 Plummer-equivalent force softneing scale of the simulation
# in units of Mpc/h
eps = "0.0012"
# change to the mass of a single particle in units of Msun/h
mp = "1.3e8"
# Number of particle files per snapshot in the Gadget output
num_snapshot_files = 16

# Change to the list of valid halo names. (You'll probably want to
# auto-generate this)
haloes = ['Halo_000', 'Halo_001', 'Halo_002', 'Halo_003', 'Halo_004', 'Halo_005', 'Halo_006', 'Halo_008', 'Halo_009', 'Halo_010', 'Halo_012', 'Halo_013', 'Halo_014', 'Halo_015', 'Halo_016', 'Halo_017', 'Halo_019', 'Halo_020', 'Halo_021', 'Halo_022', 'Halo_023', 'Halo_024', 'Halo_025', 'Halo_026', 'Halo_027', 'Halo_028', 'Halo_029', 'Halo_042', 'Halo_043', 'Halo_044', 'Halo_046', 'Halo_047', 'Halo_050']

# Change to the list of z=0 host IDs. (You'll definitely want to auto-generate
# this)
halo_ids = [10914911, 8770931, 7726254, 5038120, 8527005, 7638512, 5833438, 3892256, 8702277, 7560605, 5962510, 8693489, 4326456, 5419295, 5368875, 5975171, 4476459, 8616866, 5095304, 5879686, 5289770, 5147994, 7533094, 6091609, 17088745, 5930980, 5644044, 3303297, 6552174, 11795977, 3309491, 4234624, 8396547]

# Change to the snapshots of
halo_snaps = [199]*len(halo_ids)

# Use one string printf verb for the halo name, then two int printf verbs, one
# for the snapshot and the other for the index. Use double %% instead of single
# % for the ints and a single % for the halo name.
snapshot_format = "/sdf/group/kipac/u/ollienad/L-Cluster_zoomins/%s/cdm_output/snapshot_%%03d.%%d"
snapshot_format_2 = "/sdf/group/kipac/u/ollienad/L-Cluster_zoomins/%s/cdm_output/snaps/snapshot_%%03d.%%d"
# Change this to the directory which contains the consistent trees text files.
# Use a string printf verb for the halo name with a single %.
tree_dir = "/sdf/group/kipac/u/ollienad/L-Cluster_zoomins/%s/cdm_output/rockstar_rvmax/trees"
# Directory where the output data products will go. Use a string printf verb
# for the halo name with a single %.
data_product_dir = "/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns/SymphonyLCluster/%s/"

tree_style = "ct_rvmax"

fmt_string = "%%d %%d %s %s %d %s %s %s %s nil" % (
    eps, mp, num_snapshot_files, 
    snapshot_format, tree_dir,
    data_product_dir, tree_style
)

fmt_string_2 = "%%d %%d %s %s %d %s %s %s %s nil" % (
    eps, mp, num_snapshot_files, 
    snapshot_format_2, tree_dir,
    data_product_dir, tree_style
)

for i in range(len(haloes)):
    h = haloes[i]
    #print(fmt_string)
    if path.exists((snapshot_format % h) % (0, 0)):
        print(fmt_string % (halo_ids[i], halo_snaps[i], h, h, h))
    else:
        print(fmt_string_2 % (halo_ids[i], halo_snaps[i], h, h, h))
