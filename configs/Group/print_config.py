# Change to the z=0 Plummer-equivalent force softneing scale of the simulation
# in units of Mpc/h
eps = "0.00036"
# change to the mass of a single particle in units of Msun/h
mp = "2.25585e6"
# Number of particle files per snapshot in the Gadget output
num_snapshot_files = 8

# Change to the list of valid halo names. (You'll probably want to
# auto-generate this)
haloes = ['Halo015', 'Halo024', 'Halo029', 'Halo046', 'Halo055', 'Halo090', 'Halo175', 'Halo183', 'Halo186', 'Halo257', 'Halo286', 'Halo302', 'Halo313', 'Halo331', 'Halo336', 'Halo342', 'Halo346', 'Halo347', 'Halo352', 'Halo383', 'Halo384', 'Halo399', 'Halo428', 'Halo444', 'Halo470', 'Halo496', 'Halo501', 'Halo504', 'Halo551', 'Halo570', 'Halo579', 'Halo581', 'Halo593', 'Halo606', 'Halo656', 'Halo752', 'Halo755', 'Halo759', 'Halo774', 'Halo781', 'Halo785', 'Halo790', 'Halo806', 'Halo834', 'Halo861', 'Halo909', 'Halo927', 'Halo962', 'Halo985']

# Change to the list of z=0 host IDs. (You'll definitely want to auto-generate
# this)
halo_ids = [91159631, 39697241, 93129043, 36695949, 55122063, 85983277, 48324569, 67206804, 34636226, 84293773, 86951809, 80266135, 24205934, 94525382, 42684762, 47710960, 27910352, 31274696, 21284510, 27300026, 30131520, 38777142, 66040602, 86948024, 25185982, 30391531, 16693056, 70754137, 71641444, 22054813, 74813237, 83403126, 26132293, 20850724, 79691691, 20572707, 75778658, 86677143, 29968699, 36236185, 44542205, 52165124, 55275302, 45682125, 86881604, 79465301, 32567263, 47577365, 41260396]

# Change to the snapshots of
halo_snaps = [235]*len(halo_ids)

# Use one string printf verb for the halo name, then two int printf verbs, one
# for the snapshot and the other for the index. Use double %% instead of single
# % for the ints and a single % for the halo name.
snapshot_format = "/sdf/group/kipac/u/ollienad/Group_zoomins/%s/output_new/snapshot_%%03d.%%d"
# Change this to the directory which contains the consistent trees text files.
# Use a string printf verb for the halo name with a single %.
tree_dir = "/sdf/group/kipac/u/ollienad/Group_zoomins/%s/output_new/rockstar/trees"
# Directory where the output data products will go. Use a string printf verb
# for the halo name with a single %.
data_product_dir = "/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns/SymphonyGroup/%s/"

tree_style = "ct_rvmax"

fmt_string = "%%d %%d %s %s %d %s %s %s %s nil" % (
    eps, mp, num_snapshot_files, 
    snapshot_format, tree_dir,
    data_product_dir, tree_style
)

for i in range(len(haloes)):
    h = haloes[i]
    #print(fmt_string)
    print(fmt_string % (halo_ids[i], halo_snaps[i], h, h, h))
