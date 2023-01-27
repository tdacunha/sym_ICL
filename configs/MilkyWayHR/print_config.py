# Change to the z=0 Plummer-equivalent force softneing scale of the simulation
# in units of Mpc/h
eps = "%f" % (0.00017 / 2)
# change to the mass of a single particle in units of Msun/h
mp = "%g" % (2.81981e5 / 8)
# Number of particle files per snapshot in the Gadget output
num_snapshot_files = 16

# Change to the list of valid halo names.
haloes = ["Halo023", "Halo247", "Halo268", "Halo530", "Halo829"]

# Change to the list of z=0 host IDs. (You'll probably want to auto-generate
# this)
halo_ids = [54341170, 51673895, 57215907, 60602989, 50139819]

# Change to the snapshots of
halo_snaps = [235]*len(halo_ids)

# Use one string printf verb for the halo name, then two int printf verbs, one
# for the snapshot and the other for the index. Use double %% instead of single
# % for the ints and a single % for the halo name.
snapshot_format = "/sdf/group/kipac/u/ollienad/MW_zoomins/%s/output_16K/snapshot_%%03d.%%d"
# Change this to the directory which contains the consistent trees text files.
# Use a string printf verb for the halo name with a single %.
tree_dir = "/sdf/group/kipac/u/ollienad/MW_zoomins/%s/output_16K/rockstar/trees/"
# Directory where the output data products will go. Use a string printf verb
# for the halo name with a single %.
data_product_dir = "/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns/SymphonyMilkyWayHR/%s/"

fmt_string = "%%d %%d %s %s %d %s %s %s ct_rvmax nil" % (
    eps, mp, num_snapshot_files, 
    snapshot_format, tree_dir,
    data_product_dir
)

for i in range(len(haloes)):
    h = haloes[i]
    print(fmt_string % (halo_ids[i], halo_snaps[i], h, h, h))
