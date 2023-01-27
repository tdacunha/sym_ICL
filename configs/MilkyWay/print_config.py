# Change to the z=0 Plummer-equivalent force softneing scale of the simulation
# in units of Mpc/h
eps = "0.00017"
# change to the mass of a single particle in units of Msun/h
mp = "2.81981e5"
# Number of particle files per snapshot in the Gadget output
num_snapshot_files = 8

# Change to the list of valid halo names.
haloes = ['Halo023', 'Halo088', 'Halo119', 'Halo188', 'Halo247', 'Halo268', 'Halo270', 'Halo288', 'Halo327', 'Halo349', 'Halo364', 'Halo374', 'Halo414', 'Halo415', 'Halo416', 'Halo440', 'Halo460', 'Halo469', 'Halo490', 'Halo530', 'Halo558', 'Halo567', 'Halo570', 'Halo606', 'Halo628', 'Halo641', 'Halo675', 'Halo718', 'Halo738', 'Halo749', 'Halo797', 'Halo800', 'Halo825', 'Halo829', 'Halo852', 'Halo878', 'Halo881', 'Halo925', 'Halo926', 'Halo937', 'Halo939', 'Halo967', 'Halo9749', 'Halo9829', 'Halo990']

# Change to the list of z=0 host IDs. (You'll probably want to auto-generate
# this)
halo_ids = [7019390, 31120521, 10208174, 28839883, 6646440, 7287306, 31107790, 8697419, 8391099, 15119051, 10301677, 9405794, 9487756, 29718260, 8297694, 14783515, 8932799, 30280719, 9659071, 7714515, 9967184, 8104130, 19722077, 27371347, 12638322, 30457872, 23284353, 28833029, 28077485, 12607178, 9948707, 9113976, 9721967, 6414883, 7961010, 10422676, 42248692, 8529408, 8839742, 9785057, 8282747, 18133566, 27982424, 18701512, 11431405]

# Change to the snapshots of
halo_snaps = [235]*len(halo_ids)

# Use one string printf verb for the halo name, then two int printf verbs, one
# for the snapshot and the other for the index. Use double %% instead of single
# % for the ints and a single % for the halo name.
snapshot_format = "/sdf/group/kipac/u/ollienad/MW_zoomins/%s/output/snapshot_%%03d.%%d"
# Change this to the directory which contains the consistent trees text files.
# Use a string printf verb for the halo name with a single %.
tree_dir = "/sdf/group/kipac/u/ollienad/MW_zoomins/%s/output/rockstar/trees/"
# Directory where the output data products will go. Use a string printf verb
# for the halo name with a single %.
data_product_dir = "/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns/SymphonyMilkyWay/%s/"


um_fmt = "/sdf/group/kipac/u/ycwang/MWmass_new/%s/output/rockstar/groupcat/sfr_catalog_%%.6f.txt"

tree_style = "ct_rvmax"
fmt_string = ("%%d %%d %s %s %d %s %s %s %s %s" %
              (eps, mp, num_snapshot_files, 
               snapshot_format, tree_dir,
               data_product_dir, tree_style, um_fmt))

for i in range(len(haloes)):
    h = haloes[i]
    print(fmt_string % (halo_ids[i], halo_snaps[i], h, h, h, h))
