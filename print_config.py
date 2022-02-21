# Change to the z=0 Plummer-equivalent force softneing scale of the simulation
# in units of Mpc/h
eps = "0.00017"
# change to the mass of a single particle in units of Msun/h
mp = "2.81981e5"
# Number of particle files per snapshot in the Gadget output
num_snapshot_files = 8
# Change to the list of valid halo IDs/names.
haloes = ["023", "008", "119", "188", "247", "268", "270", "288", "327", "349",
          "364", "374", "414", "415", "416", "440", "460", "469", "490", "530",
          "558", "567", "570", "606", "628", "641", "675", "718", "738", "749",
          "797", "800", "825", "829", "852", "878", "881", "925", "926", "937",
          "939", "967", "9749", "9829", "990"]

# Use one string printf verb for the halo name, then two int printf verbs, one
# for the snapshot and the other for the index. Use double %% instead of single
# % for the ints and a single % for the halo name.
snapshot_format = "/oak/stanford/orgs/kipac/users/ycwang19/ZEUS/MWmass_new/Halo%s/output/snapshot_%%03d.%%d"
# Change this to the directory which contains the consistent trees text files.
# Use a string printf verb for the halo name with a single %.
tree_dir = "/oak/stanford/orgs/kipac/users/ycwang19/ZEUS/MWmass_new/Halo%s/output/rockstar/trees/"
# Directory where the output data products will go. Us a string printf verb
# for the halo name with a single %.
data_product_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/Symphony/MW_mass/Halo%s/"

fmt_string = "%s %s %d %s %s %s" % (eps, mp, num_snapshot_files, 
                                    snapshot_format, tree_dir,
                                    data_product_dir)

for i in range(len(haloes)):
    h = haloes[i]
    print(fmt_string % (h, h, h))
