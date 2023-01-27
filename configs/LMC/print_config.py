# Change to the z=0 Plummer-equivalent force softneing scale of the simulation
# in units of Mpc/h
eps = "0.000080"
# change to the mass of a single particle in units of Msun/h
mp = "3.52476e4"
# Number of particle files per snapshot in the Gadget output. 
num_snapshot_files = 16

# Change to the list of valid halo names. (You'll probably want to
# auto-generate this)
haloes = ['Halo032', 'Halo059', 'Halo0662', 'Halo083', 'Halo088', 'Halo097', 'Halo104', 'Halo110', 'Halo202', 'Halo208', 'Halo218', 'Halo296', 'Halo301', 'Halo303', 'Halo340', 'Halo374', 'Halo380', 'Halo391', 'Halo405', 'Halo440', 'Halo463', 'Halo4662', 'Halo511', 'Halo524', 'Halo539', 'Halo567', 'Halo575', 'Halo602', 'Halo697', 'Halo711', 'Halo721', 'Halo767', 'Halo802', 'Halo824', 'Halo850', 'Halo853', 'Halo914', 'Halo932', 'Halo933']

# Change to the list of z=0 host IDs. (You'll definitely want to auto-generate
# this)
halo_ids = [6809161, 11451612, 4714720, 2614218, 5470740, 2572433, 2104392, 8920770, 9048851, 5571100, 4797341, 9373983, 3585923, 5676033, 2705971, 2781538, 6323270, 4919362, 4865033, 9443609, 4414397, 5613596, 4210297, 4382204, 3371981, 3775946, 4263000, 9045701, 4524194, 4286404, 4959649, 3145683, 2965025, 6225728, 6473303, 5309166, 7168215, 2680560, 2599455]

# Change to the snapshots of
halo_snaps = [235]*len(halo_ids)

# Use one string printf verb for the halo name, then two int printf verbs, one
# for the snapshot and the other for the index. Use double %% instead of single
# % for the ints and a single % for the halo name.
snapshot_format = "/sdf/group/kipac/u/ollienad/LMC_zoomins/%s/output/snapshot_%%03d.%%d"
# Change this to the directory which contains the consistent trees text files.
# Use a string printf verb for the halo name with a single %.
tree_dir = "/sdf/group/kipac/u/ollienad/LMC_zoomins/%s/output/rockstar/trees"
# Directory where the output data products will go. Use a string printf verb
# for the halo name with a single %.
data_product_dir = "/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns/SymphonyLMC/%s/"

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
