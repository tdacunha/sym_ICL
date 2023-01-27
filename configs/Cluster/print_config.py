import os

# Change to the z=0 Plummer-equivalent force softneing scale of the simulation
# in units of Mpc/h
eps = "0.00325"
# change to the mass of a single particle in units of Msun/h
mp = "1.3e8"
# Number of particle files per snapshot in the Gadget output
num_snapshot_files = 4

# Change to the list of valid halo names. (You'll probably want to
# auto-generate this)
haloes = ['Halo156', 'Halo175', 'Halo200', 'Halo211', 'Halo213', 'Halo222', 'Halo225', 'Halo266', 'Halo274', 'Halo277', 'Halo282', 'Halo293', 'Halo304', 'Halo305', 'Halo306', 'Halo308', 'Halo317', 'Halo321', 'Halo324', 'Halo326', 'Halo335', 'Halo337', 'Halo339', 'Halo345', 'Halo346', 'Halo348', 'Halo349', 'Halo352', 'Halo354', 'Halo358', 'Halo360', 'Halo361', 'Halo366', 'Halo367', 'Halo377', 'Halo378', 'Halo385', 'Halo386', 'Halo387', 'Halo390', 'Halo391', 'Halo394', 'Halo400', 'Halo407', 'Halo409', 'Halo415', 'Halo416', 'Halo419', 'Halo428', 'Halo429', 'Halo436', 'Halo437', 'Halo441', 'Halo445', 'Halo447', 'Halo448', 'Halo452', 'Halo454', 'Halo455', 'Halo456', 'Halo461', 'Halo462', 'Halo465', 'Halo471', 'Halo472', 'Halo474', 'Halo475', 'Halo476', 'Halo478', 'Halo479', 'Halo480', 'Halo483', 'Halo489', 'Halo494', 'Halo502', 'Halo517', 'Halo518', 'Halo522', 'Halo529', 'Halo544', 'Halo545', 'Halo546', 'Halo561', 'Halo572', 'Halo574', 'Halo595', 'Halo600', 'Halo604', 'Halo629', 'Halo631', 'Halo639', 'Halo645', 'Halo653', 'Halo734', 'Halo372', 'Halo425']

skipped_snaps = {
    "Halo305": [174],
    "Halo518": [78, 79],
    "Halo529": [91, 92, 93, 94],
    "Halo595": [171],
    "Halo604": [175],
    "Halo631": [181],
    "Halo653": [173],
    "Halo416": [121, 122],
    "Halo419": [127, 128, 129],
    "Halo429": [126, 126, 127],
    "Halo436": [126, 127, 128, 129],
    "Halo452": [186],
    "Halo462": [95, 105],
    "Halo465": [103],
    "Halo471": [126],
}

# Change to the list of z=0 host IDs. (You'll definitely want to auto-generate
# this)
halo_ids = [26306919, 22872262, 22304649, 25823426, 21364727, 22830522, 19444013, 21923711, 17996154, 19509726, 19345993, 22511822, 17053913, 20283686, 23487669, 20709295, 18055416, 26385776, 16843658, 15093842, 12555842, 20281263, 20213371, 19806854, 19810872, 24486951, 22930905, 20113301, 22372037, 19107019, 20832433, 20378485, 18055846, 27268269, 18706415, 18735140, 15504288, 14535742, 19980448, 18078701, 13886213, 20119654, 16498354, 19811073, 17194681, 17811414, 20249110, 16614205, 17768477, 15974521, 14709200, 23877693, 16689965, 16010380, 16796052, 16409038, 23882189, 18682260, 17373356, 15680496, 13966698, 18097409, 16528779, 14774465, 17387585, 16469506, 19266616, 13698986, 15433414, 16822101, 16621939, 12937822, 20592825, 16648506, 14740441, 12110250, 18990859, 13155267, 13726813, 14898766, 15620652, 12736545, 22388339, 15134974, 19606197, 14031499, 14048626, 16280758, 15035153, 16529352, 17222345, 16150861, 13429728, 16271339, 14479982, 19214605]

# Change to the snapshots of
halo_snaps = [199]*len(halo_ids)

for i in range(len(halo_snaps)):
    if haloes[i] in skipped_snaps:
        halo_snaps[i] -= len(skipped_snaps[haloes[i]])

def all_halos(parent_dir):
    return [file for file in os.listdir(parent_dir) if 
            len(file) > 4 and file[:4] == "Halo"]

particle_dir_1 = "/sdf/group/kipac/g/cosmo/ki14/Done"
dir_1_halos = all_halos(particle_dir_1)
particle_dir_2 = "/sdf/group/kipac/g/cosmo/ki14/Rhapsody8192/Done"
dir_2_halos = all_halos(particle_dir_2)
particle_dir_3 = "/sdf/group/kipac/g/cosmo/ki20/Rhapsody_8192_Sample2"
dir_3_halos = all_halos(particle_dir_3)

halo_dir = "/sdf/group/kipac/u/ycwang/rhapsody_halo"

def is_buggy_dir(parent_dir, halo):
    files = os.listdir(parent_dir)
    return halo + "_2" in files

snapshot_format = "/sdf/group/kipac/u/ollienad/Group_zoomins/%s/output_new/snapshot_%%03d.%%d"
tree_dir = "/sdf/group/kipac/u/ollienad/Group_zoomins/%s/output_new/rockstar/trees"
data_product_dir = "/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns/SymphonyCluster/%s/"

tree_style = "ct_rhapsody"

for i in range(len(haloes)):
    h = haloes[i]
    
    if is_buggy_dir(halo_dir, h):
        h_tree = h + "_2"
    else:
        h_tree = h

    tree_dir = "%s/%s/rockstar/trees" % (halo_dir, h_tree)
    if not os.path.exists(tree_dir):
        tree_dir = tree_dir + "_old"
    
    if h in dir_1_halos:
        p_dir = particle_dir_1
    elif h in dir_2_halos:
        p_dir = particle_dir_2
    elif h in dir_3_halos:
        p_dir = particle_dir_3
    else:
        continue

    snapshot_format = "%s/%s/output_HR/snapshot_%%03d.%%d" % (p_dir, h)

    print("%d %d %s %s %d %s %s %s %s nil" % (
        halo_ids[i], halo_snaps[i], eps, mp,
        num_snapshot_files, snapshot_format,
        tree_dir, data_product_dir % h, tree_style
    ))
