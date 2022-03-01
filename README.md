# How to Generate Data Products

First, create a directory in `configs/`

``` mkdir configs/MW_mass```

Next, run the following commands to diagnose any early issues (replacing the `/oak/stanford/orgs/kipac/users/ycwang19/ZEUS/MWmass_new/Halo*` substring with one that accesses all your halo directories)

``` ls /oak/stanford/orgs/kipac/users/ycwang19/ZEUS/MWmass_new/Halo*/output/rockstar/trees/ > configs/MW_mass/tree_locations.txt
du -h /oak/stanford/orgs/kipac/users/ycwang19/ZEUS/MWmass_new/Halo*/output/rockstar/trees/tree_*.dat > configs/MW_mass/tree_sizes.txt
head -n 1 /oak/stanford/orgs/kipac/users/ycwang19/ZEUS/MWmass_new/Halo*/output/rockstar/trees/tree_0_0_0.dat > MW_mass/tree_headers.txt ```

This will create three files that you can browse to diagnose problems. Look through `tree_locations.txt` and see if there are any simulations which don't have tree files yet: these would be jobs that haven't finished running. Look through `tree_sizes.txt` to get a sense for which simulations are going to be time intensive, and check that none of the tree files are suspisciously small. Look through `tree_headers.txt` and confirm that all your simulations used the same version of Rockstar/consistent-trees. 

The code currently asusmes that your header looks like this:
``` #scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_Mvir(9) Mvir(10) Rvir(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29) Orig_halo_ID(30) Snap_idx(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Last_mainleaf_depthfirst_ID(34) Tidal_Force(35) Tidal_ID(36) Rs_Klypin Mvir_all M200b M200c M500c M2500c Xoff Voff Spin_Bullock b_to_a c_to_a A[x] A[y] A[z] b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) T/|U| M_pe_Behroozi M_pe_Diemer Halfmass_Radius RVmax ```
If you have any headers that don't look like this, contact me (Phil) and I'll change the code accordingly.

Next, copy the Python script print_config.py into your config folder

``` cp print_config.py configs/MW_mass ```

Open the file and edit all the lines which have a comment in front of them to match your simulation specifics. You'll need to be familiar with `"%s"`-style printf formatting, but might be able to pattern-match off the examples if you aren't. You'll also need to manually list the halo names: I'd recommend having your `tree_locations.txt` file open in another tab list while you do this. Next, run the Python script and pipe the output to the main config file.

``` python configs/MW_mass/print_config.py > configs/MW_mass/config.txt ```