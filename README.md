# How to Generate Data Products

First, create a directory in `configs/`

``` mkdir configs/MW_mass
```

Next, run the following commands to diagnose any early issues (replacing the `/oak/stanford/orgs/kipac/users/ycwang19/ZEUS/MWmass_new/Halo*` substring with one that accesses all your halo directories)

``` ls /oak/stanford/orgs/kipac/users/ycwang19/ZEUS/MWmass_new/Halo*/output/rockstar/trees/ > configs/MW_mass/tree_locations.txt
du -h /oak/stanford/orgs/kipac/users/ycwang19/ZEUS/MWmass_new/Halo*/output/rockstar/trees/ > configs/MW_mass/tree_sizes.txt
head -n 1 /oak/stanford/orgs/kipac/users/ycwang19/ZEUS/MWmass_new/Halo*/output/rockstar/trees/tree_0_0_0.dat > MW_mass/tree_headers.txt
```

This will create three files that you can browse to diagnose problems. Look through `tree_locations.txt` and see if there are any simulations which don't have tree files yet: these would be jobs that haven't finished running. Look through `tree_sizes.txt` to get a sense for which simulations are going to be time intensive, and check that none of the tree files are suspisciously small. Look through `tree_headers.txt` and confirm that all your simulations used the same version of Rockstar/consistent-trees. 

The code currently asusmes that your header looks like this:

> #scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_Mvir(9) Mvir(10) Rvir(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29) Orig_halo_ID(30) Snap_idx(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Last_mainleaf_depthfirst_ID(34) Tidal_Force(35) Tidal_ID(36) Rs_Klypin Mvir_all M200b M200c M500c M2500c Xoff Voff Spin_Bullock b_to_a c_to_a A[x] A[y] A[z] b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) T/|U| M_pe_Behroozi M_pe_Diemer Halfmass_Radius RVmax

If you have any headers that don't look like this, contact me (Phil) and I'll change the code accordingly.

Next, copy the Python script print_config.py into your config folder

``` cp print_config.py configs/MW_mass
```

Open the file and edit all the lines which have a comment in front of them to match your simulation specifics. You'll need to be familiar with `"%s"`-style printf formatting, but might be able to pattern-match off the examples if you aren't. You'll also need to manually list the halo names: I'd recommend having your `tree_locations.txt` file open in another tab list while you do this. Next, run the Python script and pipe the output to the main config file.

``` python configs/MW_mass/print_config.py > configs/MW_mass/config.txt ```

Now you're ready to start running the pipeline. You run a series of commands as

``` go <file_name>.go path/to/config.txt ```

You'll probably want to do test runs in an interactive session to make sure
everything is working okay. If so, you cna specify the index of the halo you
want to look at as

``` go <file_name>.go path/to/config.txt <index> ```

I'd recommend using the halo with the smallest tree files.

The steps are the following:
- `write_binary_tree.go` - Converts text consistent-trees files into binary
  depth-first files.
- `write_tree_header.go` - Identifies all main branches and annotates them with
  tons of useful information.
- `write_subhalo_files.go` - Collects major subhalos of the central halo
  into a single, easy to access file.
- `tag_particles.go` - Assigns particles to each subhalo and tags them with 
  their infall times and other useful information
- `xv.go` - Extracts x and v for all major subhalos.
- `phi.go` - Computes potentials for all major subhalos.

The first three are responsible for building up files related to merger trees
and the last three are reponsible for tracking paritcles over time. The tree
scripts don't require any particle access and are vastly cheaper to run, so
you may be interested in running them on their own. `phi.go` is by far the
most expensive and most niche, so you may not want to run it even if you're
doing particle tracking.

For the particle tracking files, you'll want to check what level your
highest-reoslution particles are in your simulations (Symphony has them in
Gadget's 1-indexed level, so that's the default).0 If you don't still have the 
original IC configuraiton files lying around, you can check by opening up the
gadget files and seeing which level is the first to have any particles
(using, say, github.com/phil-mansfield/read_gadget.) Once you know that, change
`HRLevel` at the top of those files, if needed.

You may also want to change the number of files that particles are split across.
By default, this number is 8. But you can increase it if this would lead to 
particle files that are too large to comfortably hold in RAM at one time.

