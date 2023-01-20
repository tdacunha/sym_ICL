# How to Generate Data Products

## Installing software

1. Have this repo cloned somewhere in your local directory
2. Install the numpy, scipy, and matplotlib
3. Install symlib, using `pip install symlib -U`
4. Install my plotting library, palette. Clone https://github.com/phil-mansfield/palette into a local directory and change your .bashrc file in your home directory so it has the following line in it: `export PYTHONPATH=$PYTHONPATH:/absolute/path/to/palette` (use `pwd` in the directory to get the absoluute path)
5. Install the Go compiler. Ins ome local directory, run `wget https://go.dev/dl/go1.19.5.linux-amd64.tar.gz`. Go to this webpage: https://go.dev/doc/install, switch the tab under "Go installation" to "Linux", then follow the instructions. Replace `/usr/local` with the absolute path of whatever directory you want your go compiler to be in.
6. Install my tree-code library, gravitree. Clone https://github.com/phil-mansfield/gravitree into some local directory, change the .bashrc file to have the line `export PYTHONPATH=$PYTHONPATH:/absolute/path/to/gravitree/python` in it. Go into `gravitree/python` and run `./build_script.sh`.


## Running the pipeline

First, create a directory in `configs/`

``` 
mkdir configs/MW_mass
```

Next, run the following commands to diagnose any early issues (replacing the `/oak/stanford/orgs/kipac/users/ycwang19/ZEUS/MWmass_new/Halo*` substring with one that accesses all your halo directories)

``` 
ls /oak/stanford/orgs/kipac/users/ycwang19/ZEUS/MWmass_new/Halo*/output/rockstar/trees/ > configs/MW_mass/tree_locations.txt
du -h /oak/stanford/orgs/kipac/users/ycwang19/ZEUS/MWmass_new/Halo*/output/rockstar/trees/ > configs/MW_mass/tree_sizes.txt
head -n 1 /oak/stanford/orgs/kipac/users/ycwang19/ZEUS/MWmass_new/Halo*/output/rockstar/trees/tree_0_0_0.dat > MW_mass/tree_headers.txt
```

Note that if you have very complicated, non uniform directory structure, you can
pipe arguments out of config files and around with awk and xargs like this:

```
cat config.txt | awk '{print $7}' | xargs ls > tree_locations.txt
```

This will create three files that you can browse to diagnose problems. Look through `tree_locations.txt` and see if there are any simulations which don't have tree files yet: these would be jobs that haven't finished running. Look through `tree_sizes.txt` to get a sense for which simulations are going to be time intensive, and check that none of the tree files are suspisciously small. Look through `tree_headers.txt` and confirm that all your simulations used the same version of Rockstar/consistent-trees. 

The code currently asusmes that your header looks like this:

> #scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_Mvir(9) Mvir(10) Rvir(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29) Orig_halo_ID(30) Snap_idx(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Last_mainleaf_depthfirst_ID(34) Tidal_Force(35) Tidal_ID(36) Rs_Klypin Mvir_all M200b M200c M500c M2500c Xoff Voff Spin_Bullock b_to_a c_to_a A[x] A[y] A[z] b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) T/|U| M_pe_Behroozi M_pe_Diemer Halfmass_Radius RVmax

If you have any headers that don't look like this, contact me (Phil) and I'll change the code accordingly.

Next, copy the Python script print_config.py into your config folder

``` 
cp print_config.py configs/MW_mass
```

Open the file and edit all the lines which have a comment in front of them to match your simulation specifics. You'll need to be familiar with `"%s"`-style printf formatting, but might be able to pattern-match off the examples if you aren't. You'll also need to manually list the halo names: I'd recommend having your `tree_locations.txt` file open in another tab list while you do this. Next, run the Python script and pipe the output to the main config file.

```
python configs/MW_mass/print_config.py > configs/MW_mass/config.txt
```

Now you're ready to start running the pipeline. You run a series of commands as

```
go run <file_name>.go path/to/config.txt <index>
```

Here, if <index> is -1, it means that the file runs on everything in the suite, and if it's a non-negative number, only the halo at that index will be run.

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

The first three are responsible for building up files related to merger trees
and the last three are reponsible for tracking paritcles over time. The tree
scripts don't require any particle access and are vastly cheaper to run, so
you may be interested in running them on their own.

For the particle tracking files, you'll want to check what level your
highest-reoslution particles are in your simulations (Symphony has them in
Gadget's 1-indexed level, so that's the default). If you don't still have the 
original IC configuraiton files lying around, you can check by opening up the
gadget files and seeing which level is the first to have any particles
(using, say, github.com/phil-mansfield/read_gadget.) Once you know that, change
`HRLevel` at the top of those files, if needed.

You may also want to change the number of files that particles are split across.
By default, this number is 8. But you can increase it if this would lead to 
particle files that are too large to comfortably hold in RAM at one time.

The `job_tree.sh` batch file will create tree files for one suite, and the `job_track.sh` file will parallelize particle tracking across the suite. You'll need to change the config variable so that it points to the right file, and you'll need to change the `array` variable in the header so that it corresponds to the number of halos in the suite. You'll also want to give the job a new name, and you'll want to make sure the log/error files are being written to a real directory.

Lastly, you run into the absolutely most annoying part of the pipeline, which is running the halo finder. (If you, dear reader, have any interest in making this less miserable, let me know.) There are three python scripts that need to be run. The first is find_infall_cores.py, which finds the "core" particles for each halo. The second is `print_core_catalogue.py`, which is the acutal halo finder and the thrid is `convert_core_catalogue.py`, which turns the plain-text output of the first file into a binary format. The script `job_track.sh` runs and parallelizes both of these.

The file itself looks something like this
```
config=configs/MilkyWay/config.txt
suffix=k64_n30
snap_range=235:235

python3 find_infall_cores.py ${config} ${SLURM_ARRAY_TASK_ID} &&              
	python3 print_core_catalogue.py ${config} ${SLURM_ARRAY_TASK_ID} --reset --suffix=${suffix} --snap_range=${snap_range} &&
   python3 convert_core_catalogue.py ${config} ${SLURM_ARRAY_TASK_ID} --suffix=${suffix} &&
   echo "done"
```

The `config` variable should be set to the suite's config file, `suffix` allows you to have multiple runs of the subhalo finders with different parameters, so if you don't care about that, delete this line and remove the `--suffix=...` bit from the script. snap_range lets you choose which snapshot range you run the (i.e. "first":"last"). The pipeline will know if your job crashe dhalfway through and won't recalculate snapshots a second time, so only mess with this if you're trying ot parallelize or only look at late times or something.

print_core_catalogue.py has an additional flag, --reset. When not set, the pipeline looks for half finished files, figures out where you crashed, and starts from there. If you're restarting because something was wrong and needed to be fixed, you can cause it to start over with the --reset flag.

Small simulations (e.g. MilkyWay, LMC) usually finish in about a day, while bigger sims (like the HR runs of the Cluster suite) take several days. This is because I never actually parallelized things and it's all single core other than what you can do through array jobs. What this means is that it's likely that the maximum job size SDF allows is shorter than the length of the job. To deal with this, check the log files of all the jobs that finish whnever you check in and see if they finished due to an error or due to completing the job. Rerun the job script on the ones that crashed due to timeout error. You can specify which ones to rerun by changing the SBATCH array variable to array=1,5,20 or whatever.

Making data downloadable
------------------------

If you haven't already, you'll want to add your suite to list of objects that symlib can analyze. There are enough special cases that this has to be done manually, I'm afraid.
- Add an entry to symlib.paramter_table corresponding to your suite. (symlib/lib.py)
- Add your simulation to symlib.SUITE_NAMES (symlib/util.py)
- Add an entry symlib.DEAFULT_HALO_NAMES (symlib/util.py)
- Generate scale_factors.py. To do this, call go run scale_factor_table <suite name> <config file>. This reads in the scale factors of every snapshot of every halo, so you might not want to do this on a login node. Once you've done this for your suite, add a line to symlib/download_tables/generate_scale_factor_table.py
s main function corresponding to your suite.

(TODO: add packing and auto-generating instructions)