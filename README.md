# MD Simulations Analysis 
scripts to visualize Amber molecular dynamics simulations using **VMD** and **Chimera**, analyze rmsf, hbond, secstruct ...etc. using **cpptraj** and Python/R, construct Markov model using **Pyemma**
1. `mdout_nonTCR_process.sh` to combine individual nc files and measure rmsf, hbond, and native contacts between chains
2. `vmd_msm.sh` to visualize the simulation with selected residues and protein within 5 A of the selected residue in licourie, e.g. "185 to 189"
3. `secstruct.sh`, `dihedral.sh`, `ion.sh`, `specialdistance.sh` to calculate secondary structure for selected segments, dihedral angle, distance distribution and ion density. corresponding python scripts to plot them. `dwell_time.py` to do dwell time analysis
4. `commontraj.sh` to strip all but backbone atoms to construct a common trajectory. `project.py` to construct a shared conformational space. 
