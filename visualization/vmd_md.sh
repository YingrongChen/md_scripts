res=$1

cat > add_mols_traj.tcl << EOF
# Clear the current molecules in VMD
mol delete all

# Get a list of all parm7 files in the current directory
set parm7_files [glob -nocomplain *.parm7]

# Loop through each parm7 file
foreach parm7_file \$parm7_files {
  # Extract the filename without the extension
  set basename [file rootname [file tail \$parm7_file]]

  # Load the parm7 file
  mol new \$parm7_file type parm7 waitfor all

  # Get the list of corresponding NetCDF trajectory files
  set nc_files [glob -nocomplain "\${basename}_trial*.offset_combine.nc"]

  # Loop through each NetCDF trajectory file
  foreach nc_file \$nc_files {
    # Load the NetCDF trajectory file
    mol addfile \$nc_file waitfor all first 0 last -1 step 1 filebonds 1 autobonds 1
  }

  # Get the molecule ID
  set molid [molinfo top]
  mol addrep \$molid
  mol modselect 1 \$molid same resid as (within 3 of resid $res)
  mol modstyle 1 \$molid Licorice 0.100000 10.000000 10.000000
  mol modcolor 1 \$molid Name

  mol addrep \$molid
  mol modselect 2 \$molid resid $res
  mol modstyle 2 \$molid Licorice 0.100000 10.000000 10.000000
  mol modcolor 2 \$molid Type

  mol smoothrep \$molid 0 5
  mol smoothrep \$molid 1 5
  mol smoothrep \$molid 2 5
}

# Center and zoom the view
animate goto 0

EOF

/sb/meilerapps/Linux2/x86_64/bin/vmd -e add_mols_traj.tcl
