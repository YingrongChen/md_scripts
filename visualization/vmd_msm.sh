#bash /dors/meilerlab/home/chey120/chainA_chainA/scripts/vmd_msm.sh "185 189"

res=$1
cat > add_mols_traj.tcl << EOF
# List all PDB files in the current directory
set filelist [glob *.pdb]

# Iterate over each PDB file
foreach file \$filelist {
    # Extract the filename without the extension
    set filename [file rootname \$file]

    # Load the PDB file into a new molecule
    mol new \$file type pdb waitfor all

    # Get the molecule ID
    set molid [molinfo top]

    mol modstyle 0 \$molid NewCartoon
    mol modcolor 0 \$molid ColorID \$molid

    mol addrep \$molid
    mol modselect 1 \$molid resid $res
    mol modstyle 1 \$molid Licorice
    mol modcolor 1 \$molid ColorID \$molid

    mol smoothrep \$molid 0 5
    mol smoothrep \$molid 1 5
}

EOF

/sb/meilerapps/Linux2/x86_64/bin/vmd -e add_mols_traj.tcl