import glob
import os
import chimera
import sys
#important residues to view as licorice
# res = int(sys.argv[1])

# Get a list of all parm7 files in the current directory
parm7_files = glob.glob("*.parm7")

# Loop through each parm7 file
for parm7_file in parm7_files:
    # Extract the filename without the extension
    basename = os.path.splitext(os.path.basename(parm7_file))[0]

    # Load the parm7 file
    chimera.openModels.open(parm7_file, type="AMBER7", merge="skip")

    # Get the list of corresponding NetCDF trajectory files
    nc_files = glob.glob("{}_trial*.offset_combine.nc".format(basename))

    # Loop through each NetCDF trajectory file
    for nc_file in nc_files:
        # Load the NetCDF trajectory file
        chimera.openModels.open(nc_file, type="MD", first=0, last=-1, step=1, modelID=1)

    # Get the molecule ID
    # mol = chimera.openModels.list()[0]
    # molID = mol.id

    # # Create two new representations for the molecule
    # chimera.runCommand("represent licorice #1 sel :{}-{} style 0.1 radius 10 color name".format(res-3, res+3))  
    # chimera.runCommand("represent licorice #2 sel :{} style 0.1 radius 10 color type".format(res))

    # # Smooth the representations
    # chimera.runCommand("smooth #1 5")
    # chimera.runCommand("smooth #2 5")

# Center and zoom the view
# chimera.runCommand("view center")
# chimera.runCommand("view fit")