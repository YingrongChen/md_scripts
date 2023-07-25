# import required module
import os
import pandas as pd
import sys

# assign directory
directory = sys.argv[1]
 
# iterate over files in
# that directory
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    # checking if it is a file
    if os.path.isfile(f):
        df = pd.read_csv(f, sep="\t", header=None)