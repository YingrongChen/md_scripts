#! /usr/bin/env python
from argparse import ArgumentParser

####################################
# Simple PDB clean script. Written primarily to easily separate multiple occupancies
# for residues and ligands in PDB files
#
# Written by: Benjamin P. Brown
# Date: July 04, 2020
#####################################

# Main
def main():
    # Load commandline options
    args = cmdlineparse()

    # Read in record types
    type_list = [str(r_type) for r_type in args.record_types]

    # Read in states
    state_list = [str(state) for state in args.states]

    # Prepare output file
    outfile = open(args.output_pdb, 'wb')

    # Read in the PDB file
    new_pdb = []
    with open(args.input_pdb, 'r') as infile:
        file_content = infile.readlines()
        # Go over each line individually
        last_chain_id = "QQQQQ" # Dummy value
        for line in file_content:
            # Only consider lines beginning with allowed record types (e.g. ATOM, HETATM, CONECT, etc.)
            for r_type in type_list:
                record_type = str([str(x) for x in line.strip().split(' ')][0])
                if r_type == record_type:
                    # Collect chain ID
                    current_chain_id = str([str(x) for x in line.strip()][21])
                    if last_chain_id == "QQQQQ":
                        last_chain_id = current_chain_id
                    if len(args.skip_chains):
                        bad_chain = False
                        for chain in args.skip_chains:
                            if current_chain_id == chain:
                                bad_chain = True
                        if bad_chain:
                            break;
                    # Insert TER record after chainbreak
                    if args.insert_ter_records and current_chain_id != last_chain_id:
                        new_pdb.append("TER\n")
                        last_chain_id = current_chain_id
                    # Only output lines of the desired state
                    if len(state_list):
                        for state in state_list:
                            state_type = str([str(x) for x in line.strip()][16])
                            if state_type == state:
                                new_pdb.append(line.strip() + "\n")
                        # No point going through the other record_types if this one was the one
                        break
    # Add END record
    if args.insert_end_record:
        new_pdb.append("END")

    # Rename chain IDs
    if args.rename_chain:
        chain_pairs = ParseChains( args.rename_chain)
        for pair in chain_pairs:
            new_pdb = RenameChain(new_pdb, pair[0], pair[1])

    # Rename state IDs
    if args.rename_state:
        state_pairs = args.rename_state
        for pair in state_pairs:
            new_pdb = RenameState(new_pdb, pair[0], pair[1])

    # Rename segment IDs
    if args.rename_segment:
        segment_pairs = args.rename_segment
        for pair in segment_pairs:
            new_pdb = RenameSegment(new_pdb, pair[0], pair[1])

    # Rename HETATM chain IDs
    if args.rename_hetatm_chain:
        new_chain_id = str(args.rename_hetatm_chain)
        new_pdb = RenameHetatmChain(new_pdb,  new_chain_id)

    # Rename HETATM residue ID
    if args.rename_hetatm_resid:
        new_resid = str(args.rename_hetatm_resid)
        new_pdb = RenameHetatmResID(new_pdb,  new_resid)

    # Renumber the final PDB
    if args.renumber:
        new_pdb = Renumber( new_pdb, int(args.renumber_start))

    # Write final new file
    for line in new_pdb:
        outfile.write(line)

    # Close the output file and end
    outfile.close()
    return 0

def RenameChain(NEW_PDB, OLD_CHAIN_ID, NEW_CHAIN_ID):
    new_pdb = []
    for line in NEW_PDB:
        record_type = str([str(x) for x in line.strip().split(' ')][0])
        if record_type == "ATOM" or record_type == "HETATM" or record_type == "ANISOU":
            if line[21] == OLD_CHAIN_ID:
                line = line[0:21] + NEW_CHAIN_ID + line[22:-1]
        new_pdb.append(line.strip() + "\n")
    return new_pdb

def RenameHetatmChain(NEW_PDB, NEW_CHAIN_ID):
    new_pdb = []
    for line in NEW_PDB:
        record_type = str([str(x) for x in line.strip().split(' ')][0])
        if record_type == "HETATM":
            line = line[0:21] + NEW_CHAIN_ID + line[22:-1]
        new_pdb.append(line.strip() + "\n")
    return new_pdb

def RenameHetatmResID(NEW_PDB, NEW_RESID):
    new_pdb = []
    for line in NEW_PDB:
        record_type = str([str(x) for x in line.strip().split(' ')][0])
        if record_type == "HETATM":
            line = line[0:17] + NEW_RESID + line[20:-1]
        new_pdb.append(line.strip() + "\n")
    return new_pdb

def RenameState(NEW_PDB, OLD_STATE_ID, NEW_STATE_ID):
    new_pdb = []
    for line in NEW_PDB:
        record_type = str([str(x) for x in line.strip().split(' ')][0])
        if record_type == "ATOM" or record_type == "HETATM" or record_type == "ANISOU":
            if line[16] == OLD_STATE_ID:
                line = line[0:16] + NEW_STATE_ID + line[17:-1]
        new_pdb.append(line.strip() + "\n")
    return new_pdb

def RenameSegment(NEW_PDB, OLD_SEG_ID, NEW_SEG_ID):
    new_pdb = []
    for line in NEW_PDB:
        record_type = str([str(x) for x in line.strip().split(' ')][0])
        if record_type == "ATOM" or record_type == "HETATM" or record_type == "ANISOU":
            if line[72] == OLD_SEG_ID:
                line = line[0:72] + NEW_SEG_ID + line[73:-1]
        new_pdb.append(line.strip() + "\n")
    return new_pdb

def Renumber(NEW_PDB, START=1):
    new_pdb = []
    new_resid = int(START)
    last_resid = int(START - 1)  # Dummy value
    for line in NEW_PDB:
        record_type = str([str(x) for x in line.strip().split(' ')][0])
        if record_type == "ATOM" or record_type == "HETATM" or record_type == "ANISOU":
            last_resid = int(line[24:28])
            break
    for line in NEW_PDB:
        record_type = str([str(x) for x in line.strip().split(' ')][0])
        if record_type == "ATOM" or record_type == "HETATM" or record_type == "ANISOU":
            current_resid = int(line[24:28])
            if current_resid is not last_resid:
                last_resid = current_resid
                new_resid += int( 1)
            # Hacky as shit but I suck at Python so here we are
            if new_resid < int(10):
                line = line[0:24] + "  " + str(new_resid) + line[28:-1]
            elif new_resid < int(100):
                line = line[0:24] + " " + str(new_resid) + line[28:-1]
            elif new_resid < int(1000):
                line = line[0:24] + str(new_resid) + line[28:-1]
        new_pdb.append(line.strip() + "\n")
    return new_pdb

def ParseChains( CHAIN_IDS):
    n_chains = len( CHAIN_IDS) / 2
    if len( CHAIN_IDS) % 2 == 0:
        chains = [[ CHAIN_IDS[col + (row * 2)] for col in range(2) ] for row in range(n_chains)]
    else:
        print("Old and new chain IDs must be given in pairs!")
    return chains


# Argument parser
def cmdlineparse():
    parser = ArgumentParser(description="command line arguments")
    parser.add_argument("-input_pdb", dest="input_pdb", required=True, help="Input PDB file to be cleaned", metavar="<input_pdb>")
    parser.add_argument("-output_pdb", dest="output_pdb", required=True, help="Output filename for the cleaned PDB", metavar="<output_pdb>")
    parser.add_argument("-record_types", dest="record_types", required=False, default="ATOM", nargs="+",
                        help="Specify PDB record types to keep from original file (e.g. ATOM, HETATM, CONECT, etc.)",
                        metavar="<record_type>")
    parser.add_argument("-states", dest="states", required=False, default=" ", nargs="+",
                        help="Specify states to keep as indicated by position 16 (0-indexed) in the PDB file;"
                             "by default, only includes unspecified (empty position 16) states",
                        metavar="<state>")
    parser.add_argument("-skip_chains", dest="skip_chains", required=False, default="", nargs="+",
                        help="Specify chains to skip when writing new PDB file",
                        metavar="<chain_id>")
    parser.add_argument("-insert_ter_records", dest="insert_ter_records", required=False, default=False, action="store_true",
                        help="Insert TER records after each chainbreak")
    parser.add_argument("-insert_end_record", dest="insert_end_record", required=False, default=False, action="store_true",
                        help="Insert END record at end of file")
    parser.add_argument("-rename_chain", dest="rename_chain",  required=False, nargs="+",
                        help="Pass two chain IDs, first the chain ID to be altered and second the new chain ID; "
                             "can be specified multiple times on the command-line; note - occurs before"
                             "'rename_hetatm', and thus passing 'rename_hetatm' potentially overwrites changes specified"
                             "here", metavar="<old_chain_ID new_chain_ID>")
    parser.add_argument("-rename_state", dest="rename_state",  required=False, action="append", nargs=2,
                        help="Pass two state IDs, first the state ID to be altered and second the new state ID; "
                             "can be specified multiple times on the command-line",
                        metavar="<old_state_ID new_state_ID>")
    parser.add_argument("-rename_segment", dest="rename_segment",  required=False, action="append", nargs=2,
                        help="Pass two segment IDs, first the segment ID to be altered and second the new segment ID; "
                             "can be specified multiple times on the command-line",
                        metavar="<old_segment_ID new_segment_ID>")
    parser.add_argument("-rename_hetatm_chain", dest="rename_hetatm_chain", required=False,
                        help="Set HETATM chains to this chain ID; note - occurs after 'rename_chain' and "
                             "therefore takes precedence", metavar="<rename_hetatm>")
    parser.add_argument("-rename_hetatm_resid", dest="rename_hetatm_resid", required=False,
                        help="Set HETATM residue ID to this residue ID", metavar="<rename_hetatm>")
    parser.add_argument("-renumber", dest="renumber", required=False, help="Renumber chain", action="store_true")
    parser.add_argument("-renumber_start", dest="renumber_start", required=False, default="1",
                        help="If renumbering chain, start from this value", metavar="<renumber_start>")
    args=parser.parse_args()
    return args

# Run from commandline
if __name__ == '__main__':
    main()
