import os
directory='/dors/meilerlab/home/chey120/mhcii_asyn/MD/hip/5ni9'
for file in os.listdir(directory):
    if file.endswith(".pdb"):
        fin = open(file)
        fout = open(str(file) + "_hip", "wt")
        file_content = fin.readlines()
        # Go over each line individually
        last_chain_id = "A" 
        for line in file_content:
            if (len(line) > 20):
                current_chain_id = str([str(x) for x in line.strip()][21])
                if current_chain_id != last_chain_id:
                    fout.write("TER\n")
                last_chain_id = current_chain_id
            fout.write(line)
        fout.write("TER\n")
        fin.close()
        fout.close()


# fout.write(f"{l[0]}    {l[1]}  {l[2]}  {l[3]} {l[4]}   {l[5]}      {l[6]}   {l[7]}  {l[8]}  {l[9]}  {l[10]}           {l[11]}\n") 
# for file in os.listdir(directory):
#     if file.endswith(".pdb"):
#         print(file)
#         with open(file, 'r') as f:
#             for line in f:
#                 l = line.split()
#                 if (len(l) > 9 and l[3] == "GLN"):
#                     print(l[3])
        # print(file)
        # f = open(file)
        # for line in f:
        #     l = line.split()
        #     print(l[4])
