
import os

out = "/ssd1/chey120/workspace/chainA_chainA/Y39_mut"

# the chain ID and the residue index of the peptide sequence to be replaced
chain = 'C'
ori = range(1,16)

# input sequences of chainA to replaced the original peptide
chainA = open("Y39_mutated.txt", "r")
sequences = chainA.readlines()
chainA.close()

for seq_w in sequences:
  #strip the whitespace
  seq = seq_w.strip()
  #left to right
  name = str(seq)[0:5] + ".resfile"
  print(name)
  res_file = os.path.join(out, name)
  res = open(res_file, "w") 
  res.write("NATRO \n \nstart\n") #Header: preserve the input rotamer
  for (i, j) in zip(ori, seq): #<PDBNUM> of the peptide
    res.write(str(i) + " " + chain + " PIKAA " + j + "\n") #<CHAIN> of the peptide 
  res.close()

  #right to left
  # name = "reverse_" + name
  # res_file = os.path.join(out, name)
  # res = open(res_file, "w") 
  # res.write("NATRO \n \nstart\n") #Header: preserve the input rotamer
  # for (i, j) in zip(ori, seq[::-1]): #<PDBNUM> of the peptide
  #   res.write(str(i) + " " + chain + " PIKAA " + j + "\n") #<CHAIN> of the peptide 
  # res.close()