#!/bin/bash
pdb=(`cat list_of_fastas.txt`) 					### load into array
for i in ${pdb[0]}; do						### create variable i from full array pdb
	fasta1=`grep -A1 "^${i}" rosetta.aln | tail -1`	### reads through file, takes line containing pdb id and line after which is the sequence, and then tail takes just sequence
	for j in ${pdb[@]}; do					### creates second variable with same list of pdbs because I want to make a matrix
		fasta2=`grep -A1 "^${j}" rosetta.aln | tail -1`	### same as fasta1 line but with variable j
		echo "## $i ${j}_out" > ${j}.aln			### creates header for alignment file
		echo "#" >> ${j}.aln				### necessary text
		echo "scores_from_program: 0" >> ${j}.aln	### necessary text
		echo "0 $fasta1" >> ${j}.aln			### places target sequence in first line of alignment file
		echo "0 $fasta2" >> ${j}.aln			### places template sequence in second line of alignment file
	done
done

