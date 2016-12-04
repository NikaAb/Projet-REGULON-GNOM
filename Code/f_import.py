"""
Projet REGULON - GENOM
Nika Abdollahi & Melissa Cardon

12-2016
"""

import numpy as np


#========================================================================
#                          Import functions
#========================================================================

def parse_PSSM(filename):
	"""
	# Input
	filename : string : filename of TXT file

	# Output
	PSSM_all : list of all PSSM
	"""
	f = open(filename, 'r')
	PSSM_all = []
	PSSM = []
	first_PSSM = True
	nuc = []

	l = f.readline()
	while(l):
		if l[0] != ";":
			if l[0] == "/":
				first_PSSM = False
				PSSM_all.append(PSSM)
				PSSM = []	
			else:
				weights = l.replace("\n","").split("\t")[1:len(l)]
				PSSM.append([int(w) for w in weights])
				if first_PSSM:
					nuc.append(l.split("\t")[0])
			
		l = f.readline()

	f.close()

	return(nuc, PSSM_all)



nuc, PSSM_all = parse_PSSM("./../Datas/oligo-analysis_2016-11-30.180333_2GFaRb_pssm_count_matrices.txt")

print(nuc)
for PSSM  in PSSM_all:
	print("     ")
	for l in PSSM:
		print(l)






















