# -*- coding: utf-8 -*-

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


def PSSM_pseudocount(PSSM_all, pseudocount):
    """
    pseudocount = float = number to add for pseudocount
    returns PSSM_all with pseudocounts
    """
    PSSM_all_psc = []
    for PSSM in PSSM_all:
        npPSSM = np.array(PSSM)
        PSSM_all_psc.append(npPSSM + pseudocount)
    
    return(PSSM_all_psc)





#========================================================================
#                          Compare 2 PSSM
#========================================================================


def chi_2_col(col_X, col_Y):
    """
    col = np.array
    
    Returns chi2 value for 2 columns (formula ref_4)
    """
    nb_nuc = col_X.shape()[0]
    N = sum(col_X) + sum(col_Y)
    chi2col = 0
    
    for col in [col_X,col_Y]:
        nk = sum(col)
        
        for nuc in range(nb_nuc):
            nkb = col[nuc]
            nekb = nk*nkb/float(N)
            chi2_col_nuc = ((nkb - nekb)**2)/float(nekb)
            chi2col += chi2_col_nuc
            
    return(chi2col)








nuc, PSSM_all = parse_PSSM("./Datas/Q1_PSSM.txt")

print(nuc)
for PSSM  in PSSM_all:
	print("     ")
	for l in PSSM:
		print(l)


PSSM_all_psc = PSSM_pseudocount(PSSM_all, 0.1)
for PSSM  in PSSM_all_psc:
	print("     ")
	print(PSSM)























