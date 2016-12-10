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


def PSSM_freqs(PSSM_all, pseudocount):
    """
    returns PSSM with each col as freq instead of count
    (with pseudocounts)
    """
    PSSM_all_psc = PSSM_pseudocount(PSSM_all, pseudocount)
    
    PSSM_all_f = []
    for PSSM in PSSM_all_psc:
        PSSM_colsums = np.sum(PSSM,0,dtype='float')
        PSSM_all_f.append(PSSM / PSSM_colsums)
    
    return(PSSM_all_f)
















