# -*- coding: utf-8 -*-

"""
Projet REGULON - GENOM
Nika Abdollahi & Melissa Cardon

12-2016
"""

import numpy as np
import sklearn.cluster as skclust
import math

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



#========================================================================
#                          Compare 2 columns
#========================================================================


def chi_2_col(col_X, col_Y):
    """
    col = np.array
    
    Returns chi2 value for 2 columns (formula ref_4)
    """
    nb_nuc = len(col_X)
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

##########################################################################

def SSD (col_X, col_Y):
	SSD2=0
	sum_col_X = sum(col_X)
	sum_col_Y = sum(col_Y)
	for Nucleotide in range(len(col_X)):
		SSD2+= math.pow(((float(col_X[Nucleotide])/sum_col_X)-(float(col_Y[Nucleotide])/sum_col_Y)),2)
	#print SSD2
	return 2-SSD2



##########################################################################
def PCC(col_X, col_Y):
    """
    input : PSSM with frequencies
    returns pearson correlation coef
    """
    nb_nuc = len(col_X)
    
    mu = 1/float(nb_nuc)
    col_X_centered = col_X - mu
    col_Y_centered = col_Y - mu
    numerateur = np.sum(col_X_centered * col_Y_centered )
    denominateur = (np.sum(col_X_centered**2))*(np.sum(col_Y_centered**2))
    
    return( numerateur / float(denominateur))


##########################################################################
def AKL(col_X, col_Y):
    """
    input : PSSM with frequencies (with pseudocount > no zero !!)
    returns Average Kullback–Leibler score (Ref4)
    """
    log_col_X = np.log(col_X)
    log_col_Y = np.log(col_Y)
    diff_log_XY = log_col_X - log_col_Y
    
    result = 10 - (np.sum(col_X * diff_log_XY) + np.sum(-col_Y * diff_log_XY))/float(2)
    
    return(result)
    


#========================================================================
#                         Compare 2 PSSM
#========================================================================




#========================================================================
#                          Use the functions
#========================================================================
nuc, PSSM_all = parse_PSSM("./Datas/Q1_PSSM.txt")

PSSM_all_freqs = PSSM_freqs(PSSM_all, 0.1)

#print(nuc)
print PSSM_all_freqs[0] # un PPSM 
print PSSM_all_freqs[1] # un autre PSSM ! 
"""
for PSSM  in PSSM_all:
	print("     ")
	for l in PSSM:
		print(l)


PSSM_all_psc = PSSM_pseudocount(PSSM_all, 0.1)
for PSSM  in PSSM_all_psc:
	print("     ")
	print(PSSM)
"""






A= PSSM_all_freqs[0] # un PPSM 
a = np.array(A)
for l in a:
	print l
print a[:,2]
B= PSSM_all_freqs[1] # un autre PSSM ! 
b= np.array(B)
for l in b:
	print l
print b[:,2]



col_X=a[:,2]
col_Y= b[:,2]

print("\nCompare 2 different columns")
print("SSD",SSD(col_X, col_Y))
print("PCC",PCC(col_X, col_Y))
print("AKL",AKL(col_X, col_Y))

print("\nCompare 2 identical columns")
print("SSD col_X",SSD(col_X, col_X))
print("SSD col_Y",SSD(col_Y, col_Y))
print("PCC col_X",PCC(col_X, col_X))
print("PCC col_Y",PCC(col_Y, col_Y))
print("AKL col_X",AKL(col_X, col_X))
print("AKL col_Y",AKL(col_Y, col_Y))




















