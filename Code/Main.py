# -*- coding: utf-8 -*-


"""
Projet REGULON - GENOM
Nika Abdollahi & Melissa Cardon

12-2016
"""


from f_import2 import *
from SWalign import *
from f_cluster import *




#========================================================================
#                          Use the functions
#========================================================================
#nuc, PSSM_all = parse_PSSM("./Datas/Q1_PSSM.txt")


##########################################################################
#>>>>>>> Stashed changes

#nuc, PSSM_all = parse_PSSM("./Datas/Q1_PSSM.txt")
nuc, PSSM_all = parse_PSSM("./../Datas/oligo-analysis_2016-11-30.180333_2GFaRb_pssm_count_matrices.txt")
#nuc, PSSM_all = parse_PSSM("./../Datas/test.txt")
PSSM_all_freqs = PSSM_freqs(PSSM_all, 0.1)


TF_Q1 = parse_PSSM_set("./../Datas/RegulonDB_PSSMSet.txt")
TF_Q1_f = PSSM_freqs_dict(TF_Q1, 0.1)



for PSSM1 in PSSM_all_freqs:
	best_score = 0
	best_TF = ""

	for TF in TF_Q1_f.keys():
		PSSM2 = TF_Q1_f[TF]
		score = Score_Calculator(PSSM1,PSSM2,-1,"SSD")
		if score > best_score:
			best_TF = TF
			best_score = score

	print("Best TF is" + str(best_TF))


"""
Verification de l'import sous forme de dictionnaire
count=0
for key in TF_Q1_f.keys():
	count +=1
	print(key)
	print(np.sum(TF_Q1_f[key],0))
	for l in TF_Q1_f[key]:
		print(np.round(l,3))

print(count)
"""


"""
print("\n SSD Metric")
Metric = "SSD"
print(Matrix_Score(PSSM_all_freqs,Metric))


print("\n PCC Metric")
Metric = "PCC"
print(Matrix_Score(PSSM_all_freqs,Metric))


print("\n AKL Metric")
Metric = "AKL"
print(Matrix_Score(PSSM_all_freqs,Metric))
"""

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


"""













