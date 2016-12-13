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
"""
affinity_matrix = np.load("./../Datas/Q2/Q2_affinity_matrix.npy")

clusters_Aff_prop = Aff_prop(affinity_matrix)

print(clusters_Aff_prop)
print(affinity_matrix[0:6,0:6])
"""



nuc, PSSM_all = parse_PSSM("./../Datas/oligo-analysis_2016-11-30.180333_2GFaRb_pssm_count_matrices.txt")
#PSSM_all_freqs = PSSM_freqs(PSSM_all, 0.1)

TF_Q1 = parse_PSSM_set("./../Datas/RegulonDB_PSSMSet.txt")
TF_Q1_f = PSSM_freqs_dict(TF_Q1, 0.1)

Bact = ["actinobacteria","cyanobacteria","firmicutes","proteobacteria"]
Prot = ["PhoA","PhoD","PhoX"]

PSSM_dict = create_dict_result2(Bact, Prot)
nb_PSSM = len(PSSM_dict)
print(nb_PSSM)

affinity_matrix = np.zeros((nb_PSSM,nb_PSSM))

for i in range(nb_PSSM):
	print(i)
	for j in range(nb_PSSM):
		if i<=j:
			affinity_matrix[i,j] = Score_Calculator(PSSM_dict[i][0],PSSM_dict[j][0],-10,"PCC")
			affinity_matrix[j,i] = affinity_matrix[i,j]


np.save("./../Datas/Q2/Q2_affinity_matrix_new.npy",affinity_matrix)
np.savetxt("./../Datas/Q2/Q2_affinity_matrix_new.txt",affinity_matrix)



#oligo-analysis_PhoA_actinobacteria_pssm_count_matrices.txt


"""
m = alignit2(PSSM1,PSSM2,-1,"PCC")
print(np.round(m,1))

m = alignit2(PSSM1,PSSM3,-1,"PCC")
print(np.round(m,1))

v= Score_Calculator(PSSM1,PSSM2,-1,"SSD")
print(v)
"""


"""
nuc, PSSM_all = parse_PSSM("./../Datas/oligo-analysis_2016-11-30.180333_2GFaRb_pssm_count_matrices.txt")
PSSM_all_freqs = PSSM_freqs(PSSM_all, 0.1)

TF_Q1 = parse_PSSM_set("./../Datas/RegulonDB_PSSMSet.txt")
TF_Q1_f = PSSM_freqs_dict(TF_Q1, 0.1)

dico_pssm1={}
g=0
for PSSM1 in PSSM_all_freqs:
	g+=1
	dico_Score={}
	for TF in TF_Q1_f.keys():
		PSSM2 = TF_Q1_f[TF]
		score = Score_Calculator(PSSM1,PSSM2,-10,"SSD")
		dico_Score[TF]=score
	dico_pssm1[g]=dico_Score
print "Metric : SSD"
#print dico_pssm1
for key in dico_pssm1.keys():
	d=dico_pssm1[key]
	#print max(d, key=d.get)
	print("For PSSM " +str(key)+" the best TF is " + str(max(d, key=d.get)) + " with score "+ str(d[max(d, key=d.get)]))

"""

"""
nuc, PSSM_all = parse_PSSM("./../Datas/oligo-analysis_2016-11-30.180333_2GFaRb_pssm_count_matrices.txt")
PSSM_all_freqs = PSSM_freqs(PSSM_all, 0.1)
PSSM_all_psc = PSSM_pseudocount(PSSM_all, 0.1)

TF_Q1 = parse_PSSM_set("./../Datas/RegulonDB_PSSMSet.txt")
TF_Q1_f = PSSM_freqs_dict(TF_Q1, 0.1)
TF_Q1_psc = PSSM_pseudocount_dict(TF_Q1, 0.1)

gap_penalty = -1
Results = []
#Metrics = ["SSD", "PCC", "AKL"]
Metrics = ["Chi2"]


for i in range(len(Metrics)):
	Metric = Metrics[i]
	res_all = []
	print("#####################" + Metric)

	for PSSM1 in PSSM_all_psc:

	#for PSSM1 in PSSM_all_freqs:
		res_PSSM = []
		best_score = 0
		best_TF = ""

		for TF in TF_Q1_psc.keys():
			PSSM2 = TF_Q1_psc[TF]
			score = Score_Calculator(PSSM1,PSSM2,gap_penalty,Metric)
			#print score
			res_PSSM.append((score,TF))
			if score > best_score:
				best_TF = TF
				best_score = score
			if TF == "PhoB":
				score_phob = score

		print("\nBest TF is " + str(best_TF) + " with score " + str(best_score) + " while PhoB has score = " +str(score_phob))
		res_PSSM_sorted = sorted(res_PSSM)[::-1]
		index_PhoB = [i for i in range(len(res_PSSM_sorted)) if res_PSSM_sorted[i][1] == "PhoB"]
		print( "PhoB is " + str(index_PhoB[0]+1) + " / " + str(len(res_PSSM_sorted)))
		res_all.append(res_PSSM_sorted)
		print(res_PSSM_sorted[0:5])


	Results.append(res_all)
"""

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