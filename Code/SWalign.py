# -*- coding: utf-8 -*-

"""
Projet REGULON - GENOM
Nika Abdollahi & Melissa Cardon

12-2016
"""

import numpy as np
import sys
import math


#========================================================================
#                    Metrics to compare 2 columns
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
	for Nucleotide in range(len(col_X)):
		SSD2+= math.pow((col_X[Nucleotide]-col_Y[Nucleotide]),2)
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
    if numerateur == 0:
    	return 0
    else:
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
    
#<<<<<<< Updated upstream




#========================================================================
#                    Choice of metrics
#========================================================================
##########################################################################
def ChooseMetric(col_X, col_Y,metric):
    """
    Metric : possible choices : SSD , PCC, AKL
    """
    if metric =="SSD":
        res= SSD(col_X, col_Y)
    elif metric== "PCC":
        res= PCC(col_X, col_Y)
    elif metric== "AKL":
        res= AKL(col_X, col_Y)
    else:
        print "ChooseMetric :: Error, model not found"
    return res




#========================================================================
#                         Compare 2 PSSM
#========================================================================

########################################################################

def alignit(PSSM1,PSSM2,gap,Metric):
	# matrice des distances
	m = range(len(PSSM2[0])+1)
	#print m
	# matrice des chemins
	for i in range(len(PSSM2[0])+1):
		m[i] = range(len(PSSM1[0])+1)
	#print m
	# root cell
	m[0][0] = (0, 'o')
	# first line
	for j in range(1,len(PSSM1[0])+1):
		#v=m[0][j-1][0] + gap
		m[0][j]=(0, 'g')
	# first column
	for i in range(1,len(PSSM2[0])+1):
		m[i][0]=(0, 'h')
	# tab

	for i in range(1,len(PSSM2[0])+1):
		for j in range(1,len(PSSM1[0])+1):
			distd = ChooseMetric(PSSM2[:,i-1], PSSM1[:,j-1],Metric) + m[i-1][j-1][0]
			disth = gap + m[i-1][j][0]
			distg = gap + m[i][j-1][0]
			if distd >= disth and distd >= distg:
				m[i][j] = (distd,'d')
#				c[i][j] = 'd'  # substitution
			elif disth >= distd and disth >= distg:
				m[i][j] = (disth, 'h')
#				c[i][j] = 'h' # insertion
			else:
				m[i][j] = (distg, 'g')
#				c[i][j] = 'g' # deletion
	return (m)

########################################################################
def substitue(res1, res2,subs):
	score= subs[findIndice(res1),findIndice(res2)]
	return score

########################################################################

def findIndice(res):
	if res=="A":
		Indice=0
	elif res =="T":
		Indice=1
	elif res=="G":
		Indice=2
	else:
		Indice =3
	return Indice
########################################################################
#Retrouver le chemin
def backtrack(m, imax,jmax):
	#print "distance = ",m[len(s1)][len(s2)]
	#print m
	#print c
	# backtrack
	#nb col
	j = jmax
	#len(m[0])-1
	#nb lignes
	i =imax
	#len(m)-1
	chemin = ''
	while (i != 0 or j != 0):
		if i < 0 or j < 0:
			print "backtrack:: ERROR i or j <0",i,j
			exit(1)
		#print i,j,m[i][j]
		if m[i][j][1] == 'd':
			i = i-1
			j = j-1
			chemin = 'd'+chemin
		elif m[i][j][1] == 'g':
			if len (m[i][j]) ==3:
				#print "Saut de ",m[i][j][2]
				for k in xrange(m[i][j][2]):
					j=j-1 
					chemin = 'g'+chemin
			else:
				j = j-1
				chemin = 'g'+chemin
		else:
			if len (m[i][j]) ==3:
				#print "Saut de ",m[i][j][2]
				for k in xrange(m[i][j][2]):
					i = i-1
					chemin = 'h'+chemin
			else:
				i = i-1
				chemin = 'h'+chemin
		#print i,j
	#print chemin
	return chemin


########################################################################
## Print aligment
def printali(s1,s2,sa):
	ii = 0
	#print sa
	for i in range(len(sa)):
		if sa[i] == 'd' or sa[i] == 'g':
			sys.stdout.write(s1[ii])
			ii = ii + 1
		else:
			sys.stdout.write('-')
	print ''
	ii = 0
	for i in range(len(sa)):
		if sa[i] == 'd' or sa[i] == 'h':
			sys.stdout.write(s2[ii])
			ii = ii + 1
		else:
			sys.stdout.write('-')
	print ''
########################################################################
def solve(mat):
    print " ", " ".join([str(x) for x in xrange(len(mat))])
    for i, x in enumerate(mat):
        print i, " ".join([str(y) for y in x])


########################################################################
def FindIndiceMax(matrix):
	Nb_line=len(matrix)
	Nb_colomn=len(matrix[0])
	Valmax=0
	imax=0
	jmax=0
	for i in range(Nb_line) :
		for j in range(Nb_colomn):
			#print matrix[i][j][0] 
			if matrix[i][j][0] >Valmax:
				Valmax=matrix[i][j][0]
		 		imax=i
		 		jmax=j
	return Valmax,imax,jmax
########################################################################
def Score_Calculator(PSSM1,PSSM2,gap,Metric):
	m= alignit(PSSM1,PSSM2,gap,Metric)
	v,imax,jmax=FindIndiceMax(m)
	sa=backtrack(m,imax,jmax)
	return v


########################################################################
def listCombinaton(listPSSM):
    #print list(range(0, len(listPSSM)))
    listoflist=[]
    for l in range(len(listPSSM)):
        for k in range(l,len(listPSSM)):
            listoflist.append([listPSSM[l],listPSSM[k]])
    #print listoflist
    return listoflist
########################################################################
def Matrix_Score(listPSSM, Metric):
    """
    Metric : possible choices : SSD , PCC, AKL
    """
    n=len(listPSSM)
    Matrix_Score=np.zeros((n,n))
    List_Of_PSSM_Ind=listCombinaton(list(range(0, len(listPSSM))))
    for Deux_PSSM in List_Of_PSSM_Ind:
        A=listPSSM[Deux_PSSM[0]]
        B=listPSSM[Deux_PSSM[1]]
        if len(A)>len(B):
            score= Score_Calculator(A,B,-1,Metric)
        else :
            score= Score_Calculator(B,A,-1,Metric)
        #print Deux_PSSM,score
        Matrix_Score[(Deux_PSSM[0]),(Deux_PSSM[1])]=score
        Matrix_Score[(Deux_PSSM[1]),(Deux_PSSM[0])]=score
    return  Matrix_Score


"""
def Matrix_Score(listPSSM):
    n=len(listPSSM)
    Matrix_Score=np.zeros((n,n))

    List_Of_PSSM_Ind=listCombinaton(list(range(0, len(PSSM_all_freqs))))
    for Deux_PSSM in List_Of_PSSM_Ind:
        A=listPSSM[Deux_PSSM[0]]
        B=listPSSM[Deux_PSSM[1]]
        if len(A)>len(B):
            score= alignit(A,B,-1,"SSD")
        else :
            score= alignit(B,A,-1,"SSD")
        Matrix_Score[Deux_PSSM[0]-1,Deux_PSSM[1]-1]=score
    print Matrix_Score



def Normalise_Mtrix_Score(MatixScore):
"""



"""
########################################################################
#								Main
########################################################################

print Matrix_Score(PSSM_all_freqs)
"""
