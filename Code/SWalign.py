import sys
import numpy as np
from f_import2 import *

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
def Matix_Score(listPSSM):
    n=len(listPSSM)
    Matix_Score=np.zeros((n,n))
    List_Of_PSSM_Ind=listCombinaton(list(range(0, len(PSSM_all_freqs))))
    for Deux_PSSM in List_Of_PSSM_Ind:
        A=listPSSM[Deux_PSSM[0]]
        B=listPSSM[Deux_PSSM[1]]
        if len(A)>len(B):
            score= Score_Calculator(A,B,-1,"SSD")
        else :
            score= Score_Calculator(B,A,-1,"SSD")
        #print Deux_PSSM,score
        Matix_Score[(Deux_PSSM[0]),(Deux_PSSM[1])]=score
        Matix_Score[(Deux_PSSM[1]),(Deux_PSSM[0])]=score
    return  Matix_Score
########################################################################
#								Main
########################################################################

print Matix_Score(PSSM_all_freqs)
