import sys
import numpy as np

########################################################################

def alignit(s1,s2,gap,subs):
	# matrice des distances
	m = range(len(s2)+1)
	# matrice des chemins
	for i in range(len(s2)+1):
		m[i] = range(len(s1)+1)
	# root cell
	m[0][0] = (0, 'o')
	# first line
	for j in range(1,len(s1)+1):
		#v=m[0][j-1][0] + gap
		m[0][j]=(0, 'g')
	# first column
	for i in range(1,len(s2)+1):
		m[i][0]=(0, 'h')	
	# tab
	for i in range(1,len(s2)+1):
		for j in range(1,len(s1)+1):
			distd = substitue(s2[i-1], s1[j-1],subs) + m[i-1][j-1][0]
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
def backtrack(m):
	#print "distance = ",m[len(s1)][len(s2)]
	#print m
	#print c
	# backtrack
	#nb col
	j = len(m[0])-1
	#nb lignes
	i = len(m)-1
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
	print sa
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
#								Main
########################################################################

subs=np.matrix([[3,0,1,0],[0,3,0,1],[1,0,3,0],[0,1,0,3]])

s1="AGCACACA"
s2="ACACACTA"


m=alignit(s1,s2,-1,subs)
print type(m)
sa=backtrack(m)
printali(s1,s2,sa)
printMatAli(matAli, s1, s2)