# -*- coding: utf-8 -*-
"""
TME8 SBAS: Suite du TME1, Coevolution des séquences

Melissa Cardon, Nathalie Lehmann & Chloe Quignot
"""
import matplotlib.pyplot as plt
import numpy as np
import random
import math

#==============================================================================
#                                 READ FILES
#==============================================================================
#import des donnees
def read_fasta(filename):
    """
    filename: name of the fasta file
    
    Returns a dictionary containing all the sequences and their names { name : sequence }
    """
    data = open(filename,'r')
    sequences = data.read()
    data.close()
    sequences = sequences.split('>')
    sequences = sequences[1:len(sequences)]

    # dictionary of sequences
    dict_DNA = {}
    for seq in sequences:
        name = str(seq[0:seq.find('\n')])
        DNA = str(seq[seq.find('\n'):len(seq)])
        DNA = DNA.replace('\n','')
        dict_DNA[name] = DNA
    
    return(dict_DNA)


def read_distance(filename, nb_aa_total):
    """
    filename : string, name of the file containig the distances between each amino acid pair
    nb_aa_total : int, length of the sequences (number of amino acids in each sequence)
    
    Returns np.array distance matrix
    """
    data = open(filename,'r')
    distances_file = data.read()
    data.close()
    distances_file = distances_file.split('\n')
    
    mat_dist = np.zeros((nb_aa_total,nb_aa_total))
    
    for i in range(len(distances_file)):
        distline = distances_file[i]
        dist = distline.split(' ')
        if len(dist) == 3:
            mat_dist[dist[0],dist[1]] = dist[2]
            mat_dist[dist[1],dist[0]] = dist[2]

        else:
            print("line " + str(i) + " : " + str(distline) + " empty : not added to the distance matrix") #i.e. empty line

    return mat_dist


#==============================================================================
#                                   FUNCTIONS
#==============================================================================
# creation de la matrice de comptage
def matrice_ni(dict_aa2nb, sequences):
    """
    dict_aa2nb : { amino acid : associated index }
    sequences : dictionary of sequences { name : sequence } i.e. Dtrain
    
    Returns the matrix (lines = amino acid, column = position of the amino acid in the sequence) containing the number of each amino acid at each position in the sequences
    """
    #initialisation of the matrix
    random_key = random.choice(sequences.keys())
    matrix_ni = np.zeros(( len(dict_aa2nb),len(sequences[random_key]) ))
    
    #counting the amino acids at each position
    for seq in sequences.values() :
        for i in range(len(seq)):
            aa = seq[i]
            matrix_ni[dict_aa2nb[aa],i] +=1
    
    return(matrix_ni)


# creation de la matrice PSWM
def matrice_wi(dict_aa2nb, sequences, ps):
    """
    dict_aa2nb : { amino acid : associated index }
    sequences : dictionary of sequences { name : sequence } i.e. Dtrain
    ps: pseudo-count, value added to each count, e.g. 1
    
    Returns the PSWM matrix (normalised ni matrix with pseudo-count)
    """
    matrix_wi = matrice_ni(dict_aa2nb, sequences)
    matrix_wi += ps
    div = np.sum(matrix_wi, axis=0)
    matrix_wi /= div
    
    return(matrix_wi)


# creation de la matrice de comptage avec correlation
def matrice_nij(combi_aa2nb, combi_ij2idx, sequences):
    """
    combi_aa2nb: { (aa1, aa2) : line index } dictionary that associates a pair of amino acids to a specific line in the nij matrix
    combi_ij2idx:  { (i, j) : column index } dictionary that associates a pair of positions i and j to a specific column in the nij matrix (NB: i < j and i != j)
    sequences: dictionary of sequences { name : sequence } i.e. Dtrain
    
    Returns the matrix (lines = all possible combinations of 2 amino acids, column = position of the amino acid pairs in the sequence) containing the number of each amino acid pair at each position in the sequences
    """    
    #initialisation of nij
    nij = np.zeros(( len(combi_aa2nb), len(combi_ij2idx) ))    
    
    #filling in the matrix
    for seq in sequences.values():
        liste_i = [0]
        for j in range(1,len(seq)):
            for i in liste_i: # i < j
                nij[combi_aa2nb[(seq[i],seq[j])], combi_ij2idx[(i,j)]] += 1
            liste_i.append(j)       
        
    return nij


# creation de la matrice PSWM avec correlation
def matrice_wij(combi_aa2nb, combi_ij2idx, sequences, q):
    """
    combi_aa2nb: { (aa1, aa2) : line index } dictionary that associates a pair of amino acids to a specific line in the nij matrix
    combi_ij2idx:  { (i, j) : column index } dictionary that associates a pair of positions i and j to a specific column in the nij matrix (NB: i < j and i != j)
    sequences: dictionary of sequences { name : sequence } i.e. Dtrain
    q: pseudo-count /!\ it has to be the same q as in matrice_wi() where q = 1*len(dict_aa2nb), i.e 21, in order to guarantee that wi(a) = Sum over b { wij(a,b) }
    
    Returns the PSWM matrix (normalised nij matrix with pseudo-count) with correlation (lines = all possible combinations of 2 amino acids, column = position of the amino acid pairs in the sequence)
    """    
    wij = matrice_nij(combi_aa2nb, combi_ij2idx, sequences)
    M = np.sum(wij, axis=0)       #sum of each column before adding the pseudo-count
    wij += 1./q                   #we add the pseudo-count  
    wij = wij/(M+q)               #we normalise

    return wij


def Mij_information_mutuelle(wi, wij, dict_aa2nb, combi_aa2nb, combi_ij2idx):
    """
    """
    nbpos = len(wi[0])
    M = np.zeros((nbpos, nbpos))
    
    for (i,j) in combi_ij2idx.keys(): #i = pos of aa a, j = pos of aa b
        for a,ia in dict_aa2nb.items():
            wi_a = wi[ia,i]
            for b,jb in dict_aa2nb.items():
                wj_b = wi[jb,j]
                wij_ab = wij[combi_aa2nb[(a,b)], combi_ij2idx[(i,j)]]
                M[i,j] += wij_ab * math.log(wij_ab/(wi_a*wj_b),2)
                M[j,i] = M[i,j]
    
    return M


def fraction_close_contact(dist_mtx, dist_min):
	count_all = 0.
	count_close = 0.
	for i in range(len(dist_mtx)):
		for j in range(len(dist_mtx)):
			if j >i :
				count_all +=1
				if dist_mtx[i,j] < dist_min :
					count_close +=1

	return count_close / float(count_all)


def sortMij(Mij, combi_ij2idx, dist_mtx, nb):
    """
    Mij: mutal information matrix for position i & j in the sequences
    combi_ij2idx: { (i, j) : column index } dictionary that associates a pair of positions i and j to a specific column in the nij matrix (NB: i < j and i != j)
    dist_mtx: distances between positions i and j, imported from the file distances.txt    
    nb: number of best positions (i,j) to extract from the set
    
    Returns a table of the 'nb' positions [i,j] with the best M value (from highest to lowest), their distances and M values
    """
    Mij_pos = combi_ij2idx.keys()   #list of all positions (i,j)
    
    M_list = []                     #list of all M values corresponding to the positions in Mij_pos
    for (i,j) in Mij_pos:
        M_list.append(Mij[i,j])
    
    sorted_idx = list(reversed(np.argsort(M_list))) #sorted_idx[:50] = idx of the 50 positions with the highest M value
    
    Mpos = np.array(Mij_pos)[sorted_idx[:nb]]     #positions [i,j] with the nb best M values
    Mval = np.array(M_list)[sorted_idx[:nb]]      #nb best M values (highest to lowest)
    
    dist = []                                     #distance between positions i and j
    for [i,j] in Mpos:
        dist.append(dist_mtx[i,j])    
    
    table = np.zeros((4,nb))
    table[0] = Mpos[:,0]  # position i
    table[1] = Mpos[:,1]  # position j
    table[2] = dist       # distance between i and j
    table[3] = Mval       # M value of (i,j)
    
    return table 


def calc_fraction(dist, lim_dist):
    """
    dist: list of distances between positions i and j aranged in the right order (i.e. from highest to lowest M value associated)
    lim_dist: a distance, all (i,j) with distances under this limit are considered "neighbour positions"
    
    Returns the fraction of distances under 'lim_dist' for increasing numbers of position pairs taken into account (i.e. from only one pair to len(dist) pairs)
    """    
    tot = len(dist)
    cumulated_fraction = np.zeros(tot)
    
    for i, dist in enumerate(dist):
        if dist < lim_dist:
            cumulated_fraction[i] = 1
    
    #fraction of positions within the 1, 2, 3, ..., 'tot' best that have a distance under 'lim_dist'    
    cumulated_fraction = np.cumsum(cumulated_fraction)
    cumulated_fraction = cumulated_fraction/ range(1,tot+1)
    
    return cumulated_fraction



def fraction_graph(cum_frac, overall_contact):
    """
    cum_frac: fraction of distances under 8 for increasing numbers of position pairs taken into account, as returned by calc_fraction()
    
    Plots the fraction of positions with a distance under 8 vs the number of position pairs considered
    """
    plt.figure()
    plt.plot(np.arange(1,len(cum_frac)+1), cum_frac)
    plt.axhline(y=overall_contact, color = 'k')
    plt.ylim([0.0,1.005])
    plt.title("Fraction of positions with a distance under 8\nvs number of position pairs considered")
    plt.xlabel("total number of top position pairs (i.e. pairs with best M value)\ntaken into account")
    plt.ylabel("fraction of pairs that have a distance\nunder 8 (i.e. neighbours) within the top\nposition pairs taken into account") 
    plt.show()
    plt.close()
    
    return


def savewij(filename, wij, combi_nb2aa, combi_idx2ij):
    """
    filename : name of the text file the results are to be saved in (ex: results.txt)
    wij: coevolution PSWM matrix 
    combi_nb2aa: { line index : (aa1, aa2) } dictionary that associates a pair of amino acids to a specific line in the nij matrix
    combi_idx2ij: { column index : (i, j) } dictionary that associates a pair of positions i and j to a specific column in the nij matrix (NB: i < j and i != j)

    Creates a text file named filename containing the wij matrix.
    """
    file = open(filename,"w")
    for i in range(len(combi_idx2ij)):
        file.write("\t"+str(combi_idx2ij[i][0])+","+str(combi_idx2ij[i][1]))
    for i in range(len(combi_nb2aa)):
        file.write("\n"+combi_nb2aa[i][0]+" "+combi_nb2aa[i][1])
        for j in range(len(combi_idx2ij)):
            file.write("\t%.3f" %(wij[i,j]))
    file.close()
    
    return

#==============================================================================
#                           IMPORTING THE DATA
#==============================================================================

sequences = read_fasta('Dtrain.txt')               #dicitonary of sequences (of same length) { name : sequence }
lgseq = len(sequences[sequences.keys()[0]])        #length of a sequence
dist_mtx = read_distance('distances.txt', lgseq)   #distance matrix

#==============================================================================
#                       CREATION OF THE DICTIONARIES
#==============================================================================

dict_nb2aa = {1:'A', 2:'C', 3:'D', 4:'E', 5:'F', 6:'G', 7:'H', 8:'I', 9:'K', 10:'L', 11:'M', 12:'N', 13:'P', 14:'Q', 15:'R', 16:'S', 17:'T', 18:'V', 19:'W', 20:'Y', 0:'-'}
dict_aa2nb =  {v: k for k, v in dict_nb2aa.items()}

combi_nb2aa = {} #{ idx: (aa1,aa2) }
k = 0
for i, aa1 in dict_nb2aa.iteritems():
    for j, aa2 in dict_nb2aa.iteritems():
        combi_nb2aa[k] = (aa1,aa2)
        k += 1

combi_aa2nb = {v: k for k, v in combi_nb2aa.items()} #{ (aa1,aa2): idx }

combi_idx2ij = {} #{ idx: (i,j) } with i < j & i != j
k = 0
for i in range(lgseq):
    for j in range(lgseq):
        if i < j:
            combi_idx2ij[k] = (i,j)
            k += 1
    
combi_ij2idx = {v: k for k, v in combi_idx2ij.items()} #{ (i,j): idx }


#==============================================================================
#                           USING THE FUNCTIONS
#==============================================================================
ps = 1                  #pseudo-count
q = ps*len(dict_aa2nb)  #weight

#wi matrix:
wi = matrice_wi(dict_aa2nb, sequences, ps)

np.savetxt("wi_matrix.txt",wi)

#wij matrix:
wij = matrice_wij(combi_aa2nb, combi_ij2idx, sequences, q)

#to save the wij matrix in a text file:
savewij("wij.txt", wij, combi_nb2aa, combi_idx2ij)

#to check wij:
print "\nw0(-): %.4f" % wi[0,0]
print "∑aa w0,0(-,aa): %.4f\n" % np.sum(wij[:21,0])   #wi[0,0] = np.sum(wij[:21,i]) for all i in [0,48[

print "w1(-): %.4f" % wi[0,1]
print "∑aa w1,0(-,aa): %.4f\n" % np.sum(wij[:21,48])  #wi[0,1] = np.sum(wij[:21,i]) for all i in [48,96[

print "w0(A): %.4f" % wi[1,0]
print "∑aa w0,0(A,aa): %.4f\n" % np.sum(wij[21:42,0]) #wi[0,1] = np.sum(wij[21:42,i]) for all i in [0,48[

#etc.

#Mij matrix containing the M values of each pair of position i,j in the sequences
Mij = Mij_information_mutuelle(wi, wij, dict_aa2nb, combi_aa2nb, combi_ij2idx)
print "M(0,1): %.3f" % Mij[0,1]
np.savetxt("Mij_matrix.txt", Mij)

#sorted_pos is an array (4, 50) containing 50 positions (i,j) with the best M value (from highest to lowest):
#1st line = position i
#2nd line = position j
#3rd line = distance between the aa at positions i and j
#4th line = M value associated to the position pair (i,j)
sorted_pos = sortMij(Mij, combi_ij2idx, dist_mtx, 50)

#fraction of position pairs with a distance under 8 and for an increasing number of pairs taken into account
cum_frac = calc_fraction(sorted_pos[2], 8)

# for pair_i in range(len(sorted_pos[0])):
# 	print(pair_i,sorted_pos[:,pair_i], cum_frac[pair_i])



overall_contact = fraction_close_contact(dist_mtx, 8)
print("Il y a en moyenne " + str(overall_contact*100) + "pourcent de positions en contact ")
#to plot the graph fraction vs number of pairs
fraction_graph(cum_frac, overall_contact)

print(Mij.shape)

np.save("Mij_matrix.npy",Mij)
np.save("Dist_matrix.npy",dist_mtx)

plt.pcolor(Mij, cmap=plt.cm.BuPu)
plt.colorbar()
plt.title("Mutual information")
plt.xlabel("position in linear peptide")
plt.ylabel("position in linear peptide") 
plt.show()

plt.pcolor(dist_mtx, cmap=plt.cm.BuPu_r)
plt.colorbar()
plt.title("Distance in folded protein")
plt.xlabel("position in linear peptide")
plt.ylabel("position in linear peptide") 
plt.show()


CR_Mij = (Mij - np.mean(Mij)) / np.std(Mij)
CR_dist_mtx = ( dist_mtx - np.mean(dist_mtx) ) / np.std(dist_mtx)
Comp = CR_Mij - CR_dist_mtx

plt.pcolor(Comp, cmap=plt.cm.PiYG_r)
plt.colorbar()
plt.title("Mutual info - Distance (normalised values)")
plt.xlabel("position in linear peptide")
plt.ylabel("position in linear peptide") 
plt.show()







