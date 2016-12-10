# -*- coding: utf-8 -*-


"""
Projet REGULON - GENOM
Nika Abdollahi & Melissa Cardon

12-2016
"""

import numpy as np
import sklearn.cluster as skclust


#========================================================================
#                         Clusterings
#========================================================================
"""
Affinity Propagation can be interesting as it chooses the number of clusters based 
on the data provided. For this purpose, the two important parameters are the preference, 
which controls how many exemplars are used, and the damping factor.

The main drawback of Affinity Propagation is its complexity. 
The algorithm has a time complexity of the order O(N^2 T), 
where N is the number of samples and T is the number of iterations until convergence. 
Further, the memory complexity is of the order O(N^2) if a dense similarity matrix is used, 
but reducible if a sparse similarity matrix is used. 
This makes Affinity Propagation most appropriate for small to medium sized datasets.

http://stackoverflow.com/questions/35494458/affinity-propagation-in-python

Attention le input doit etre des distance et pas des similarite (diagonale = 0)
"""

## Fake distance matrix
sim_matrix = np.ones((10,10))
for i in range(10):
    for j in range(10):
        if j<i:
            rnumb = np.random.rand()
            sim_matrix[i,j] = rnumb
            sim_matrix[j,i] = rnumb

calcul_clusters_Aff_prop = skclust.AffinityPropagation(affinity='precomputed').fit(sim_matrix)
clusters_Aff_prop = calcul_clusters_Aff_prop.labels_