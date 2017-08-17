### Peyton Greenside
### Script with various functions to cluster grammars
### 4/11/16
#################################################################################

import os, sys
import numpy as np
from scipy.signal import correlate2d
import pandas as pd

import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

import pickle
scriptsDir = os.environ.get("UTIL_SCRIPTS_DIR");
if (scriptsDir is None):
    raise Exception("Please set environment variable UTIL_SCRIPTS_DIR");
sys.path.insert(0,scriptsDir);
import pathSetter;
import fileProcessing as fp;
from importDataPackage import importData;
import importDataPackage
from synthetic import synthetic as sn;
import keras;
scriptsDir = os.environ.get("ENHANCER_SCRIPTS_DIR");
if (scriptsDir is None):
    raise Exception("Please set environment variable ENHANCER_SCRIPTS_DIR");
sys.path.insert(0,scriptsDir+"/featureSelector/deepLIFFT/");
from plottingUtilitiesPackage import matplotlibHelpers as mplh;
import scipy.misc
import h5py
import util

import copy
import random
import time
import yaml
import argparse
from scipy.signal import correlate2d
import copy
from collections import OrderedDict

#  sys.path.insert(0,scriptsDir+"/featureSelector/deepLIFFT/kerasBasedBackprop");
#  import deepLIFTonGPU
#  import deepLIFTutils
#  from deepLIFTonGPU import ScoreTypes, Activations_enum, OutLayerInfo, getScoreFunc
import criticalSubsetIdentification as csi
# sys.path.insert(0,"/users/pgreens/git/xcor");
# import xcor;

# sys.path.append('/users/pgreens/git/disease_projects/deltaDeepLift/')
# import deltaDeepLift_w_grammarDetection_Graph
# reload(deltaDeepLift_w_grammarDetection_Graph)
# from deltaDeepLift_w_grammarDetection_Graph import agglomerative_clustering

allowedMethods = ['kmeans_w_agg', 'tsne', 'iter_kmeans_w_agg', 'kmedoids_w_agg', 'DBSCAN', 'phenograph']

def grammarDetection(seqlets, seqletCorrMat, flankSize=0, grammarClusterMethod='kmedoids_w_agg',
                     cc_threshold=0.9, return_seqlet_indices=False, trackNameToUse='multipliers',
                     trim=True, n_init=10, **kwargs):
    ### Use t-SNE for getting labels
    if grammarClusterMethod == 'tsne':
        mergedSeqlets = cluster_by_tsne(seqlets, subtracksToInclude, accountForRevComp)
        cluster_label = 'cluster_by_tsne'
    # Iterative k-means starting from small number of clusters, followed by agglomerative clustering
    elif grammarClusterMethod == 'iter_kmeans_w_agg':
        (mergedSeqlets, seqlet_indices) = cluster_by_iterative_kmean(seqlets, seqlets, cc_threshold)
        cluster_label = 'cluster_by_iterative_kmeans'
    # Iterative k-medoids
    elif grammarClusterMethod == 'kmedoids_w_agg':
        (mergedSeqlets, seqlet_indices) = cluster_by_kmedoids(np.max(seqletCorrMat)-seqletCorrMat, 
                                                              seqlets, 
                                                              cc_threshold=cc_threshold, 
                                                              trackNameToUse=trackNameToUse,
                                                              n_init=n_init)
        cluster_label = 'cluster_by_kmedoids'
    # agglomerative clustering for labels
    elif grammarClusterMethod == 'kmeans_w_agg':
        (mergedSeqlets, seqlet_indices) = cluster_by_agg(seqletCorrMat, seqlets)
        # First cluster by Kmeans
        cluster_label = 'cluster_by_agg'
    elif grammarClusterMethod == 'DBSCAN':
        (mergedSeqlets, seqlet_indices) = cluster_by_DBSCAN(seqletCorrMat, seqlets, min_samples=10, cc_threshold=0.92)
        cluster_label = 'cluster_by_DBSCAN'
    elif grammarClusterMethod == 'phenograph':
        (mergedSeqlets, seqlet_indices) = cluster_by_phenograph(seqletCorrMat, seqlets, min_samples=10, cc_threshold=None, **kwargs)
        cluster_label = 'cluster_by_phenograph'
    else:
        raise Exception("Use only allowed clustering methods in grammarClustering.py")
    # trimmingFunc = csi.TrimArrayColumnsToNumUnderlyingObs(0.2) # takes sequences iwth 
    reload(csi)

    if trim:
        min_seqlet_size = np.min([s.summedCoreDeepLIFTtrack.shape[-1] for s in seqlets])

        trimmingFunc = csi.TrimArrayColumnsToPeak(slidingWindowSizeForPeak=min_seqlet_size, flanksToExpand=flankSize,
                                                  trackNameToUse=trackNameToUse, useRangeNotSum=True)
        #once again, subtracksToInclude indicates the subtracks to consider for merging. Should be
        #the same as what you supplied for the cross-correlation
        mergedGrammars = csi.adjustGrammarsUsingTrimmingCriterion(
                            {i: seqlet for (i, seqlet) in enumerate(mergedSeqlets)},
                             trimmingFunc=trimmingFunc);

    else:
        mergedGrammars = mergedSeqlets

    if return_seqlet_indices:
        return (mergedGrammars, seqlet_indices)
    else:
        return mergedGrammars

# Function from CSFoo
# Agglomerative clustering to find grammars
def agglomerative_clustering(grammars_list, gradient_track, cc_threshold, trackNameToUse='coreDeepLIFTtrack'):
    indices_list = [set([i]) for i in range(len(grammars_list))]
    while True:
        # The PerPosNormFuncs are used for ensuring that the corrrelations are
        # between -1 and 1
        grammars_cc = csi.getCorrelationMatrix(
            grammars_list,
            subtracksToInclude=[trackNameToUse],
            accountForRevComp=True,
            numThreads=1,
            secondsBetweenUpdates=6,
            xcorBatchSize=10,
            smallerPerPosNormFuncs=[],
            largerPerPosNormFuncs=[])

        np.fill_diagonal(grammars_cc, 0)
        max_grammars_cc = np.max(grammars_cc)

        print("max cc: ", max_grammars_cc)
        if max_grammars_cc < cc_threshold:
            break

        max_cc_idx1, max_cc_idx2 = \
            np.unravel_index(np.argmax(grammars_cc), grammars_cc.shape)

        merged_grammar = grammars_list[max_cc_idx1].merge(
            grammars_list[max_cc_idx2],
            subtracksToInclude=[gradient_track],
            subtrackNormaliseFunc=util.CROSSC_NORMFUNC.meanAndTwoNorm,
            normaliseFunc=util.CROSSC_NORMFUNC.none,
            smallerPerPosNormFuncs=[],
            largerPerPosNormFuncs=[],
            revComp=True)

        removeset = {max_cc_idx1, max_cc_idx2}
        merged_indices = \
            indices_list[max_cc_idx1].union(indices_list[max_cc_idx2])
        indices_list = [indices_set
                        for i, indices_set in enumerate(indices_list)
                        if i not in removeset]
        indices_list.append(merged_indices)
        grammars_list = [grammar
                         for i, grammar in enumerate(grammars_list)
                         if i not in removeset]
        grammars_list.append(merged_grammar)

    return grammars_list, indices_list


# hierarchical clustering
def cluster_by_hier(distMatrix, seqletDictTask, linkage_method='complete', cc_threshold=0.93, trackNameToUse='multipliers'):
    # import scipy.cluster.hierarchy as hac
    # z = hac.linkage(distMatrix, method=linkage_method)

    import fastcluster
    z = fastcluster.linkage(distMatrix, method=linkage_method) ## !!! WAY WAY TOO SLOW

    cluster_idx = scipy.cluster.hierarchy.fcluster(z, t=4)

    # ! Stopped here
    # Merge seqlets from Kmeans clustering
    seqlet_list = csi.createMergedGrammars(cluster_idx,
                                       seqletDictTask,
                                       subtracksToInclude=[trackNameToUse],
                                       accountForRevComp=True)

    # Prune overly similar ones
    (mergedSeqlets, indices_list) = agglomerative_clustering(grammars_list=seqlet_list.values()
                                    , gradient_track=trackNameToUse, cc_threshold=cc_threshold)

    # Get the index of the original seqlets
    reverse_indices_list = {val: key for (key, listy) in enumerate(indices_list) for val in listy}
    seqlet_indices = [reverse_indices_list[idx] for idx in cluster_idx]

    return (mergedSeqlets, seqlet_indices)

# Kmedoids clustering
def cluster_by_kmedoids(distMatrix, seqletDictTask, cc_threshold=0.93, trackNameToUse='multipliers', n_init=10):
    KMEDOIDS_NUM_CLUSTERS = 300
    sys.path.append('/users/pgreens/git/enhancer_prediction_code/featureSelector/deepLIFFT/Applied/Chromputer/DnaseMnaseSequence')
    import kmedoids
    reload(kmedoids)
    (cluster_idx, best_medoids, min_cost) = kmedoids.kmedoids(distMatrix, n_clusters=KMEDOIDS_NUM_CLUSTERS
                                                            , n_init=n_init, max_iter=300, init='kmeans++'
                                                            , random_search_tries=3)


    # Merge seqlets from Kmeans clustering
    seqlet_list = csi.createMergedGrammars(cluster_idx,
                                           seqletDictTask,
                                           subtracksToInclude=[trackNameToUse],
                                           accountForRevComp=True)

    # Prune overly similar ones
    (mergedSeqlets, indices_list) = agglomerative_clustering(grammars_list=seqlet_list.values(), 
                                                gradient_track=trackNameToUse, cc_threshold=cc_threshold, trackNameToUse=trackNameToUse)

    # Get the index of the original seqlets
    reverse_indices_list = {val: key for (key, listy) in enumerate(indices_list) for val in listy}
    seqlet_indices = [reverse_indices_list[idx] for idx in cluster_idx]

    return (mergedSeqlets, seqlet_indices)

# kmeans agglomerative clustering
def cluster_by_agg(distMatrix, seqletDictTask, cc_threshold = 0.92, trackNameToUse='multipliers'):
    KMEANS_NUM_CLUSTERS = 100
    KMEANS_INITS = 50
    from sklearn.cluster import MiniBatchKMeans
    clusterer = MiniBatchKMeans(n_clusters=KMEANS_NUM_CLUSTERS,
                    n_init=KMEANS_INITS)
    kmeans_labels = clusterer.fit_predict(distMatrix)

    # Merge seqlets from Kmeans clustering
    seqlet_list = csi.createMergedGrammars(kmeans_labels,
                                       seqletDictTask,
                                       subtracksToInclude=[trackNameToUse],
                                       accountForRevComp=True)

    # Prune overly similar ones
    (mergedSeqlets, indices_list) = agglomerative_clustering(grammars_list=seqlet_list.values()
                                    , gradient_track=trackNameToUse, cc_threshold=cc_threshold)

    # Get the index of the original seqlets
    reverse_indices_list = {val: key for (key, listy) in enumerate(indices_list) for val in listy}
    seqlet_indices = [reverse_indices_list[idx] for idx in kmeans_labels]


    return (mergedSeqlets, seqlet_indices)

### Cluster by tSNE
def cluster_by_tsne(distMatrix, subtracksToInclude, accountForRevComp):
    numClusters = 50    
    embedding = csi.getTsneEmbeddingOfGrammars(distMatrix, perplexity=30, verbose=2)
    labels=csi.colorTSNEembeddingBySpectralClustering(distMatrix, embedding, n_clusters=numClusters)
    mergedSeqlets = csi.createMergedGrammars(labels, distMatrix
                                  , subtracksToInclude=subtracksToInclude
                                  , accountForRevComp=accountForRevComp)
    return mergedSeqlets

# Cluster by iterative kmeans
def cluster_by_iterative_kmean(distMatrix, seqletDictTask, trackNameToUse='multipliers'):
    KMEANS_NUM_CLUSTERS_START=30
    KMEANS_NUM_CLUSTERS_END=120
    SUBKMEANS_NUM_CLUSTERS=2
    KMEANS_INITS=50
    from sklearn.cluster import MiniBatchKMeans
    clusterer = MiniBatchKMeans(n_clusters=KMEANS_NUM_CLUSTERS_START,
                    n_init=KMEANS_INITS)
    kmeans_labels = clusterer.fit_predict(distMatrix)
    # While fewer than desired cluster, split clusters (way to asses cluster purity)
    iteration=0
    while len(np.unique(kmeans_labels)) < KMEANS_NUM_CLUSTERS_END:
        print ('Subclustering iteration {0}'.format(iteration))
        new_kmeans_labels = copy.deepcopy(kmeans_labels)
        for label in np.unique(kmeans_labels):
            print(label)
            subClusterer = MiniBatchKMeans(n_clusters=SUBKMEANS_NUM_CLUSTERS,
                        n_init=KMEANS_INITS)
            subKmeans_labels = subClusterer.fit_predict(distMatrix[kmeans_labels==label])
            new_kmeans_labels[kmeans_labels==label] = [label if sublabel==0 else np.max(
                            np.unique(new_kmeans_labels))+sublabel for sublabel in subKmeans_labels]
        kmeans_labels = new_kmeans_labels
        iteration+=1

    # Merge seqlets from Kmeans clustering
    seqlet_list = csi.createMergedGrammars(kmeans_labels,
                                           seqletDictTask,
                                           subtracksToInclude=[trackNameToUse],
                                           accountForRevComp=True)

    # Prune overly similar ones
    (mergedSeqlets, indices_list) = agglomerative_clustering(grammars_list=seqlet_list.values(),
                                                             gradient_track=trackNameToUse, cc_threshold=cc_threshold)

    # Get the index of the original seqlets
    reverse_indices_list = {val: key for (key, listy) in enumerate(indices_list) for val in listy}
    seqlet_indices = [reverse_indices_list[idx] for idx in cluster_idx]

    return (mergedSeqlets, seqlet_indices)

# Cluster by DBSCAN
def cluster_by_DBSCAN(distMatrix, seqletDictTask, eps=0.3, min_samples=10, cc_threshold=0.92, trackNameToUse='multipliers'):
    from sklearn.cluster import DBSCAN
    print("Clustering using DBSCAN")
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(distMatrix)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    dbscan_labels = db.labels_

    # Merge seqlets from Kmeans clustering
    seqlet_list = csi.createMergedGrammars(dbscan_labels,
                                           seqletDictTask,
                                           subtracksToInclude=[trackNameToUse],
                                           accountForRevComp=True)
    (mergedSeqlets, indices_list) = agglomerative_clustering(grammars_list=seqlet_list.values(),
                                    gradient_track=trackNameToUse, cc_threshold=cc_threshold)

    # Get the index of the original seqlets
    reverse_indices_list = {val: key for (key, listy) in enumerate(indices_list) for val in listy}
    seqlet_indices = [reverse_indices_list[idx] for idx in dbscan_labels]

    return (mergedSeqlets, seqlet_indices)


def cluster_by_phenograph(distMatrix, seqletDictTask, min_samples=10, 
                          num_k_neighbors=30, q_tol=0.01,
                          cc_threshold=None, trackNameToUse='multipliers',
                          recursive=False, max_recursions=1):

    import phenograph
    communities, graph, Q = phenograph.cluster(distMatrix, 
                                               k=num_k_neighbors, 
                                               directed=True,
                                               n_jobs=4,
                                               min_cluster_size=2)

    if recursive:

        MIN_GROUP_SIZE=10

        # Break into sub distance matrix and append to end
        # while the number of unique communities is greater than 1
        # sub divide and append
        # remove the entry from the dictionary
        # if the number of subclusters is 1, then assign the indices to the communities 

        current_assignment = [str(el) for el in communities]

        # Store the current distance matrix and 
        comm_dist_dict = OrderedDict()
        community_index_dict = OrderedDict()
        for c in current_assignment:
            comm_ind  = np.where(communities==int(c))[0]
            comm_dist_dict[c] = distMatrix[comm_ind, :][:, comm_ind]
            community_index_dict[c] = comm_ind # Index into the length of communities


        print('Embedding in phenograph')
        #from IPython import embed; embed()

        num_recurses = 0
        # while np.max([len(el) for el in community_index_dict.values()]) > 1:
        while len(community_index_dict.values()) > 0 and num_recurses < max_recursions:

            print('Recursion: %s'%str(num_recurses))
            print('Partioning clusters: %s'%community_index_dict.keys())

            for c in community_index_dict.keys():
                if c == -1:
                    ### This means that all the -1s will probably be in one cluster (maybe re-think)
                    continue
                print('On community: %s'%c)
                sub_indices_list = np.where(np.array(current_assignment) == c)[0] # Index into the current assignment
                if len(sub_indices_list) <= MIN_GROUP_SIZE:
                    del comm_dist_dict[c]
                    del community_index_dict[c]
                    continue
                # sub_current_assignment = [community_index_dict[c][i2] for i2 in sub_indices_list] # Index into the length of communities
                sub_corr_mat = distMatrix[sub_indices_list, :][:, sub_indices_list]
                sub_communities, sub_graph, sub_Q = phenograph.cluster(sub_corr_mat, 
                                                                       k=np.min([num_k_neighbors, len(sub_indices_list)-4]),  # -1 and -2 can throw bugs
                                                                       directed=True,
                                                                       n_jobs=4,
                                                                       min_cluster_size=2)
                print('Broke cluster %s of size %s into %s sub-clusters'%(c, len(sub_communities), np.max(np.unique(sub_communities))+1))
                if np.max(np.unique(sub_communities)) == 0: 
                    print("Stopping for cluster %s, no further sub-clusters"%c)
                    # If there is only one sub-clustered community then keep current assignments and delete from dictionary
                    del comm_dist_dict[c]
                    del community_index_dict[c]
                    continue

                for sub_c in np.unique(sub_communities):
                    if sub_c == -1:
                        ### This means that all the -1s will probably be in one cluster (maybe re-think)
                        continue
                    new_name = str(c) + '_' + str(sub_c)
                    sub_comm_ind  = np.where(np.array(sub_communities) == sub_c)[0]
                    comm_dist_dict[new_name] = sub_corr_mat[sub_comm_ind, :][:, sub_comm_ind] 
                    community_index_dict[new_name] = [community_index_dict[c][i] for i in sub_comm_ind]
                    for i in community_index_dict[new_name]:
                        current_assignment[i] = new_name

                # current_assignment = []

                # Delete the larger clusters
                del comm_dist_dict[c]
                del community_index_dict[c]

            # Add to the number of recursions
            num_recurses += 1


        # Assign final communities
        # assign_dict = {c_assign: num for (num, c_assign) in np.unique(current_assignment)}
        # communities = [assign_dict[ca] for ca in current_assignment]
        communities = current_assignment
        # Possible convert to NUMBERS AGAIN


    seqlet_list = csi.createMergedGrammars(communities,
                                           seqletDictTask,
                                           subtracksToInclude=[trackNameToUse],
                                           accountForRevComp=True)
    # Get rid of -1, not assigned seqlets
    if -1 in seqlet_list.keys():
        seqlet_list.pop(-1)
    # Merge
    if cc_threshold is not None:
        (mergedSeqlets, indices_list) = agglomerative_clustering(grammars_list=seqlet_list.values(),
                                        gradient_track=trackNameToUse, cc_threshold=cc_threshold)
        # Get the index of the original seqlets
        reverse_indices_list = {val: key for (key, listy) in enumerate(indices_list) for val in listy}
        seqlet_indices = [reverse_indices_list[idx] for idx in communities]
        return (mergedSeqlets, seqlet_indices)

    # seqlet list is a dictionary, communities is a list
    return (seqlet_list.values(), communities)


    # import deeplift
    # from deeplift.visualization import viz_sequence
    # for key in seqlet_list.keys():
    #     f = viz_sequence.plot_weights(seqlet_list[key].summedCoreDeepLIFTtrack)
    #     f.savefig('/users/pgreens/test_grammar_%s.png'%key)




