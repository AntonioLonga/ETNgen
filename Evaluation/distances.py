
import numpy as np
from numpy import linalg as la
import networkx as nx
from scipy import sparse as sps
from scipy.sparse import issparse
from scipy import stats



_eps = 10**(-10) # a small parameter
def _canberra_dist(v1,v2):
    eps = 10**(-15)
    v1,v2 = [_flat(v) for v in [v1,v2]]
    d_can = 0
    for u,w in zip(v1,v2):
        if np.abs(u)<eps and np.abs(w)<eps:
            d_update = 1
        else:
            d_update = np.abs(u-w) / (np.abs(u)+np.abs(w))
        d_can += d_update
    return d_can


def aggregate_features(feature_mat,row_var=False,as_matrix=False):
    axis = int(row_var) # 0 if column-oriented, 1 if not
    description = np.array([feature_mat.mean(axis=axis),
                            np.median(feature_mat,axis=axis),
                            np.std(feature_mat,axis=axis),
                            stats.skew(feature_mat,axis=axis),
                            stats.kurtosis(feature_mat,axis=axis)])
    if not as_matrix:
        description = description.flatten()
    return description

def _flat(D):
    if issparse(D):
        raise ValueError('Cannot flatten sparse matrix.')
    d_flat = np.array(D).flatten()
    return d_flat


def get_features_temporal(A,As):
    try:
        G = nx.from_scipy_sparse_matrix(A)
    except AttributeError:
        G = nx.from_numpy_matrix(A)
    n = len(G)
    # degrees, array so we can slice nice
    d_vec = np.array(list(G.degree().values()))
    # list of clustering coefficient
    clust_vec = np.array(list(nx.clustering(G).values()))
    neighbors = [G.neighbors(i) for i in range(n)]
    # average degree of neighbors (0 if node is isolated)
    neighbor_deg = [d_vec[neighbors[i]].sum()/d_vec[i]
                    if d_vec[i]>_eps else 0 for i in range(n)]
    # avg. clustering coefficient of neighbors (0 if node is isolated)
    neighbor_clust = [clust_vec[neighbors[i]].sum()/d_vec[i] 
                    if d_vec[i]>_eps else 0 for i in range(n)]
    egonets = [nx.ego_graph(G,i) for i in range(n)]
    # number of edges in egonet
    ego_size = [G.number_of_edges() for G in egonets]
    # number of neighbors of egonet
    ego_neighbors = [len(set.union(*[set(neighbors[j])
                                     for j in egonets[i].nodes()]) -
                         set(egonets[i].nodes()))
                     for i in range(n)]
    # number of edges outgoing from egonet
    outgoing_edges = [len([edge for edge in G.edges(egonets[i].nodes()) 
                           if edge[1] not in egonets[i].nodes()]) 
                      for i in range(n)]

    # Codice mio
    # !!!!!!!!!!!!!!!!!!!!!!!!!!
    v = get_As_features(As) 
    # use mat.T so that each node is a row (standard format)
    feature_mat = np.array([d_vec,clust_vec,neighbor_deg,neighbor_clust,
                            ego_size,ego_neighbors,outgoing_edges,v]).T


    return feature_mat


def get_As_features(As):
    n = As[0].shape[0]
    interactions = dict()
    for i in np.arange(n):
        interactions[i] = 0 

    for A in As:
        try:
            G = nx.from_scipy_sparse_matrix(A)
        except AttributeError:
            G = nx.from_numpy_matrix(A)

        for i,j in G.edges():
            interactions[i] = interactions[i] + 1
            interactions[j] = interactions[j] + 1
    return(list(interactions.values()))


def netsimile2(G1,G2,graphs1,graphs2):

    As1 = [nx.adjacency_matrix(g) for g in graphs1]
    A1 = nx.adjacency_matrix(G1)
    As2 = [nx.adjacency_matrix(g) for g in graphs2]
    A2 = nx.adjacency_matrix(G2)
    feat_A1 = get_features_temporal(A1,As1)
    feat_A2 = get_features_temporal(A2,As2)

    agg_A1,agg_A2 = [aggregate_features(feat) for feat in [feat_A1,feat_A2]]
    # calculate Canberra distance between two aggregate vectors
    d_can = _canberra_dist(agg_A1,agg_A2)
    return d_can











# weighted laplacian
from scipy import sparse as sps
import numpy as np
from scipy.sparse import linalg as spla
from numpy import linalg as la
from scipy.sparse import issparse


def _eigs(M,which='SR',k=None):
    n,_ = M.shape
    if k is None:
        k = n
    if which not in ['LR','SR']:
        raise ValueError("which must be either 'LR' or 'SR'.")
    M = M.astype(float)
    if issparse(M) and k < n-1:
        evals,evecs = spla.eigs(M,k=k,which=which)
    else:
        try: M = M.todense()
        except: pass
        evals,evecs = la.eig(M)
        # sort dem eigenvalues
        inds = np.argsort(evals)
        if which == 'LR':
            inds = inds[::-1]
        else: pass
        inds = inds[:k]
        evals = evals[inds]
        evecs = np.matrix(evecs[:,inds])
    return np.real(evals),np.real(evecs)



def weighted_laplacian(G1,graphs1):
    As1 = [nx.adjacency_matrix(g) for g in graphs1]
    W = sum(As1)

    A1 = nx.adjacency_matrix(G1)
    n = len(A1.A[0])
    D = np.zeros([n,n])
    degs = np.sum(A1.A,axis=1)
    for i in range(len(degs)):
        D[i,i] = degs[i]

    L = D - W
    
    return(L)



def lambda_dist_weighted_L(G1,G2,graphs1,graphs2,k=None,p=2):
    # norma alla p, usando i primi K eigs
    
    L1 = weighted_laplacian(G1,graphs1)
    L2 = weighted_laplacian(G2,graphs2)
    
    if k == None:
        k = min(L1.shape[0],L2.shape[0])
        
    # get eigenvalues, ignore eigenvectors
    evals1,evals2 = [_eigs(L)[0] for L in [L1,L2]]
    dist = la.norm(evals1[:k]-evals2[:k],ord=p)
    
    return dist
    



# etmm distances
from ETN import *
from ETMM import *
import math
from scipy.stats import pearsonr

def compute_extra_correlation(structures,names):
    etns = []
    for i in structures:
        etns.extend(np.array(i)[:,0])
    etns = np.unique(etns)
    etns = sorted(etns, key=lambda x: int(x, 2))

    dict_etns = dict()
    for i in etns:
        dict_etns[i] = np.zeros(len(structures))

    for i in range(len(structures)):
        for string,count in structures[i]:
            dict_etns[string][i]=count

    counts = []
    for i in range(len(structures)):
        tmp = np.array(list(dict_etns.values()))[:,i]
        counts.append(tmp/np.max(tmp))

        
    return counts




#### extra correlation
def load_structures(names,structure_type,gap,k,label,alpha,beta,gamma):

    structures = []
    for name in names:

        if structure_type == "ETM":
            c = load_etm_count(name,gap,k,label)
            structures.append(get_ETM(c,alpha=alpha,beta=beta,gamma=gamma))

        elif structure_type == "ETN":
            c = load_etm_count(name,gap,k,label)
            ETN = []
            for i,j in c.items():
                ETN.append([i,j[0]])
            structures.append(ETN)
        else:
            print("specify ETM or ETN")

    return structures

def etmm_distance(names,structure_type,gap,k,label,alpha=0.01,beta=0.1,gamma=50):

    structures = load_structures(names,structure_type,gap,k,label,alpha,beta,gamma)
    distribution_1,distribution_2 = compute_extra_correlation(structures,names)
    dist = from_corr_to_distance(distribution_1,distribution_2)
    
    return(dist)
    
    
def from_corr_to_distance(a,b):
    v = np.abs(pearsonr(a,b)[0])
    
    return math.sqrt(1-v)
