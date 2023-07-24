import numpy as np
from ETN import *
import pickle

def store_etm_counts(ETM_counts,file_name,gap,k,label):
    if label:
        name="gap_"+str(gap)+"_k_"+str(k)+"_LABEL.json"
    else:
        name="gap_"+str(gap)+"_k_"+str(k)+".json"

    directory = "res/"+file_name+"/ETM_counts/"
    if not os.path.exists(directory):
        os.makedirs(directory)

    a_file = open(directory+name, "w")
    json.dump(ETM_counts, a_file,indent=1)
    a_file.close()

    print("file stored in: \t"+directory+name)


def load_etm_count(file_name,gap,k,label):
    directory = "res/"+file_name+"/ETM_counts/"
    if label:
        name="gap_"+str(gap)+"_k_"+str(k)+"_LABEL.json"
    else:
        name="gap_"+str(gap)+"_k_"+str(k)+".json"

    with open(directory+name) as json_file:
        ETM_counts = json.load(json_file)       

    return ETM_counts


def get_ETM(counts,alpha,beta,gamma):
    over = over_representation(counts,alpha)
    mdev = minimum_deviation(counts,beta)
    mfrq = minimum_frequency(counts,gamma)
    valid = mfrq * over * mdev

    ETM = []
    ETNS = list(counts.keys())
    for k in range(len(ETNS)):
        if(valid[k]==1):
            ETM.append([ETNS[k],counts[ETNS[k]][0]])

    print("number of etns:\t",len(ETNS),"\nnumber of etm: \t",len(ETM))
    return(ETM)



def minimum_deviation(counts,beta):
    N_G = np.array(list(counts.values()))[:,0]
    N_G0 = np.array(list(counts.values()))[:,1:]

    N_G0_mean = np.mean(N_G0,-1)
    valid = []
    for i in range(len(N_G)):

        if(N_G[i] - N_G0_mean[i] > beta * N_G0_mean[i]):
            valid.append(1)
        else:
            valid.append(0)
            
    return(np.array(valid))
    
    
def over_representation(counts,alpha=0.01):
    N_G = np.array(list(counts.values()))[:,0]
    N_G0 = np.array(list(counts.values()))[:,1:]
    alpha = 0.01
    valid = []
    for j in range(len(N_G)):
        P = 0
        for i in N_G0[j]:
            if (i > N_G[j]):
                P = P + 1
        P = P/len(N_G0[j])
        if(P<alpha):
            valid.append(1)
        else:
            valid.append(0)
    return(np.array(valid))

def minimum_frequency(counts,gamma=5):
    
    N_G = np.array(list(counts.values()))[:,0]
    valid = []
    for i in N_G:
        if(i >= gamma):
            valid.append(1)
        else:
            valid.append(0)
    return(np.array(valid))



def counts_ETN_null_models(null_models,S,k,label,meta=None,verbose=False):
    if label:
        if meta == None:
            print("error meta is none and label is true")
    counts = dict()
    for i in list(S.keys()):
        counts[i] = [S[i]]
    c = 1
    for null_model in null_models:
        S_i = count_ETN_null_model(S,null_model,k,label,meta=meta)
        for j in list(S_i.keys()):
            counts[j].append(S_i[j])
        if (verbose):
            print("done",c)
        c = c + 1
    return (counts)



def count_ETN_null_model(S_in,graphs,k,label=False,meta=None):
    S = S_in.copy()
    for i in S:
        S[i] = 0
    for i in range(len(graphs)-k + 1):
        for v in graphs[i].nodes():

            etn = build_ETN(graphs[i:i+k+1],v)
            if not etn == None:
                if label:
                    etns = get_ETNS(etn,meta)
                else:
                    etns = get_ETNS(etn)
                if etns in S.keys():
                    S[etns] = S[etns] + 1
    return(S)


def shuffle_graphs(graphs_in,n,seed): # Shuffle array of graphs
    graphs = graphs_in.copy()
    null_models = []
    for i in range(n):
        np.random.seed(seed+i)
        np.random.shuffle(graphs)
        null_models.append(graphs.copy())
        
    return(null_models)