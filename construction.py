import networkx as nx
import numpy as np
import pandas as pd


def load_data(path,sep=" "):
    """load the cvs file representing the temporal graph

    Parameters:
    path (string): path of the input dataset

    Returns:
    np.array: a np array with the loaded data

    """
    data = pd.read_csv(path,sep,names=["t","a","b"])
    return(data)


def individuals(data):
    res = []
    res.extend(np.unique(data.a))
    res.extend(np.unique(data.b))
    return(np.unique(res))

def build_graphs(data,gap=19,with_labels=False,meta_path=None):
    graphs = []
    G=nx.Graph()
    nodes = individuals(data)
    G.add_nodes_from(nodes)
    
    splitted_data = split_input_data(data,gap)
    for t in splitted_data:
        g = G.copy()
        for _,i,j in t:
            if not i == j:
                g.add_edge(i,j)
        graphs.append(g)  
    return(graphs)

def split_input_data(data, gap=19):
    times = [int(x/(gap+1)) for x in data.t]
    data.t = times
    splitted_data = []
    c = 0 

    for i in range(max(times)+1):
        tmp = data[data.t == i].to_numpy()
        if tmp.shape[0] == 0:
            tmp = [[i,0,0]]

        splitted_data.append(tmp)
    return splitted_data
