from platform import node
import construction as cs
import networkx as nx
import numpy as np

def conn_nodes(graphs):
    conn_nodes = set()
    for g in graphs[0:100]:

        for a,b in list(g.edges()):
            conn_nodes.add(a)
            conn_nodes.add(b)

    return list(conn_nodes)
    
def coverage(graphs,K,T):
    nodes = conn_nodes(graphs)
    if nodes == []:
        nodes= [0]
    res = []
    for i in range(K):

        t0 = np.random.randint(0,len(graphs) - T -1)
        n0 = np.random.choice(nodes)

        rw = random_walk(graphs,t0,n0,T)
        res.append(len(np.unique(rw)))
    return res


def random_walk(graphs,t0,n0,T=100):
    RW = []
    for i in range(T):
        RW.append(n0)
        t0 = t0 + 1 
        g0 = graphs[t0]
        neig = list(nx.neighbors(g0,n0))
        if not neig == []:
            n0 = np.random.choice(neig)
            
    return RW


def len_random_walk_MFPT(graphs,t0,n0,n_end):
    RW = []
    for i in range(t0,len(graphs)-1):
        RW.append(n0)
        t0 = t0 + 1
        g0 = graphs[t0]
        neig = list(nx.neighbors(g0,n0))
        if not neig == []:
            n0 = np.random.choice(neig)
            if n0 == n_end:
                return len(RW)
            
    return -1



def MFPT(graphs,K=5):
    nodes = list(graphs[0].nodes())
    res = []
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            t0 = 0
            n0 = nodes[i]
            n_end = nodes[j]
            tmp = []
            for k in range(K):
                t = len_random_walk_MFPT(graphs,t0,n0,n_end)
                if t == -1:
                    tmp.append(len(graphs))
                else:
                    tmp.append(t)

            res.extend(tmp)
    return res






def SIR_model(graphs,n0,lambd = 0.5,mu = 0.01):
    
    Is = []
    I = set()
    R = set()
    S = set(graphs[0].nodes())

    S.remove(n0)
    I.add(n0)
    r0 = set()
    c = 0 
    for g in graphs:
        c = c + 1 
        add_and_remove = set()
        for i in I:    
            neigs = list(g.neighbors(i))
            for n in neigs:
                if not n in R and not n in I:
                    r = np.random.rand()
                    if lambd >= r: # became infected
                        add_and_remove.add(n)
                        
            if i == n0:
                for jj in add_and_remove:
                    r0.add(jj)
                
            
        S = S.difference(add_and_remove)
        I = I.union(add_and_remove)
        
        
        add_and_remove = set()
        for i in I:
            r = np.random.rand()
            if mu >= r:
                add_and_remove.add(i)

        I = I.difference(add_and_remove)
        R = R.union(add_and_remove)

        Is.append(list(I))
        #test
        assert(I.intersection(R) == set())
        assert(I.intersection(S) == set())
        assert(S.intersection(R) == set())
        
    return I,S,R,Is,len(r0)



def connected_initial_nodes(graphs):
    
    res = []
    c = 0
    while res == []:

        g = graphs[c]

        res = []
        for i,j in g.edges():
            if i not in res:
                res.append(i)
            if j not in res:
                res.append(j)
        c = c + 1   
    return res

def compute_r0(K,graphs,lambd,mu):
    res = []
    nodes = connected_initial_nodes(graphs)
    for i in range(K):
        
        n0 = np.random.choice(nodes)
        _,_,_,_,r0 = SIR_model(graphs,n0,lambd = lambd,mu = mu)

        res.append(r0)
    return res