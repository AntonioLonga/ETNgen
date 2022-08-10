import construction as cs
from ETN import *

# given:
# "0b001001001" return ["001","001","001"]
def split_etns(etns,k):
    splitted_etns = []
    for i in range(0,len(etns),k):
        splitted_etns.append(etns[i:i+k])

    return(splitted_etns)

# given
# ETNS, return
# dict 10x _--> 101 100001001001 etc
def get_dict(ETNS,k,return_statistics=False):
    '''It generates the dictionary. It takes as input k and ETNS, which is a
    dictionary (referring for instance to one hour of the day, or better to
    local_split) with keys the signatures of the original graph and for values
    their nb of occurrences. The dictionary that it creates puts together the
    ETNSs that have the same key. Then it transforms frequencies into
    probabilities.
    New dictionary example:
    key: 11x
    value: [[110,111],[0.7,0.3]]
    '''
    # new_node == 001,01,0001,00001
    new_node = "0"*k+"1"
    ETNS_list = list(ETNS.keys())
    diz = dict()
    for etns in ETNS_list:
        etns = etns[2:]
        if len(etns) % (k+1) == 0:
            splitted_etns = split_etns(etns,k+1)
            key = ""
            for etn in splitted_etns:
                if not etn == new_node: # per il caso 001 non serve creare la key, che è semplicemente ""
                    key = key + etn[0:k]+"x"
            if key in diz: # append to diz the new pruned signature with its nb. of occurrences
                diz[key][0].append("0b"+etns)
                diz[key][1].append(ETNS["0b"+etns])
            else:
                diz[key] = [["0b"+etns],[ETNS["0b"+etns]]]
    statistics = 0
    for key,value in diz.items():
        summ = sum(diz[key][1]) # from frequences to probabilities
        statistics += summ
        c = 0
        for val in value[1]:
            diz[key][1][c] = diz[key][1][c]/summ
            c = c + 1

    if return_statistics:
        return diz,statistics
    else:
        return diz




# from etns to key
# input: 101-011 --> 01x-11x
def create_key(etns,k):
    if len(etns) == 0:
        return ""
    else:
        return "x".join(split_etns(etns,k))+"x"



# counts etns only for a given node
def count_ETN_given_node(graphs,k,node):
    '''It returns all the etns (at each timestamp) of a given node'''
    v = node
    etn = build_ETN(graphs[:k+1],v) # etn is a nx.graph representing the motif with node v as ego
    if not etn == None:
        etns,node_encoding = get_ETNS_with_encoding(etn) # etns is the corresponding signature and node_encoding tells the node identity

    return (etns[2:],node_encoding)


# get ETNS and node generating etns
def get_ETNS_with_encoding(ETN):
    '''It translates the ETN (a nx.graph) in signature. It return the signature etns and
    node_encoding, which tells at which node the signature is referring.
    For instance:
    ETNS,node_encoding = 1011 {'65': '10', '56': '11'}
    means that node n has had an interaction 10 with node 65 and
    an interaction 11 with node 56.
    '''
    nodes = list(ETN.nodes())

    nodes_no_ego = []
    ids_no_ego = []
    lenght_ETNS = 0
    for n in nodes:
        if not ("*" in n):
            nodes_no_ego.append(n)
            if not(n.split("_")[0] in ids_no_ego):
                ids_no_ego.append(n.split("_")[0])
        else:
            ego = int(n.split("*")[0])
            lenght_ETNS = lenght_ETNS + 1

    node_encoding = get_node_encoding(ids_no_ego,nodes_no_ego,lenght_ETNS)
    for k in node_encoding.keys():
        node_encoding[k] = ''.join(str(e) for e in node_encoding[k])


    binary_node_encodings = list(node_encoding.values())
    binary_node_encodings.sort()

    ETNS = '0b'+''.join(e for e in binary_node_encodings)

    return(ETNS,node_encoding)




import random


# givend diz and key:
# return a key according to the probablity.
def get_random_etns(diz,key,k):
    '''It returns a new etns (like 101), given a key (like 10x), according to the probability stored in the dictionary.'''
    if (key in diz):
        etns = diz[key][0]
        prob = diz[key][1]
        cumuative = [sum(prob[0:i+1]) for i in range(len(prob))]
        cumuative = [0]+cumuative
        r = random.random()
        for i in range(len(cumuative)-1):
            if r >= cumuative[i]  and r < cumuative[i+1]:
                return etns[i]
    elif (key[:-(k+1)] in diz): # approximate key if the original key is not found
        key = key[:-(k+1)]
        etns = diz[key][0]
        prob = diz[key][1]
        cumuative = [sum(prob[0:i+1]) for i in range(len(prob))]
        cumuative = [0]+cumuative
        r = random.random()
        for i in range(len(cumuative)-1):
            if r >= cumuative[i]  and r < cumuative[i+1]:
                return etns[i]
    else:
        return None





# create edge_list_g2:
# given etns2 and etns3 merge them to create edge_list in g2
def create_edge_g2(n,etns3,node_encoding,k):
    '''
    It returns the list of edges to add at the (k+1)-th layer, given node n and the
    egocentric neighborhood etns3.
    Example:
    n = 36
    etns3 = '0b100100100111'
    node_enc = {'99': '10', '56': '10', '65': '11', '32': '10'}
    Returns [(36, 65)]
    '''
    new_node = k*"0"+"1"
    edges = []
    for split in split_etns(etns3[2:],k+1):
        if split == new_node:
            edges.append((n,"x"))
        elif split[-1] == "1":
            # search node correspondensy in node_embedding and remove:
            for key,value in node_encoding.items():
                if value == split[:k]:
                    edges.append((n,int(key)))
                    node_encoding[key] = "USED"

    return edges


# split ("x","a") and ("a","b")
def split_stub(e):
    '''Given a list of edges, it discerns between directed edges and stubs'''
    edges = []
    stubs = []

    for i,j in e:
        if i == "x" or j == "x":
            stubs.append((i,j))
        else:
            edges.append((i,j))

    return edges,stubs


########### originale
'''
def get_edges_to_keep(edges,alpha=0.5):

    nb_edges = int(len(edges)*alpha)
    edges_to_keep = []
    # check double directions
    for i,j in edges:
        if (j,i) in edges:
            if (i,j) not in edges_to_keep and (j,i) not in edges_to_keep:
                edges_to_keep.append((i,j))

    old = nb_edges
    nb_edges = nb_edges - len(edges_to_keep)
    print(old,nb_edges)
    # remove taken edges:
    for i,j in edges_to_keep:
        edges.remove((i,j))
        edges.remove((j,i))

    np.random.shuffle(edges)

    edges_to_keep.extend(edges[0:nb_edges])
    # take only 0.5 edges
    np.random.shuffle(edges_to_keep)
    edges_to_keep = edges_to_keep[0:nb_edges]

    return edges_to_keep
'''

################ TEST
def get_edges_to_keep(edges,alpha=0.5):
    '''Given the list of directed edges, it returns only the bidirectional ones
    + a fraction alpha of the others.
    '''

    edges_to_keep = []
    #print('Initial nb of edges:', len(edges))
    for i,j in edges:
        if (j,i) in edges:
            if (i,j) not in edges_to_keep and (j,i) not in edges_to_keep:
                edges_to_keep.append((i,j))
                edges.remove((i,j))
                edges.remove((j,i))

    nb_edges = len(edges)
    #print('Bidirectional edges to keep:',len(edges_to_keep))
    #print('Nb of remaining edges:', len(edges))
    nb_edges = int(len(edges)*alpha) # si approssima per difetto, che succede se si approssima per eccesso?
    #print('Nb of remaining edges we want to keep:', nb_edges)

    np.random.shuffle(edges)
    edges_to_keep.extend(edges[0:nb_edges])
    #print('Total nb of edges to keep:',len(edges_to_keep))

    return edges_to_keep


# merge stubs
# take in input [(91, 'x'), (42, 'x'), (60, 'x'), (78, 'x'), (77, 'x'), (21, 'x')]
# return (91,42),(60,21) etc
def merge_stubs(stubs):
    np.random.shuffle(stubs)
    edges = []

    node_seq = []
    for i in range(len(stubs)):
        u = stubs[i][0]
        node_seq.append(u)


    flag = True
    while flag:
        if len(node_seq) >= 2:
            u = np.random.choice(node_seq)
            node_seq.remove(u)
            v = np.random.choice(node_seq)
            node_seq.remove(v)
            if not u == v:
                edges.append((u,v))

            else:
                node_seq.append(u)
                node_seq.append(v)
                if len(np.unique(node_seq))==1:
                    flag = False
        else:
            flag = False

    return edges


# merge edges and stubs
def get_edges_g2(edges,alpha):
    edges,stubs = split_stub(edges)
    edges_g2 = get_edges_to_keep(edges,alpha)
    edges_g2.extend(merge_stubs(stubs))

    return edges_g2




# given node ed edges build graph g3
def build_graph_g2(edges,nodes):

    g2 = nx.Graph()
    g2.add_nodes_from(nodes)
    g2.add_edges_from(edges)

    return g2



# given g0,g1,diz and nodes retrung a nx graph
def generate_graph_g2(nodes,graphs,diz,k,alpha):
    '''It generates a new temporal layer, given the previous k, the dictionary and the nodes list'''
    edges = []
    for n in graphs[0].nodes(): #build provisional layer
        etns2,node_enc = count_ETN_given_node(graphs,k,n) # etns of the last k layers for node n
        key = create_key(etns2,k)
        etns3 = get_random_etns(diz,key,k) #i.e. etns3 = '0b100'
        if not etns3 == None:
            edges.extend(create_edge_g2(n,etns3,node_enc,k))

    edges_g2 = get_edges_g2(edges,alpha)
    g2 =build_graph_g2(edges_g2,nodes)
    return g2



#given seed and keys generate temporal graph
def generate_temporal_graph(nb_graphs, graph_seed,diz,k,alpha):
    nodes = list(graph_seed[0].nodes())
    tg = graph_seed
    for i in range(nb_graphs-1):
        graphs_in = tg[i:i+k]
        g_new = generate_graph_g2(nodes,graphs_in,diz,k,alpha)
        tg.append(g_new)

    return tg





# generate a single graph_seed given the previous graphs
def generate_seed_graph(graphs,graphs_seed,k,alpha):
    nodes = list(graphs_seed[0].nodes())
    ETNS = count_ETN(graphs,k)
    ETNS = {k: v for k, v in sorted(ETNS.items(),reverse=True, key=lambda item: item[1])}
    ETNS_list = list(ETNS.keys())
    diz = get_dict(ETNS,k)

    new_g = generate_graph_g2(nodes,graphs_seed,diz,k,alpha)
    graphs_seed.append(new_g)
    return graphs_seed


# generate k graph seed given g0
def generate_seed_graphs(g0,graphs,k,alpha):
    graphs_seed = [g0]

    for i in range(k-1):
        graphs_seed = generate_seed_graph(graphs,graphs_seed,i+1,alpha)

    return graphs_seed
