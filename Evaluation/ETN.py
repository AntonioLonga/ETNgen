from copy import deepcopy
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import itertools
import json
import os

def store_etns(S,file_name,gap,k,label):

    if label:
        name="gap_"+str(gap)+"_k_"+str(k)+"_LABEL.json"
    else:
        name="gap_"+str(gap)+"_k_"+str(k)+".json"

    directory = "res/"+file_name+"/ETNS/"
    if not os.path.exists(directory):
        os.makedirs(directory)

    

    a_file = open(directory+name, "w")
    json.dump(S, a_file,indent=1)
    a_file.close()

    print("file stored in: \t"+directory+name)


def load_etns(file_name,gap,k,label):
    directory = "res/"+file_name+"/ETNS/"
    if label:
        name="gap_"+str(gap)+"_k_"+str(k)+"_LABEL.json"
    else:
        name="gap_"+str(gap)+"_k_"+str(k)+".json"

    with open(directory+name) as json_file:
        S = json.load(json_file)       

    return S

def count_ETN(graphs,k,meta=None):
    S = dict()
    for i in range(len(graphs)-k + 1):
        for v in graphs[i].nodes():
            etn = build_ETN(graphs[i:i+k+1],v)
            if not etn == None:
                etns = get_ETNS(etn,meta)
                if etns in S.keys():
                    S[etns] = S[etns] + 1
                else:
                    S[etns] = 1
    return(S)





def get_node_encoding_labeled(meta,node_encoding,ego):
    categories = np.sort(list(np.unique(list(meta.values())))+["0"])

    meta_binary = list(itertools.product([0, 1], repeat=round(len(categories)**(1/2)+0.5)))
    meta_dict = dict()
    for i in range(len(categories)):    
        meta_dict[categories[i]] = list(meta_binary[i])


    new_node_encoding = dict()
    for i in node_encoding:
        tmp = []
        for v in node_encoding[i]:
            if(v==0):
                tmp.extend(meta_dict["0"])
            else:
                tmp.extend(meta_dict[meta[int(i)]])

        new_node_encoding[i]=tmp

    ego_encoding = meta_dict[meta[ego]]

    return(new_node_encoding,ego_encoding)
        
    


def from_ETNS_to_ETN(s,k,meta=None):

    n = k + 1
    if not(meta==None):
        categories = np.sort(list(np.unique(list(meta.values())))+["0"])
        meta_binary = list(itertools.product([0, 1], repeat=round(len(categories)**(1/2)+0.5)))
        meta_dict = dict()
        for i in range(len(categories)):
            value = "".join(str(e) for e in meta_binary[i])
            meta_dict[value] = categories[i]

        s = s[2:] # update s
        ego_encoding = s[:len(meta_binary[0])]
        s = s[len(meta_binary[0]):] # update s 

        node_encoding = [s[i:i+len(meta_binary[0])] for i in range(0,len(s),len(meta_binary[0]))]


        new_node_encoding = []
        for i in node_encoding:
            tmp = []
            for j in range(0, len(i), len(meta_binary[0])):
                value = i[j:j+len(meta_binary[0])]
                tmp.append(meta_dict[value])

            new_node_encoding.append(tmp)
        node_encoding_labels = new_node_encoding
    else:
        node_encoding = [s[2:][i:i+1] for i in range(0, len(s[2:]))]

    egos = [("0*_"+str(i),"0*_"+str(i+1)) for i in range(n-1)]



    ETN = nx.Graph()
    ETN.add_edges_from(egos)
    # add ego labels 
    if not(meta == None):
        for nod in list(ETN.nodes()):
            ETN.nodes()[nod]["label"] = meta_dict[ego_encoding]


    new_node_encoding = []

    for j in node_encoding:
        if("1" in j):
            new_node_encoding.append(1)
        else:
            new_node_encoding.append(0)

    node_encoding = new_node_encoding

    for i in range(0,len(node_encoding),n):
        same_person = []
        for j in range(n):
            if not(node_encoding[i+j] == 0):
                ETN.add_edge("0*_"+(str(j)),str(i+1)+"_"+str(j))
                same_person.append(str(i+1)+"_"+str(j))
                if not (meta == None):
                    ETN.nodes()[str(i+1)+"_"+str(j)]["label"] = node_encoding_labels[i+j][0]

        if len(same_person) > 1:
            for k in range(len(same_person)-1):
                ETN.add_edge(same_person[k],same_person[k+1])



    return(ETN)

       
    
def get_ETNS(ETN,meta=None):
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
    if not(meta == None):
        node_encoding,ego_encoding = get_node_encoding_labeled(meta,node_encoding,ego)
    
    for k in node_encoding.keys():
        node_encoding[k] = '0b'+''.join(str(e) for e in node_encoding[k])

    binary_node_encodings = list(node_encoding.values())
    binary_node_encodings.sort()
    
    ETNS = '0b'+''.join(e[2:] for e in binary_node_encodings)
    
    # add ego label encoding
    if not (meta == None):
        ETNS  = '0b'+''.join(str(e) for e in ego_encoding)+ETNS[2:]
    
    return(ETNS)



def get_node_encoding(ids_no_ego,nodes_no_ego,lenght_ETNS):

    node_encoding = dict()
    for n in ids_no_ego:
        enc = []
        for k in range(lenght_ETNS):
            if (str(n)+"_"+str(k) in nodes_no_ego):
                enc.append(1)
            else:
                enc.append(0)
        node_encoding[n]=enc
        
    return(node_encoding)




def get_egocentric_neighborhood(g,v):
    return([str(v)+"*"]+list(g.neighbors(v)))

def build_ETN(graphs,v):
    if len(list(graphs[0].neighbors(v))) > 0 :
        en_list = []
        en_ids = []
        for i in graphs:
            en_list.append(get_egocentric_neighborhood(i,v))
            en_ids.append(get_egocentric_neighborhood(i,v))

        # add temporal step to en_list
        for i in range(len(en_list)):
            for j in range(len(en_list[i])):
                en_list[i][j] = str(en_list[i][j])+"_"+str(i)


        #buld graphs en
        en_graph = []
        for en in en_list:
            g = nx.Graph()
            for i in range(len(en)-1):
                g.add_edge(en[0],en[i+1])
            en_graph.append(g)

        # compose an to get disconnected ETN 
        ETN = nx.Graph()
        for g in en_graph:
            ETN = nx.compose(ETN,g)


        # merge en
        en_list_long = []
        for en in en_list:
            en_list_long_tmp = []
            for n in en:
                en_list_long_tmp.append(n.split("_")[0])

            en_list_long.append(en_list_long_tmp)

        for k in range(len(en_list_long)-1):
            for n in en_list_long[k]:
                for en in range(len(en_list_long[k+1:])):
                    add = False
                    if (n in en_list_long[k+en+1]):
                        add = True 
                        t = k + 1 + en
                        break
                if (add == True):
                    u = str(n)+"_"+str(k)
                    v = str(n)+"_"+str(t)
                    ETN.add_edge(u,v)
        return(ETN)
    else:
        return(None)



def draw_ETN(ETN,multiple=False):
    ids,k = get_ids_and_k(ETN)
    pos = dict()
    id_ego = []
    for t in range(k+1):
        for i in ids:
            if "*" in i:
                id_ego.append(str(i)+"_"+str(t))
                pos[str(i)+"_"+str(t)] = [t,int(i[0])]
            else:
                pos[str(i)+"_"+str(t)] = [t,int(i)]
                
    node_label = dict() 
    nodes_data = dict(ETN.nodes(data=True))
    for i in list(ETN.nodes()):
        if not(nodes_data[i] == {}):
            node_label[i] = nodes_data[i]["label"]
            
            
    if (node_label == {}):
        nx.draw(ETN,pos=pos,node_size=100,alpha=0.9,with_labels=True)
        nx.draw_networkx_nodes(ETN, pos, nodelist=id_ego, node_size=300, node_color='red',alpha=0.5)
    else:
        nx.draw(ETN,pos=pos,node_size=100,alpha=0.5)
        nx.draw_networkx_nodes(ETN, pos, nodelist=id_ego, node_size=300, node_color='red',alpha=0.5)
        nx.draw_networkx_labels(ETN, pos, labels=node_label, font_size=12)

    if not multiple:
        plt.show()


# get unique ids (no time consideration)
def get_ids_and_k(ETN):
    nodes = list(ETN.nodes())
    ids = []
    k = 0
    for n in nodes:
        id_n = (n.split("_")[0])
        k_tmp = (n.split("_")[1])
        if not (id_n in ids):
            ids.append(id_n)
        if int(k_tmp) > k:
            k = int(k_tmp)
    return(ids,k)