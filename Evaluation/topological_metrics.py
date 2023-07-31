import construction as cs
import networkx as nx
import numpy as np
from scipy import stats


def compute_average_ks2(metric,graphs_in,graphs_gen,weight=None):
    if not weight == None:
        metrics_in = metric(graphs_in,weight)
        metrics_gen = metric(graphs_gen,weight)
    else:
        metrics_in = metric(graphs_in)
        metrics_gen = metric(graphs_gen)

    res = [[],[]]
    if metrics_in == [] or metrics_gen == []:
        s = 0
        p = 0
    else:
        s,p = stats.ks_2samp(metrics_in,metrics_gen)
    res[0].append(s)
    res[1].append(p)

    s,p = np.mean(res,-1)
    return s,p


def compute_all_metrics_comphy(graphs_in,graphs_gen):
    #v_density,p_density = compute_average_ks2(closeness,graphs_in,graphs_gen)
    #print("closeness \t %.3f %f"% (v_density,p_density))
    #v_glo_c,p_glo_c = compute_average_ks2(betweenness,graphs_in,graphs_gen)
    #print("betweenness\t %.3f %f"% (v_glo_c,p_glo_c)) 
    #v_aspl,p_aspl = compute_average_ks2(betweenness,graphs_in,graphs_gen,True)
    #print("W betweenness p\t %.3f %f"% (v_aspl,p_aspl))
    v_p,p_p = compute_average_ks2(hour_closeness,graphs_in,graphs_gen)
    print("hour_closeness \t %.3f %f"% (v_p,p_p))
    v_n,p_n = compute_average_ks2(hour_betweenness,graphs_in,graphs_gen,False)
    print("hour_betweenness \t %.3f %f"% (v_n,p_n))
    v_ass,p_ass = compute_average_ks2(hour_betweenness,graphs_in,graphs_gen,True)
    print("W hour_betweenness \t %.3f %f"% (v_ass,p_ass))
    
    v_sm,p_sm = compute_average_ks2(nb_connected_components,graphs_in,graphs_gen)
    print("nb_connected_components\t %.3f %f"% (v_sm,p_sm))   
    v_f,p_f = compute_average_ks2(hour_modularity,graphs_in,graphs_gen)
    print("hour_modularity f \t %.3f %f"% (v_f,p_f))
    v_str,p_str = compute_average_ks2(hour_clustering,graphs_in,graphs_gen)
    print("hour_clustering \t %.3f %f"% (v_str,p_str))
    v_dur,p_dur = compute_average_ks2(hour_assort,graphs_in,graphs_gen)
    print("hour_assort \t %.3f %f"% (v_dur,p_dur))

         
    
    return [(v_p,p_p),(v_n,p_n),(v_ass,p_ass),
            (v_sm,p_sm),(v_f,p_f),(v_str,p_str),(v_dur,p_dur)]
#compute_all_metrics
def compute_all_metrics(graphs_in,graphs_gen):
    v_density,p_density = compute_average_ks(density,graphs_in,graphs_gen)
    print("density \t %.3f %f"% (v_density,p_density))
    v_glo_c,p_glo_c = compute_average_ks(global_clustering,graphs_in,graphs_gen)
    print("global clu\t %.3f %f"% (v_glo_c,p_glo_c)) 
    v_aspl,p_aspl = compute_average_ks(average_shortest_path,graphs_in,graphs_gen)
    print("avg short p\t %.3f %f"% (v_aspl,p_aspl))
    v_p,p_p = compute_average_ks(dist_number_of_individuals,graphs_in,graphs_gen)
    print("metric p \t %.3f %f"% (v_p,p_p))
    v_n,p_n = compute_average_ks(dist_number_of_new_conversations,graphs_in,graphs_gen)
    print("metric n \t %.3f %f"% (v_n,p_n))
    v_ass,p_ass = compute_average_ks(get_ass,graphs_in,graphs_gen)
    print("assortat \t %.3f %f"% (v_ass,p_ass))
    v_sm,p_sm = compute_average_ks(s_metric,graphs_in,graphs_gen)
    print("s metric\t %.3f %f"% (v_sm,p_sm))   
    v_f,p_f = compute_average_ks(dist_frequency_of_interactions,graphs_in,graphs_gen)
    print("metric f \t %.3f %f"% (v_f,p_f))
    v_str,p_str = compute_average_ks(dist_strength_of_nodes,graphs_in,graphs_gen)
    print("strenght \t %.3f %f"% (v_str,p_str))
    v_dur,p_dur = compute_average_ks(dist_duration,graphs_in,graphs_gen)
    print("duration \t %.3f %f"% (v_dur,p_dur))

     
    
    return [(v_density,p_density),(v_glo_c,p_glo_c),(v_aspl,p_aspl),(v_p,p_p),(v_n,p_n),(v_ass,p_ass),
            (v_sm,p_sm),(v_f,p_f),(v_str,p_str),(v_dur,p_dur)]



def compute_average_ks(metric,graphs_in,graphs_gen,weight=False):
    if weight:
        metrics_in = metric(graphs_in,True)
        metrics_gen = metric(graphs_gen,True)
    else:
        metrics_in = metric(graphs_in)
        metrics_gen = metric(graphs_gen)

    res = [[],[]]
    if metrics_in == [] or metrics_gen == []:
        s = 0
        p = 0
    else:
        s,p = stats.ks_2samp(metrics_in,metrics_gen)
    res[0].append(s)
    res[1].append(p)

    s,p = np.mean(res,-1)
    return s,p




# <p>
def dist_number_of_individuals(graphs):
    nb_individuals = []
    for g in graphs:
        individuals = number_of_individuals(g)
        nb_individuals.append(len(individuals))        
    return(nb_individuals)

# used for <p>
def number_of_individuals(g):
    individuals = []
    conn_comp = [x for x in list(nx.connected_components(g)) if len(x)>1]
    for comp in conn_comp:
        for node in comp:
            individuals.append(node)
    return individuals



# assortativity
def get_ass(graphs):
    ass = []
    for g in graphs:
        if len(g.edges)>0:
            assort = nx.degree_assortativity_coefficient(g)
            if not np.isnan(assort):
                ass.append(assort)
    return ass

# <f>
def dist_frequency_of_interactions(graphs):
    nb_interactions = []
    for g in graphs:
        nb_interactions.append(len(g.edges))
    return(nb_interactions)






# <n>
def dist_number_of_new_conversations(graphs):
    
    nb_new_convs = []
    for i in range(len(graphs)-1):
        a1 = nx.adjacency_matrix(graphs[i])
        a2 = nx.adjacency_matrix(graphs[i+1])

        unique, counts = np.unique(np.asarray((a2-a1).A).reshape(-1), return_counts=True)
        d = dict(zip(unique, counts))
        if 1 in d:
            nb_new_convs.append(d[1]/2)
    return(nb_new_convs)
# matricxe di adj, nota che è sugli archi non punto di vista dei nodi




# strength_of_nodes
def dist_strength_of_nodes(graphs):
    G = get_weighted_graph(graphs)
    strength = [e[2]["weight"] for e in G.edges(data=True)]
    return(strength)

def get_weighted_graph(graphs):
    G = nx.Graph()
    for g in graphs:
        for u,v in g.edges():
            if (u,v) in G.edges():
                if u > v:
                    G.edges()[(u,v)]["weight"] = G.edges()[(u,v)]["weight"] + 1
                else:
                    G.edges()[(u,v)]["weight"] = G.edges()[(u,v)]["weight"] + 1
            else:
                G.add_edge(u,v,weight=1)                
    return G






def dist_duration(graphs):
    dict_edges = dict()

    for g in graphs:
        for e in g.edges():
            if not e in dict_edges:
                dict_edges[e] = []


    for g in graphs:
        edges = list(g.edges())
        for k in dict_edges:
            if k in edges:
                dict_edges[k].append(1)
            else:
                dict_edges[k].append(0)


    for k,v in dict_edges.items():
        old = v[0]
        if old == 1:
            array = [1]
        else:
            array = []

        for i in v[1:]:
            if not i == 0:
                if old == 0:
                    array.append(1)
                if old == 1:
                    array[-1] += 1
            old = i

        dict_edges[k] = array


    res = [np.mean(x)for x in list(dict_edges.values())]
    return res








# deg density
def density(graphs):
    res = []
    for g in graphs:
        res.append(nx.density(g))
    return res

# local clustering
# medialo sul tempo
def local_clustering(graphs):
    res = []
    for g in graphs:
        tmp = nx.clustering(g)
        res.append(np.mean(list(tmp.values())))
    return res

# da rifare! 
# global_clustering
def global_clustering(graphs):
    res = []
    for g in graphs:
        res.append(nx.transitivity(g))
    return res
# da rifare! 

# average shortest path
def average_shortest_path(graphs):
    res = []
    for g in graphs:
        largest_cc = max(nx.connected_components(g), key=len) # get the biggest connected components
        sub_G = g.subgraph(largest_cc).copy()  # get subgrph
        res.append(nx.average_shortest_path_length(sub_G)) # get average shortest_path

    res = [i for i in res if i > 0]
    return res


# s metric come in 16 dymond
def s_metric(graphs):
    res = []
    for g in graphs:
        res_in = 0
        for i,j in g.edges():
            d_i = g.degree(i)
            d_j = g.degree(j)
            res_in = res_in + d_i*d_j
        res.append(res_in)

    return res




def closeness(graphs):
    G_aggr = get_weighted_graph(graphs)
    c = nx.closeness_centrality(G_aggr)
    return c

def betweenness(graphs,weighted):
    G_aggr = get_weighted_graph(graphs)
    if weighted:
        b = nx.betweenness_centrality(G_aggr,weight='weight')
    else:
        b = nx.betweenness_centrality(G_aggr,weight=None)
    return b


def hour_closeness(graphs):
    c = []
    tot_hours_dataset = len(graphs)*300/3600
    for hour_n in range(int(tot_hours_dataset)):
        G_aggr = get_weighted_graph(graphs[int(3600/300)*hour_n:int(3600/300)*(hour_n+1)])
        c_dict = nx.closeness_centrality(G_aggr)
        if c_dict == {}:
            c.append(0)
        else:
            c.append(np.mean(list(c_dict.values())))

    return c


def hour_betweenness(graphs,weighted):
    c = []
    tot_hours_dataset = len(graphs)*300/3600
    for hour_n in range(int(tot_hours_dataset)):
        G_aggr = get_weighted_graph(graphs[int(3600/300)*hour_n:int(3600/300)*(hour_n+1)])
        if weighted:
            c_dict = nx.betweenness_centrality(G_aggr,weight='weight')
        else:
            c_dict = nx.betweenness_centrality(G_aggr,weight=None)
        if c_dict == {}:
            c.append(0)
        else:
            c.append(np.mean(list(c_dict.values())))

    return c

def nb_connected_components(graphs):
    nc = []
    for g in graphs:
        nc.append(nx.number_connected_components(g))
    return nc


def hour_modularity(graphs):
    m_list = []
    tot_hours_dataset = len(graphs)*300/3600

    for hour_n in range(int(tot_hours_dataset)):
        G_aggr = get_weighted_graph(graphs[int(3600/300)*hour_n:int(3600/300)*(hour_n+1)])
        if list(G_aggr.nodes()) != []:
            part = nx.community.greedy_modularity_communities(G_aggr,weight='weight')
            m = nx.community.modularity(G_aggr, part,weight='weight')
            m_list.append(m)

    return m_list


def hour_clustering(graphs):
    tot_hours_dataset = len(graphs)*300/3600
    clust = []
    for hour_n in range(int(tot_hours_dataset)):
        G_aggr = get_weighted_graph(graphs[int(3600/300)*hour_n:int(3600/300)*(hour_n+1)])
        clust.append(nx.transitivity(G_aggr))
    return clust

def hour_assort(graphs):
    tot_hours_dataset = len(graphs)*300/3600
    ass = []
    for hour_n in range(int(tot_hours_dataset)):
        G_aggr = get_weighted_graph(graphs[int(3600/300)*hour_n:int(3600/300)*(hour_n+1)])
        assort = nx.degree_assortativity_coefficient(G_aggr)
        if not np.isnan(assort):
            ass.append(assort)
    return ass

