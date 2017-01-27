"""
Generators for BA homophilic graph.

The algorithm is based on paper 'scale-free homophilic network'
by almedia et al (2013).

written by: Fariba Karimi
Date: 27-07-2016
"""

import networkx as nx
from collections import defaultdict
import random
import bisect
import copy


def homophilic_barabasi_albert_graph_asym(N, m , minority_fraction, h_ab , h_ba):
    """Return homophilic random graph using BA preferential attachment model.

    A graph of n nodes is grown by attaching new nodes each with m
    edges that are preferentially attached to existing nodes with high
    degree. The connections are established by linking probability which 
    depends on the connectivity of sites and the similitude (similarities).
    similitude varies ranges from 0 to 1.

    Parameters
    ----------
    N : int
        Number of nodes
    m : int
        Number of edges to attach from a new node to existing nodes
    seed : int, optional
        Seed for random number generator (default=None).

    minority_fraction : float
        fraction of minorities in the network

    similitude: float
        value between 0 to 1. similarity between nodes. if nodes have same attribute
        their similitude (distance) is smaller.

    Returns
    -------
    G : Graph

    Notes
    -----
    The initialization is a graph with with m nodes and no edges.

    References
    ----------
    .. [1] A. L. Barabasi and R. Albert "Emergence of scaling in
       random networks", Science 286, pp 509-512, 1999.
    """
    # a = min 
    #h_ab = 0.8
    #h_ba = 0.2

    h_aa = 1 - h_ab
    h_bb = 1 - h_ba
    # Add m initial nodes (m0 in barabasi-speak)
    G = nx.Graph()

    minority = int(minority_fraction * N)

    minority_nodes = random.sample(range(N),minority)
    node_attribute = {}
    for n in range(N):
        if n in minority_nodes:
            G.add_node(n , color = 'red')
            node_attribute[n] = 'minority'
        else:
            G.add_node(n , color = 'blue')
            node_attribute[n] = 'majority'


    dist = defaultdict(int)

    #create homophilic distance ### faster to do it outside loop ###
    for n1 in range(N):
        n1_attr = node_attribute[n1]
        for n2 in range(N):
            n2_attr = node_attribute[n2]
            if n1_attr == n2_attr:
            	if n1_attr == 'minority':
            		dist[(n1,n2)] = h_aa
            	else:
            		dist[(n1,n2)] = h_bb
            else:
            	if n1_attr == 'minority':
            		dist[(n1,n2)] = h_ab
            	else:
                	dist[(n1,n2)] = h_ba



    target_list=list(range(m))
    source = m

    time_step = 1

    #selective_nodes = random.sample(range(minority),40) + random.sample(range(minority,N),40) # to evaluate degree dynamic
    degree_dynamic = defaultdict(dict)
    sum_degree_time = defaultdict(float) #for calculating sum(qk), keys = time, values = total prob. sum

    existing_nodes = list(range(m))
    while source < N:


        targets , prob_sum = _pick_targets(G,source,target_list,dist,m)
        
        sum_degree_time[time_step] = prob_sum

        if targets != set(): #if the node does  find the neighbor
            G.add_edges_from(zip([source]*m,targets))

        ##### calculate degree dynamic ####
        for n_ in existing_nodes:
            degree_dynamic[n_][time_step] = G.degree(n_)

        target_list.append(source)
        existing_nodes.append(source)

        time_step += 1
        source += 1

    #print('#################################')
    return G , degree_dynamic , sum_degree_time

def _pick_targets(G,source,target_list,dist,m):
    
    prob_sum = 0.0
    target_prob_dict = {}
    for target in target_list:
        target_prob = (1-dist[(source,target)])* (G.degree(target)+0.00001)
        #print target_prob
        #prob_sum += target_prob
        target_prob_dict[target] = target_prob
    
    prob_sum = sum(target_prob_dict.values())

    targets = set()
    target_list_copy = copy.copy(target_list)
    count_looking = 0
    if prob_sum == 0:
        return targets
    while len(targets) < m:
        count_looking += 1
        if count_looking > len(G): # if node fails to find target
            break
        rand_num = random.random()
        cumsum = 0.0
        for k in target_list_copy:
            cumsum += float(target_prob_dict[k]) / prob_sum
            if rand_num < cumsum:
                targets.add(k)
                target_list_copy.remove(k)
                break
    return targets , prob_sum


if __name__ == '__main__':
    graph = homophilic_barabasi_albert_graph_assym(N = 100, m = 2 , minority_fraction = 0.5, h_ab=0.5 , h_ba = 0.5)







