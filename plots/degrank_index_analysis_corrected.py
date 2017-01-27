from __future__ import print_function
import sys
import random
from collections import defaultdict
from collections import Counter
from generate_homophilic_graph import homophilic_barabasi_albert_graph
import numpy as np
import datetime
import bisect
import networkx as nx
import operator

# reading global parametres



if __name__ == "__main__":

    ITTERATION = 100
    N_NODES = 5000
    
    MINORITY_FRACTION = float(sys.argv[1])
    SIMILITUDE = float(sys.argv[2])

    CONNECT_PROB = 2

    current_datetime = datetime.datetime.now() # current date and time
    DateTime_stamp = current_datetime.strftime('%d-%b-%Y, %H:%M:%S') # timestamp inserted in the output file

    for similitude in [SIMILITUDE]:
    
    #for similitude in np.arange(0.0,1.1,0.1):

        OUTPUT_FILENAME = 'minority_degrank_index_%s_%s.txt'%(MINORITY_FRACTION,similitude)
        outputStream = open(OUTPUT_FILENAME, 'w')

        outputStream.write("# date-time: %s\n" % DateTime_stamp)
        outputStream.write("# nr. of nodes: %d\n" % N_NODES)
        outputStream.write("# connect_prob: %d\n" % CONNECT_PROB)
        outputStream.write("# simulation itteration: %d\n" % ITTERATION)
        outputStream.write("# minority fraction: %s\n" % MINORITY_FRACTION)
        outputStream.write("# similitude: %s\n" % similitude)

        outputStream.write("# degrank_index , Nr_minority, Nr_majority  \n" )

        minorities = int(MINORITY_FRACTION * N_NODES ) # give indices for minorities

        #minorities_page_rank = []
        #total_page_rank = 0
        #rank_dict = defaultdict(list)


        rank_dict = defaultdict(list)

        for i in range(ITTERATION):


            G = homophilic_barabasi_albert_graph(N_NODES,CONNECT_PROB , MINORITY_FRACTION , similitude )

            #pagerank_list = nx.pagerank(graph , alpha = 0.85) # alpha = 0.85 ; for undirected graph it assumes two directional edges
            degree_rank = G.degree()
            sorted_degrnk = sorted(degree_rank.items(), key=operator.itemgetter(1) , reverse = True)

            rank_index = 0
            rank_deg = 0
            rank_index_dict = {}
            

            count_all_min = 0 ; count_all_maj = 0
            for node_index , _deg in sorted_degrnk:

                if _deg in rank_index_dict.keys():
                    rank_deg = rank_index_dict[_deg]
                else:
                    rank_deg += 1
                    rank_index_dict[_deg] = rank_deg

                node_color = G.node[node_index]['color']
                if node_color =='red':
                    rank_dict[rank_deg].append('min')
                if node_color =='blue':
                    rank_dict[rank_deg].append('maj')

        for rank_index in rank_dict.keys():
            count_min = rank_dict[rank_index].count('min')
            count_maj = rank_dict[rank_index].count('maj')

            outputStream.write(str(rank_index)+','+str( count_min )+','+str(count_maj) + '\n')

        outputStream.close()
