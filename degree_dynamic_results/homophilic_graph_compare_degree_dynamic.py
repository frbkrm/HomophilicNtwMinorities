##################################################
# plot the time evolution of degree for homophilc BA graph
# 
#
##################################################

import networkx as nx
from collections import defaultdict,Counter
import matplotlib
import matplotlib.pyplot as plt
import random
import bisect
import math
import numpy as np
from generate_homophilic_graph_degree_dynamic import homophilic_barabasi_albert_graph


def analytical_prediction(similitude_param):
    # group a is the minority
    pa = 0.2 #equal ro min fraction
    pb = 1-pa

    daa = similitude_param
    dab = 1- daa

    A = -1
    B = 2 + 2*pa*daa - pa*dab - dab*pb + 2*dab 
    C = -4*pa*daa + 2*pa*dab -4*dab - 2*daa*dab*pa + 2 * (dab**2)*pb
    D = +4*daa*dab*pa

    p = [A,B,C,D]
    print np.roots(p)

    ca = np.roots(p)[1] #only one root is acceptable

    qa = (daa*pa*(2-ca) + ca *dab*pb) / (ca*(2-ca))
    qb = 1 - (pb*(1-qa))/(2-2*qa-pa)

    return qa,qb


if __name__ == '__main__':
    

    m = 2 ; N = 5000 ; N_itteration = 20

    min_fraction_list = [0.2]
    similitude_list = [0.005,0.2,0.5,0.8,0.999,1.0]

    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 20}

    matplotlib.rc('font', **font)


    plot_counter = 0
    for similitude in similitude_list: 
        for min_fraction in min_fraction_list:

            time_dict_minority = defaultdict(list)
            time_dict_majority = defaultdict(list)

            sum_degree_dict = defaultdict(list) #to keep track of C, key = time, value = q*k / t


            for i in range(N_itteration):
                
                G , degree_dynamic , sum_degree_time = homophilic_barabasi_albert_graph(N,m, min_fraction, similitude)

                #print degree_dynamic[10]


                for n in degree_dynamic.keys():
                    #time_list = {k: v for k, v in degree_dynamic[n].iteritems() if v != 0} #remove all degree=0 values
                    #sorted_keys = sorted(time_list)
                    sorted_keys = degree_dynamic[n].keys() #time
                    if len(sorted_keys) < 10: continue
                    offset = 10*sorted_keys[0]
                    if  offset in sorted_keys  : #first 100 nodes

                        for t in sorted_keys:
                            if G.node[n]['color']== 'red':
                                time_dict_minority[t/sorted_keys[0]].append(float(degree_dynamic[n][t]) / degree_dynamic[n][offset])
                            
                            else:
                                time_dict_majority[t/sorted_keys[0]].append(float(degree_dynamic[n][t]) / degree_dynamic[n][offset] )
                            
                for t_,total_deg in sum_degree_time.iteritems():
                    sum_degree_dict[t_].append(float(total_deg)) 

            #print 'min_fraction',min_fraction , 'similitude',similitude 
            #for t_ in time_dict_minority.keys():
            #    print t_ ,  np.mean(time_dict_minority[t_])


            ########### degree dynamic ###########
            #fig, ax = plt.subplots()
            #plot_counter += 1
            #plt.figure(plot_counter)
            #binwidth = 10 ; bins=np.arange(min(x), max(x) + binwidth, binwidth)
            ##### minority####            
            y = []
            x = time_dict_minority.keys()
            for k in time_dict_minority.keys():
                v_list = time_dict_minority[k]
                y.append(np.mean(v_list))
                with open('minority_deg_f_%s_sim_%s.txt'%(min_fraction , similitude),'w') as min_out:
                    min_out.write(str(k)+','+str(np.mean(v_list))+'\n')


            #plt.scatter(x,y ,color = 'r', alpha = 0.4  )

            ##### majority ####            

            y = []
            x = time_dict_majority.keys()
            for k in time_dict_majority.keys():
                v_list = time_dict_majority[k]
                y.append(np.mean(v_list))
                with open('majority_deg_f_%s_sim_%s.txt'%(min_fraction , similitude),'w') as maj_out:
                    maj_out.write(str(k)+','+str(np.mean(v_list))+'\n')

            # plt.scatter(x,y , color = 'b' , alpha = 0.4 )

            # ##### prediction####            
            # q_min , q_maj = analytical_prediction(similitude)


            # x = np.arange(1,N)
            # beta_min = q_min 
            # #offset_min = np.mean(time_dict_minority[10.0]) - (10**beta_min)
            # offset_min = 10**(-beta_min)

            # predicted_value = (x**beta_min) * offset_min

            # plt.plot(x, predicted_value , label=r' $\beta_min$ = %.2f'%beta_min , linestyle ='--' , color = '#f1a340' , linewidth=2)



            # beta_maj = q_maj 
            # offset_maj = 10**(-beta_maj)

            # predicted_value = (x**beta_maj) * offset_maj

            # plt.plot(x, predicted_value , label=r'analytical $\beta_maj$ = %.2f'%beta_maj , linestyle ='--' , color = '#998ec3' , linewidth=2)

            # plt.xlim(10,100)
            # plt.ylim(1,50)

            # plt.title("minority fraction = %s ; similitude = %s"%(min_fraction , similitude))
            # plt.ylabel(r"$k(t)/k(10 t0)$")
            # plt.xlabel(r"$t/t0$")
            # plt.legend(loc='upper left')

            # plt.yscale('log')
            # plt.xscale('log')
            # plt.savefig('degree_dynamic_f_%s_sim_%s.pdf'%(min_fraction , similitude))
            # plt.savefig('degree_dynamic_f_%s_sim_%s.svg'%(min_fraction , similitude))



