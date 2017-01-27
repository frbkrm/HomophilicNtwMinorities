#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import powerlaw
from collections import Counter,defaultdict
import matplotlib as mpl
import networkx as nx
import operator

from matplotlib import rc
#rc('font', family='sans-serif')
rc('font', size=20.0)
#rc('text', usetex=False)

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('legend',**{'fontsize':20})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

#from matplotlib.font_manager import FontProperties

#panel_label_font = FontProperties().copy()
#panel_label_font.set_weight("bold")
#panel_label_font.set_size(26.0)
#panel_label_font.set_family("Sans-Serif")

params = {
   'axes.labelsize': 22,
   'text.fontsize': 22,
   'legend.fontsize': 18,
   'xtick.labelsize': 18,
   'ytick.labelsize': 18,
   'font.family': 'Serif',
   'font.serif': ['Times'],
   'figure.figsize': [20, 5.5]
   }
mpl.rcParams.update(params)

maj_color = '#347aa3'
min_color = '#e25115'

def read_graphs():
    G_dblp = nx.Graph()
    G_dblp = nx.read_gexf('../DBLP_ntw/DBLP_graph.gexf')

    for n in G_dblp.nodes():
        if G_dblp.node[n]['gender'] == 'f':
            G_dblp.node[n]['color'] = 'min'

        if G_dblp.node[n]['gender'] == 'm':
            G_dblp.node[n]['color'] = 'maj'

    G_brazil = nx.Graph()

    fileIn = open('../Brazil/brazil.txt','r')
    for line in fileIn:
        m,f,t = line.split(' ')
        G_brazil.add_node(m,color='min')
        G_brazil.add_node(f,color='maj')
        G_brazil.add_edge(m,f)

    G_aps = nx.Graph()
    fileIn = '../APS_data/sampled_APS_pacs052030.gexf'
    G_aps = nx.read_gexf(fileIn)

    for n in G_aps.nodes():
        if G_aps.node[n]['pacs'] == '05.20.-y':
            G_aps.node[n]['color'] = 'min'

        if G_aps.node[n]['pacs'] == '05.30.-d':
            G_aps.node[n]['color'] = 'maj'

    print 'read graphs'
    return G_brazil,G_dblp,G_aps

def rank_degree_function(G):
    degree_rank = G.degree()
    sorted_degrnk = sorted(degree_rank.items(), key=operator.itemgetter(1) , reverse = True)

    rank_index = 0
    rank_deg = 0
    rank_index_dict = {}
    rank_dict = defaultdict(list)
    count_all_min = 0 ; count_all_maj = 0
    for node_index , _deg in sorted_degrnk:

        if _deg in rank_index_dict.keys():
            rank_deg = rank_index_dict[_deg]
        else:
            rank_deg += 1
            rank_index_dict[_deg] = rank_deg

        node_color = G.node[node_index]['color']
        if node_color =='min':
            count_all_min += 1
            rank_dict[rank_deg].append('min')
        if node_color =='maj':
            count_all_maj += 1
            rank_dict[rank_deg].append('maj')

    x_list = [] ; y_list = []
    N = len(rank_dict.keys())

    minority_fraction = float(count_all_min) /(count_all_min+count_all_maj)
    print 'minority_fraction',minority_fraction

    for rank in [5,10,20,30,40,50,60,70,80,90,100]:
    
        percentage = rank / 100.
        x_list.append(percentage)
        count_min = 0
        count_maj = 0
        for rank_index in rank_dict.keys():
            count_min += rank_dict[rank_index].count('min')
            count_maj += rank_dict[rank_index].count('maj')
            if int(rank_index) > (percentage * N):
                break

        y_list.append(float(count_min) / (count_min + count_maj))

    return x_list , y_list , minority_fraction

def main():
# # Plot it
    G_brazil,G_dblp,G_aps = read_graphs()

    itter_num = 10
    fig, (ax1,ax2,ax3) = plt.subplots(nrows=1, ncols = 3) #sharex = True sharey=True

    colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired   
    colors = [colormap(i) for i in np.linspace(0, 1,6)]

    #######################  brazil #############################

    x_list , y_list , minority_fraction = rank_degree_function(G_brazil)
    ax1.plot(x_list,y_list , label = '%s'%minority_fraction , linewidth=2  , color = min_color )
    ax1.axhline(minority_fraction, linestyle='dashed' , color = 'gray', linewidth=3 )
    analytic_file = open('brazil_topk_std.txt','r')
    y_list_analytic = []
    std_list = []
    x_list_plot = []
    for line in analytic_file:
        line = line.strip()
        x,y,y_std = line.split(',')
        if float(x) < 0.05: continue
        y_list_analytic.append(float(y))
        std_list.append(float(y_std))
        x_list_plot.append(float(x))

    print len(y_list_analytic) , len(std_list)
    ax1.errorbar(x_list_plot,y_list_analytic, std_list,linestyle = 'dashed' , color = min_color , linewidth = 2)

    #######################  dblp #############################

    x_list , y_list , minority_fraction = rank_degree_function(G_dblp)
    ax2.plot(x_list,y_list , label = '%s'%minority_fraction ,linewidth=2  , color = min_color )
    ax2.axhline(minority_fraction, linestyle='dashed' , color = 'gray', linewidth=3 )
    analytic_file = open('dblp_topk_std.txt','r')
    y_list_analytic = []
    std_list = []

    for line in analytic_file:
        line = line.strip()
        x,y,y_std = line.split(',')
        if float(x) < 0.05: continue

        y_list_analytic.append(float(y))
        std_list.append(float(y_std))

    ax2.errorbar(x_list_plot,y_list_analytic, std_list,linestyle = 'dashed' , color = min_color , linewidth = 2)

    #######################  aps #############################

    x_list , y_list, minority_fraction = rank_degree_function(G_aps)
    ax3.plot(x_list,y_list ,  linewidth=2 , color = min_color , label = 'data' )
    ax3.axhline(minority_fraction, linestyle='dashed' , color = 'gray', linewidth=3 )
    analytic_file = open('aps_topk_std.txt','r')
    y_list_analytic = []
    std_list = []
    for line in analytic_file:
        line = line.strip()
        x,y,y_std = line.split(',')
        if float(x) < 0.05: continue

        y_list_analytic.append(float(y))
        std_list.append(float(y_std))

    ax3.errorbar(x_list_plot,y_list_analytic, std_list,linestyle = 'dashed' , color = min_color , linewidth = 2 , label = 'model')



    ax1.set_xlabel(r"top $d\%$ degree rank" , fontsize= 22)
    ax2.set_xlabel(r"top $d\%$ degree rank"  , fontsize= 22)
    ax3.set_xlabel(r"top $d\%$ degree rank"  , fontsize= 22)

    ax1.set_ylabel(r"fraction of minorities in top $d\%$" , fontsize= 22)
    #ax2.set_ylabel("p (top k%)" , fontsize= 18)
    #ax3.set_ylabel("p (top k%)" , fontsize= 18)

    ax1.grid(True)

    #ax1.legend()
    ax1.text(.1,.05,'A) Sexual contacts',verticalalignment='bottom' , fontsize= 22 , transform=ax1.transAxes)
    ax2.text(1.2,.9,'B) Scientific collaboration',verticalalignment='bottom' , fontsize= 22 , transform=ax1.transAxes)
    ax3.text(2.4,.9,'C) Scientific citation',verticalalignment='bottom' , fontsize= 22 , transform=ax1.transAxes)

    ax1.set_ylim([0,1])
    ax2.set_ylim([0,1])
    ax3.set_ylim([0,1])

    for ax in [ax1,ax2,ax3]:
        ax.grid(True)

    #ax2.set_xlabel("homophily " r'$(h)$'  , fontsize= 18)
    #ax2.set_ylabel("degree exponent " r'$\beta(h)$'  , fontsize= 18)
    #ax2.grid(True)
    #ax2.legend()
    #ax2.text(.8,.05,'B',verticalalignment='bottom' , fontsize= 20 , transform=ax2.transAxes)
    #ax2.set_xlim([0,1])


    plt.tight_layout()

    ax3.legend(loc='upper right',fancybox=True , fontsize = 20)

    fig.savefig('empirical_analytical_topk_std.svg')
    fig.savefig('empirical_analytical_topk_std.pdf')




if __name__ == '__main__':
    main()
