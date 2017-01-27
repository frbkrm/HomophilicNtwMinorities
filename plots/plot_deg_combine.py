import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import sys
import matplotlib as mpl

params = {
   'axes.labelsize': 16,
   'font.size': 16,
   'legend.fontsize': 18,
   'xtick.labelsize': 18,
   'ytick.labelsize': 18,
   'text.usetex': False,
   'font.family': 'serif',
   #'font.serif': ['Sans-Serif'],
   'figure.figsize': [15, 6]
   }
mpl.rcParams.update(params)


fig, (ax1,ax2) = plt.subplots(nrows=1, ncols = 2) #sharex = True sharey=True

colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired   


############################################################################################################################
# first plot

colors = [colormap(i) for i in np.linspace(0, 1,6)]

color_index = 0
for minority_fraction in [0.1,0.2,0.3,0.4,0.5]:
	x_list = [] ; y_list = []

	for i in [1,2,3,4,5,6,7,8,9,10]:
		homophily = float(i)/10

		min_deg_list = []
		maj_deg_list = []
		fileIn_min = open('../deg_dist_results/deg_dist_f_%s_s_%s_min.txt'%( minority_fraction,homophily) , 'r')

		for line in fileIn_min:
			line = line.strip()
			deg,cnt = line.split(' ')
			for i in range(int(cnt)):
				min_deg_list.append(int(deg))

		fileIn_maj = open('../deg_dist_results/deg_dist_f_%s_s_%s_maj.txt'%( minority_fraction,homophily) , 'r')

		for line in fileIn_maj:
			line = line.strip()
			deg,cnt = line.split(' ')
			for i in range(int(cnt)):
				maj_deg_list.append(int(deg))


		x_list.append(homophily)
		total_deg = sum(min_deg_list) + sum(maj_deg_list)
		mean_deg_min = float(sum(min_deg_list))/total_deg
		y_list.append(mean_deg_min)

	ax1.plot(x_list, y_list , label = '%s'%minority_fraction , color = colors[color_index] , lw = 2 )

	ax1.axhline(minority_fraction, linestyle='dashed' , color = 'gray', linewidth=2)

	color_index += 1
#plt.legend( loc = 1 , prop={'size':12}) 

ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.),
          fancybox=True, ncol=5)

#ax1.set_ylim([0,0.55])
ax1.set_xlabel(r'homophily $(h)$')
ax1.set_ylabel('fraction of minorities total degree')
ax1.text(.15,.05,'A',verticalalignment='bottom' , fontsize= 20)
ax1.text(.15,0.64,'minority fraction',verticalalignment='top' , fontsize= 20)

############################################################################################################################
############################################################################################################################
# second plot

# diverging colors from colorbrewer

colors = [
'#543005',
'#8c510a',
'#bf812d',
'#dfc27d',
'#f6e8c3',
'#c7eae5',
'#80cdc1',
'#35978f',
'#01665e',
'#003c30'
]
#for minority_fraction in [0.1,0.2,0.3,0.4,0.5]: 
for minority_fraction in [0.2]:

    #N = 5000 #number of nodes
    color_index = 0
    for i in [1,2,3,4,5,6,7,8,9,10]:
        homophily = float(i)/10
        print homophily

        x_list = [] ; y_list = [] ; y_list_maj = []

        rank_dict = defaultdict(dict)
        fileIn = open('deg_rank_corrected/minority_degrank_index_%s_%s.txt'%(minority_fraction , homophily) , 'r')

        for line in fileIn:
            line = line.strip()
            if not line.startswith('#'):
                line = line.strip()
                rnk_index,nr_min , nr_maj  = line.split(',')
                rank_dict[int(rnk_index)]['min'] = int(nr_min)
                rank_dict[int(rnk_index)]['maj'] = int(nr_maj)

        N = len(rank_dict.keys())

        for rank in [10,15,20,25,30,35,40,45,50,60,70,80,90,100]:
            percentage = rank / 100.
            print 'percentage',percentage
            x_list.append(percentage)
            count_min = 0
            count_maj = 0
            for rank_index in rank_dict.keys():
                count_min += rank_dict[rank_index]['min']
                count_maj += rank_dict[rank_index]['maj']
                if int(rank_index) > (percentage * N):
                    print 'max rank',rank_index
                    break

            y_list.append(float(count_min) / (count_min + count_maj))
            y_list_maj.append(float(count_maj) / (count_min + count_maj))


        ax2.plot(x_list, y_list , label = r'$h$ = %s'%homophily , color = colors[color_index] , lw = 2 )

        #ax2.plot(x_list, y_list_maj , label = r'$h$ = %s'%homophily , color = colors[color_index] , lw = 2 )

        color_index += 1
    ax2.axhline(minority_fraction, linestyle='dashed' , color = 'gray', linewidth=1)

    
    ax2.legend( loc = 1 , fontsize = 'small' ) 

    ax2.set_xlabel(r"top $d\%$ degree rank")
    ax2.set_ylabel(r"fraction of minorities in top $d\%$")

    #ax2.text(.80,.10,'B',verticalalignment='bottom' , fontsize= 20)

    ax2.text(0.3,.80,r'$f_{a}$ = %s'%minority_fraction,verticalalignment='bottom' , fontsize= 18)

    ax2.set_ylim([0,1])



fig.tight_layout()

fig.savefig('compare_degree_rank_all.pdf')
fig.savefig('compare_degree_rank_all.svg')
