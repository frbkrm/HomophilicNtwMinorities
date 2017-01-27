import networkx as nx
import powerlaw

import pylab
pylab.rcParams['xtick.major.pad']='8'
pylab.rcParams['ytick.major.pad']='8'
#pylab.rcParams['font.sans-serif']='Arial'

from matplotlib import rc
#rc('font', family='sans-serif')
rc('font', size=20.0)
#rc('text', usetex=False)

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('legend',**{'fontsize':20})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

from matplotlib.font_manager import FontProperties

panel_label_font = FontProperties().copy()
panel_label_font.set_weight("bold")
panel_label_font.set_size(26.0)
panel_label_font.set_family("Sans-Serif")

import matplotlib.pyplot as plt
params = {
 'font.family': 'Sans-Serif',
 'font.serif': ['Times'],
 'xtick.labelsize':24,
 'ytick.labelsize':24

}

plt.rcParams.update(params)


maj_color = '#347aa3'
min_color = '#e25115'



def main():

	G = nx.Graph()
	fileIn = 'sampled_APS_pacs052030.gexf'
	G = nx.read_gexf(fileIn)

	degree_sequence_maj =sorted([d for n,d in G.degree().iteritems() if G.node[n]['pacs']== '05.30.-d'  ], reverse=True) # degree sequence
	degree_sequence_min =sorted([d for n,d in G.degree().iteritems() if G.node[n]['pacs']== '05.20.-y'  ], reverse=True) # degree sequence

	maj_fraction = float(len(degree_sequence_maj))/(len(degree_sequence_maj) + len(degree_sequence_min))


	maj_maj = maj_min = min_min = 0

	for e in G.edges():
	    n1,n2 = e
	    g1 = G.node[n1]['pacs']
	    g2 = G.node[n2]['pacs']

	    if g1 == g2:
	        if g1 == '05.20.-y':
	            min_min += 1
	        else:
	            maj_maj += 1
	    else:
	        maj_min += 1

	h_aa = 0.88 ; h_bb = 1 ; e_min = 3.4 ; e_maj = 2.8

	data_min = degree_sequence_min
	data_pdf = powerlaw.pdf(data_min)
	x = list(data_pdf[0])
	y = data_pdf[1]
	del x[-1]
	#scatter(x,y,color='r' )
	#plot(range(10),range(10))

	data_maj = degree_sequence_maj

	maj_pdf = powerlaw.pdf(data_maj)
	x_maj = list(maj_pdf[0])
	del x_maj[-1]
	y_maj = maj_pdf[1]

	data_min = degree_sequence_min
	data_pdf = powerlaw.pdf(data_min)

	plt.scatter(x,y,color=min_color )
	#figPDF = powerlaw.plot_pdf(data_min,color='r', linewidth=2  )
	fit_min = powerlaw.Fit(data_min ,xmax = 30 , xmin = 4 )
	alpha_min = fit_min.power_law.alpha

	figFIT = fit_min.power_law.plot_pdf(color=min_color, linestyle='--' , label=r"(CSM) Fit: %.2f, Model: %.2f"%(alpha_min,e_min)  )

	data_maj = degree_sequence_maj

	plt.scatter(x_maj,y_maj,color=maj_color )

	#powerlaw.plot_pdf(data_maj,color='b', linewidth=2  )
	fit_maj = powerlaw.Fit(data_maj ,xmax = 80 , xmin = 5)
	alpha_maj = fit_maj.power_law.alpha

	fit_maj.power_law.plot_pdf(color=maj_color, linestyle='--' , label=r"(QSM) Fit: %.2f , Model: %.2f"%(alpha_maj,e_maj)  )


	handles, labels = figFIT.get_legend_handles_labels()
	leg = figFIT.legend(handles, labels, loc=3 )
	leg.draw_frame(False)

	plt.xlim([1,40])

	figFIT.text(7,0.18,'C) Scientific citation',fontsize = 24)
	figFIT.text(2,0.01,r'$f_{min} = 0.38$',fontsize = 26)


	figFIT.set_ylabel(r"$p(k)$" , size = 26)
	figFIT.set_xlabel(r"degree $(k)$" , size = 26)

	figname = 'deg_dist_APS_05'
	plt.savefig(figname+'.eps', bbox_inches='tight')
	plt.savefig(figname+'.svg', bbox_inches='tight')

if __name__ == '__main__':
	main()