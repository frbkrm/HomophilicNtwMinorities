import networkx as nx
import powerlaw
import numpy as np
import pylab
pylab.rcParams['xtick.major.pad']='8'
pylab.rcParams['ytick.major.pad']='8'
#pylab.rcParams['font.sans-serif']='Arial'

from matplotlib import rc
#rc('font', family='sans-serif')
rc('font', size=10.0)
#rc('text', usetex=False)

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('legend',**{'fontsize':20})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


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

def analytic_exponent(G,maj_fraction):

    # a is majorities and b is minorities
    assortativity = nx.attribute_assortativity_coefficient(G,'gender')
    daa = 0.5*(assortativity+1)
    dab = 1 - daa

    pa = maj_fraction
    pb = 1- pa

    A = -1
    B = 2 + 2*pa*daa - pa*dab - dab*pb + 2*dab 
    C = -4*pa*daa + 2*pa*dab -4*dab - 2*daa*dab*pa + 2 * (dab**2)*pb
    D = +4*daa*dab*pa
    
    p = [A,B,C,D]
    ca = np.roots(p)[1]

    qa = (daa*pa*(2-ca) + ca *dab*pb) / (ca*(2-ca))
    qb = 1 - (pb*(1-qa))/(2-2*qa-pa)

    exponent_a = (1+(1./qa))
    exponent_b = (1+(1./qb))

    return exponent_a , exponent_b


def main():

	G = nx.read_gexf('DBLP_graph.gexf')


	degree_sequence_min =sorted([d for n,d in G.degree().iteritems() if G.node[n]['gender']== 'f'  ], reverse=True) # degree sequence
	degree_sequence_maj =sorted([d for n,d in G.degree().iteritems() if G.node[n]['gender']== 'm'  ], reverse=True) # degree sequence


	maj_fraction = float(len(degree_sequence_maj))/(len(degree_sequence_min) + len(degree_sequence_maj))
	#e_maj,e_min = analytic_exponent(G,maj_fraction)
	e_maj = 2.94 #measured analytically
	e_min = 3.24

	pdf_female = powerlaw.pdf(degree_sequence_min)
	x_fe =  list(pdf_female[0]) ; y_fe = list(pdf_female[1])
	del x_fe[-1]

	pdf_male = powerlaw.pdf(degree_sequence_maj)
	x_ma =  list(pdf_male[0]) ; y_ma = list(pdf_male[1])
	del x_ma[-1]

	data_min = degree_sequence_min
	plt.scatter(x_fe,y_fe,color=min_color )
	fit_min = powerlaw.Fit(data_min ,xmax = 60 , xmin = 6 )
	alpha_min = fit_min.power_law.alpha

	figFIT = fit_min.power_law.plot_pdf(color=min_color, linestyle='--' , label=r"(women) Fit: %.2f, Model: %.2f"%(alpha_min,e_min) )

	data_maj = degree_sequence_maj

	plt.scatter(x_ma,y_ma,color=maj_color )
	fit_maj = powerlaw.Fit(data_maj ,xmax = 80 , xmin = 6)
	alpha_maj = fit_maj.power_law.alpha

	fit_maj.power_law.plot_pdf(color=maj_color, linestyle='--' , label=r"(men) Fit: %.2f , Model: %.2f"%(alpha_maj,e_maj) )


	plt.ylim([10**(-5),1])
	plt.xlim([1,120])

	handles, labels = figFIT.get_legend_handles_labels()
	leg = figFIT.legend(handles, labels, loc=3 )
	leg.draw_frame(False)


	figFIT.text(7,0.2,'B) Scientific collaboration',fontsize = 24)
	figFIT.text(2,0.01,r'$f_{min} = 0.23$',fontsize = 26)


	figFIT.set_ylabel(r"$p(k)$" , size = 26)
	figFIT.set_xlabel(r"degree $(k)$" , size = 26)

	figname = 'deg_dist_DBLP_2010_scatter'
	plt.savefig(figname+'.eps', bbox_inches='tight')
	plt.savefig(figname+'.svg', bbox_inches='tight')

if __name__ == '__main__':
	main()