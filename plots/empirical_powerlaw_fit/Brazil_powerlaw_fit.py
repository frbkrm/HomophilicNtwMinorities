import networkx as nx
import powerlaw
import numpy as np
import pylab


from matplotlib import rc
#rc('font', family='sans-serif')
rc('font', size=20.0)
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

	G = nx.Graph()
	fileIn = open('brazil.txt','r')
	for line in fileIn:
	    m,f,t = line.split(' ')
	    G.add_node(m,gender='m')
	    G.add_node(f,gender='f')
	    G.add_edge(m,f)

	degree_sequence_female =sorted([d for n,d in G.degree().iteritems() if G.node[n]['gender']== 'f'  ], reverse=True) # degree sequence
	degree_sequence_male =sorted([d for n,d in G.degree().iteritems() if G.node[n]['gender']== 'm'  ], reverse=True) # degree sequence



	male_fraction = float(len(degree_sequence_male))/(len(degree_sequence_male) + len(degree_sequence_female))
	print male_fraction
	e_male,e_female = analytic_exponent(G,male_fraction)

	pdf_female = powerlaw.pdf(degree_sequence_female)
	x_fe =  list(pdf_female[0]) ; y_fe = list(pdf_female[1])
	del x_fe[-1]

	pdf_male = powerlaw.pdf(degree_sequence_male)
	x_ma =  list(pdf_male[0]) ; y_ma = list(pdf_male[1])
	del x_ma[-1]

	data_fe = degree_sequence_female

	plt.scatter(x_fe,y_fe,color=maj_color  )
	fit_fe = powerlaw.Fit(data_fe, discrete=True ,xmax = 200 , xmin = 20 )
	alpha_fe = fit_fe.power_law.alpha

	figFIT = fit_fe.power_law.plot_pdf(color=maj_color, linestyle='--' , label=r"(workers) Fit: %.2f, Model: %.2f"%(alpha_fe,e_female) )

	data_ma = degree_sequence_male

	plt.scatter(x_ma,y_ma,color=min_color )
	fit_ma = powerlaw.Fit(data_ma, discrete=False ,xmax = 200 , xmin = 20)
	alpha_ma = fit_ma.power_law.alpha

	fit_ma.power_law.plot_pdf(color=min_color, linestyle='--' , label=r"(buyers) Fit: %.2f, Model: %.2f"%(alpha_ma,e_male) )

	plt.xlim([1,200])


	#plt.tight_layout()

	handles, labels = figFIT.get_legend_handles_labels()
	leg = figFIT.legend(handles, labels, loc=3 )
	leg.draw_frame(False)


	figFIT.text(9,0.18,'A) Sexual contacts',fontsize = 24)
	figFIT.text(2,0.001,r'$f_{min} = 0.40$',fontsize = 26)


	figFIT.set_ylabel("$p(k)$" , size = 26)
	figFIT.set_xlabel("degree $(k)$ " , size = 26)

	figname = 'deg_dist_brazil_scatter'
	plt.savefig(figname+'.eps', bbox_inches='tight')
	plt.savefig(figname+'.svg', bbox_inches='tight')

if __name__ == '__main__':
	main()
