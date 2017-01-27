"""
The code estimate the analytical homophily based on data that is produced for synthetic networks
by varying homophily and minority size. The analytical estimation is then plotted versus the numerical results.

written by: Fariba Karimi
Date: 22-01-2016
"""
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
import sys
import pylab
from scipy.optimize import fsolve
import math

params = {
 'font.family': 'serif',
 'font.serif': ['Serif'],
 'xtick.labelsize':22,
 'ytick.labelsize':22

}
plt.rcParams.update(params)

plt.gcf().subplots_adjust(bottom=0.15)

#plt.rcParams['xtick.labelsize'] = label_size 


def analytical_prediction(min_fraction , d_ab , d_ba):

    f = min_fraction

    fa = f
    fb = 1 - fa

    h_ab = d_ab
    h_ba = d_ba

    h_aa = 1 - h_ab
    h_bb = 1- h_ba

    A = (h_aa - h_ab)*(h_ba - h_bb)
    B = ((2*h_bb - (1-f) * h_ba)*(h_aa - h_ab) + (2*h_ab - f*(2*h_aa - h_ab))*(h_ba - h_bb))
    C = (2*h_bb*(2*h_ab - f*(2*h_aa - h_ab)) - 2*f*h_ab*(h_ba - h_bb) - 2*(1-f)*h_ba * h_ab)
    D = - 4*f*h_ab*h_bb
    p = [A,B,C,D]

    for c in np.roots(p):
        if 0<c and c < 2:
            ca = c 
    #ca = np.roots(p)[1]
    b_a = float(fa*h_aa)/ (h_aa*ca + h_ab * (2 - ca)) + float(fb*h_ba)/(h_ba*ca + h_bb*(2-ca))

    b_b = float(fb*h_bb)/ (h_ba*ca + h_bb * (2 - ca)) + float(fa*h_ab)/(h_aa*ca + h_ab*(2-ca))

    #print 'Mathieu',np.roots(p) , beta_a , (1+(1./beta_a)) , (1+(1./beta_b)) 

    return ca, b_a , b_b

def equations(p):
    h_aa,h_bb,ca,beta_a,beta_b = p

    M = maj_maj + maj_min + min_min
    m_bb = maj_maj
    m_ab = maj_min
    m_aa = min_min

    pbb = float(m_bb)/(m_bb+m_ab + m_aa)
    paa = float(m_aa)/(m_aa + m_ab + m_bb)
    pba = float(m_ab)/(m_aa+m_bb+m_ab)
    pab = pba


    h_ab = 1- h_aa
    h_ba = 1- h_bb

    A = (h_aa - h_ab)*(h_ba - h_bb)
    B = ((2*h_bb - (1-fa) * h_ba)*(h_aa - h_ab) + (2*h_ab - fa*(2*h_aa - h_ab))*(h_ba - h_bb))
    C = (2*h_bb*(2*h_ab - fa*(2*h_aa - h_ab)) - 2*fa*h_ab*(h_ba - h_bb) - 2*(1-fa)*h_ba * h_ab)
    D = - 4*fa*h_ab*h_bb
    P = [A,B,C,D]



    K = beta_b/(beta_a+beta_b)
    Z =1 - K


    return ( (pbb* ((fb * h_bb * K)+ (fa*(1-h_bb)*Z) )) - (fb**2 * h_bb * K ) ,  (paa* ((fa * h_aa * Z)+ (fb*(1-h_aa)*K) ) ) - (fa**2 * h_aa * Z ) 
       , beta_a - float(fa*h_aa)/ (h_aa*ca + h_ab * (2 - ca)) - float(fb*h_ba)/(h_ba*ca + h_bb*(2-ca)) 
        , beta_b - float(fb*h_bb)/ (h_ba*ca + h_bb * (2 - ca)) - float(fa*h_ab)/(h_aa*ca + h_ab*(2-ca)) , ca - np.roots(P)[root_num] )





x_majmaj=[] ; x_minmin = [] ;y_majmaj = [] ; y_minmin = []; y_majmaj_analytic = [] ; y_minmin_analytic = [] 


minority_fraction = float(sys.argv[1])

h_ab_fix = 0.5 #we fix homophily for one group so that we can plot 2 dimentional plots.
h_ba_fix = 0.5

fileName = 'homophily_edges_f_%s_all.txt'%minority_fraction
fileIn = open(fileName,'r')
for line in fileIn:

	line = line.strip()
	h_ab_empiric , h_ba_empiric , r , maj_maj , maj_min , min_min = line.split(' ')
	h_ba_empiric = float(h_ba_empiric)
	h_ab_empiric = float(h_ab_empiric)
	#if h_ab == h_ab_fix:

	maj_maj = float(maj_maj)
	maj_min = float(maj_min)
	min_min = float(min_min)

	total_edge = maj_min + maj_maj + min_min

	#x_minmin.append(float(min_min)/total_edge)

	h_bb_empiric = 1 - h_ba_empiric
	h_aa_empiric = 1- h_ab_empiric

	#y_minmin.append(h_aa_empiric)

	fb = 1 - minority_fraction
	fa = minority_fraction

	h_aa_fix = 1 - h_ab_fix
	# measure analytical h_aa and h_bb
	if h_ab_empiric == h_ab_fix:
		x_majmaj.append(float(maj_maj)/total_edge)
		y_majmaj.append(h_bb_empiric)

		h_bb_list = []
		print '#########################################'
		print 'empirical h_bb', h_bb_empiric

		ca_anal,beta_a_anal,beta_b_anal = analytical_prediction(minority_fraction , h_ab_empiric , h_ba_empiric)
		print 'empirical ca', ca_anal,'beta_a',beta_a_anal,'beta_b',beta_b_anal

		for root_num in [0,1,2]:
			print 'root_num',root_num
			try:
				h_aa_anal,h_bb_anal,ca,beta_a,beta_b = fsolve(equations,(0.1,0.8,0.5,0.5,0.5))
			except IndexError:
				print 'indexError'
				continue
			print 'h_aa_anal',h_aa_anal,'h_bb_anal',h_bb_anal,'ca',ca,'beta_a',beta_a,'beta_b',beta_b

			if beta_b > 1 or beta_a > 1:
				# not acceptable solution
				continue
			#elif h_bb_anal > 1:
			#	continue					
			else:
				h_bb_list.append(h_bb_anal)

		distance_previous = 10
		for h_bb_item in h_bb_list:
			dist = math.fabs(h_bb_item - h_bb_empiric)
			if dist < distance_previous:
				h_bb_accept = h_bb_item
				distance_previous = dist
				print h_bb_item,dist,distance_previous

		y_majmaj_analytic.append(h_bb_accept)

	if h_ba_empiric == h_ba_fix:
		x_minmin.append(float(min_min)/total_edge)
		y_minmin.append(h_aa_empiric)

		h_aa_list = []


		ca_anal,beta_a_anal,beta_b_anal = analytical_prediction(minority_fraction , h_ab_empiric , h_ba_empiric)
		print 'empirical ca', ca_anal,'beta_a',beta_a_anal,'beta_b',beta_b_anal

		for root_num in [0,1,2]:
			print 'root_num',root_num
			try:
				h_aa_anal,h_bb_anal,ca,beta_a,beta_b = fsolve(equations,(0.1,0.8,0.5,0.5,0.5))
			except IndexError:
				print 'indexError'
				continue

			if beta_b > 1 or beta_a > 1:
				# do not accept this solution
				continue
				
			else:
				h_aa_list.append(h_aa_anal)

		distance_previous = 10
		h_aa_accept = h_aa_anal

		for h_aa_item in h_aa_list:
			dist = math.fabs(h_aa_item - h_aa_empiric)
			if dist < distance_previous:
				h_aa_accept = h_aa_item
				distance_previous = dist
				print h_aa_item,dist,distance_previous

		y_minmin_analytic.append(h_aa_accept)


plt.scatter(y_majmaj ,x_majmaj,color = '#347aa3')
plt.plot(y_majmaj_analytic , x_majmaj,color = '#347aa3', ls = '--' , linewidth = 2 , label = r'$h_{aa} = 0.5$')

plt.scatter(y_minmin , x_minmin, color = '#e25115')
plt.plot(y_minmin_analytic,x_minmin,color ='#e25115', ls = '--' , linewidth = 2 , label = r'$h_{bb} = 0.5$')

#####################
# prediction


plt.text(0.08,0.7,r'$E) f_a = %s$'%minority_fraction,fontsize = 26)

plt.legend(loc = 'upper right' , fontsize= 26 , frameon = False)

plt.xlim(0,1.1)
plt.ylim(0,1.1)

#plt.tight_layout()


plt.xlabel(r'$h_{bb} ( h_{aa} )$' ,size = 28)
#plt.ylabel(r'$m_{bb} ( m_{aa} )$ ', size = 28)

plt.savefig('edges_f_%s_h_%s_compare.pdf'%(int(minority_fraction*10),int(h_ab_fix*10)),bbox_inches='tight')


