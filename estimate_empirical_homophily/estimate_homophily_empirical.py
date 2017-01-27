"""
The code estimate empirical asymmetric homophily.

The input is number of links amoung minorities and majorities and the fraction of minorities. 

written by: Fariba Karimi
Date: 22-08-2016
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import sys
import pylab
from scipy.optimize import fsolve
import math
import networkx as nx

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
    b_a = float(fa*h_aa)/ (h_aa*ca + h_ab * (2 - ca)) + float(fb*h_ba)/(h_ba*ca + h_bb*(2-ca))

    b_b = float(fb*h_bb)/ (h_ba*ca + h_bb * (2 - ca)) + float(fa*h_ab)/(h_aa*ca + h_ab*(2-ca))


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


#### example ####
min_min = 2294
maj_min = 34620
maj_maj = 335806 

minority_fraction = 0.2



fb = 1 - minority_fraction
fa = minority_fraction

for root_num in [0,1,2]:
	h_aa_anal,h_bb_anal,ca,beta_a,beta_b = fsolve(equations,(1,1,0.5,0.5,0.5))
    degree_exponent_a = float(1)/beta_a + 1
    degree_exponent_b = float(1)/beta_b + 1

	print h_aa_anal,h_bb_anal,ca, beta_a , beta_b, degree_exponent_a , degree_exponent_b

