#!/usr/bin/python
"""
The code estimates the analytical exponent
and also measures the degree exponent using power_law fit from simulations.

It plots the exponents versus homophily and group size.

written by: Fariba Karimi
Date: 22-01-2016
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import powerlaw
from collections import Counter
import matplotlib as mpl
params = {
   'axes.labelsize': 18,
   'text.fontsize': 18,
   'legend.fontsize': 16,
   'xtick.labelsize': 16,
   'ytick.labelsize': 16,
   'text.usetex': False,
   'font.family': 'serif',
   #'font.serif': ['DejaVu Serif'],
   'figure.figsize': [15, 6]
   }
mpl.rcParams.update(params)

def analytic_exponent(min_fraction,similitude):

    # a majoirity ; b minority
    f = 1-min_fraction
    q = similitude

    alpha = -(2*q-1)**2 
    beta  = (2*q-1)*(2*f*q + 5*q - 3)
    gamma = -2*(4*f*q**2 - 2*f*q + 3*q**2 -4*q + 1) 
    delta = -4*f*q*(1-q) 

    Ca,Cb = [],[]
    Ga,Gb = [],[]

    R = []
    roots = np.roots([alpha,beta,gamma,delta])
    for r in roots:
        if (r > 0) and (r < 2):
            c = r

    ca = f*q/(q*c + (1-q)*(2-c)) + (1-f)*(1-q)/((1-q)*c + q*(2-c))
    cb = (1-f)*q/((1-q)*c + q*(2-c)) + f*(1-q)/(q*c + (1-q)*(2-c))

    #Ca.append(f*q/(q*c + (1-q)*(2-c)) + (1-f)*(1-q)/((1-q)*c + q*(2-c)))
    #Cb.append((1-f)*q/((1-q)*c + q*(2-c)) + f*(1-q)/(q*c + (1-q)*(2-c)))

    Ga =  (float(1)/ca) + 1 
    Gb =  (float(1)/cb) + 1 

    exponent_maj = Ga
    exponent_min = Gb

    return exponent_maj , exponent_min



def main():
# # Plot it
    itter_num = 10
    fig, (ax1,ax2) = plt.subplots(nrows=1, ncols = 2) #sharex = True sharey=True

    colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired   
    colors = [colormap(i) for i in np.linspace(0, 1,6)]

    color_index = 0
    for minority_fraction in [0.1,0.2,0.3,0.4,0.5]:
        y_list = [] #analytical

        y_list_maj = [] #analytical


        x_list = np.linspace(0, 1,20)

        for similitude in  np.linspace(0, 1,20):

            analytic_maj , analytic_min = analytic_exponent(minority_fraction,similitude)
            y_list.append(-analytic_min)
            y_list_maj.append(-analytic_maj)

        ax1.plot(x_list,y_list , label = '%s'%minority_fraction , linestyle='--' , color = colors[color_index] ,linewidth=2 )
        ax2.plot(x_list,y_list_maj , linestyle='--' , color = colors[color_index] ,linewidth=2 )

        color_index += 1

    color_index = 0
    count_fig = 0
    for minority_fraction in [0.1,0.2,0.3,0.4,0.5]:
        x_list = []
        z_list = [] #numerical

        z_list_maj = [] #numerical
        
        
        for i in [0,1,2,3,4,5,6,7,8,9,10]:
            count_fig += 1
            #plt.figure(count_fig)

            homophily = float(i)/10

            min_deg_list = []
            maj_deg_list = []
            fileIn_min = open('../deg_dist_results/deg_dist_f_%s_s_%s_min.txt'%( minority_fraction,homophily) , 'r')

            for line in fileIn_min:
                line = line.strip()
                deg,cnt_all = line.split(' ')
                #cnt = int(cnt_all)/itter_num
                #if cnt > 0:
                for i in range(int(cnt_all)):
                    min_deg_list.append(int(deg))

            fileIn_maj = open('../deg_dist_results/deg_dist_f_%s_s_%s_maj.txt'%( minority_fraction,homophily) , 'r')

            for line in fileIn_maj:
                line = line.strip()
                deg,cnt_all = line.split(' ')
                cnt = int(cnt_all)/itter_num
                for i in range(int(cnt_all)):
                    maj_deg_list.append(int(deg))


            #to find xmax
            max_min = Counter(min_deg_list).keys()[-1]
            min_min = Counter(min_deg_list).keys()[14] # 14



            if homophily > 0.5:
                min_min = Counter(min_deg_list).keys()[10] #10

            if minority_fraction == 0.3:
                if homophily == 0.8:
                    min_min = Counter(min_deg_list).keys()[4]
            #    min_min = Counter(min_deg_list).keys()[6]


            if minority_fraction == 0.1:
                min_min = Counter(min_deg_list).keys()[4]

            fit_min = powerlaw.Fit(min_deg_list , xmax = max_min , xmin = min_min  ) #discrete=True
            alpha_min = fit_min.power_law.alpha

            ####################### for majority fitting is oposit ################
            if homophily < 0.5:
                min_maj = Counter(maj_deg_list).keys()[6] #10
            else:
                min_maj = Counter(maj_deg_list).keys()[14] #10

            if minority_fraction == 0.1:
                min_maj = Counter(maj_deg_list).keys()[2] #10

                
            max_maj = Counter(maj_deg_list).keys()[-1]
            fit_maj = powerlaw.Fit(maj_deg_list , xmin = min_maj  )
            alpha_maj = fit_maj.power_law.alpha

            total_deg = sum(Counter(min_deg_list).values())
            normalized_y = [float(x)/ total_deg for x in Counter(min_deg_list).values()]
            #plt.scatter(Counter(min_deg_list).keys(), normalized_y)

            #fit_min.power_law.plot_pdf(min_deg_list,color='r', linestyle='--'  )
            #plt.savefig('deg_dist_f%s_h%s.pdf'%(minority_fraction,homophily))

            x_list.append(homophily)
            z_list.append(-round(alpha_min,4))
            z_list_maj.append(-round(alpha_maj,4))

            #print(minority_fraction,homophily,alpha_min)
        ax1.scatter(x_list,z_list , color = colors[color_index] , s = 20)
        ax2.scatter(x_list,z_list_maj , color = colors[color_index] , s = 20)

        color_index += 1

    ax1.set_xlabel("homophily " r'$(h)$' , fontsize= 18)
    ax1.set_ylabel("degree exponent " r'$\gamma(h)$' , fontsize= 18)
    ax1.grid(True)
    ax1.legend()
    ax1.text(.08,.05,'A) Minority',verticalalignment='bottom' , fontsize= 20 , transform=ax1.transAxes)
    ax1.set_xlim([0,1])

    ax2.set_xlabel("homophily " r'$(h)$'  , fontsize= 18)
    ax2.set_ylabel("degree exponent " r'$\gamma(h)$'  , fontsize= 18)
    ax2.grid(True)
    #ax2.legend()
    ax2.text(.7,.05,'B) Majority',verticalalignment='bottom' , fontsize= 20 , transform=ax2.transAxes)
    ax2.set_xlim([0,1])


    plt.tight_layout()

    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
              fancybox=True, ncol=5)
    fig.savefig('numerical_analytic_exponent_all.svg')
    fig.savefig('numerical_analytic_exponent_all.pdf')

    #ax2.set_xlabel("homophily " r'(h)')
    #ax2.set_ylabel("exponent of degree distribution " r'$1 + 1/ \beta(\bar{h})$')
    #ax2.grid(True)
    #ax2.legend(loc="lower left")
    #fig2.savefig('minority_degree_exponent.eps')


if __name__ == '__main__':
    main()
