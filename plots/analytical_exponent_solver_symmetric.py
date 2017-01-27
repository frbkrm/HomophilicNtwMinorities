#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


pa = 0.2 #fraction of minority
pb = 1-pa

daa = 0.8 #distance of minority to minority
dab = 1- daa

A = -1
B = 2 + 2*pa*daa - pa*dab - dab*pb + 2*dab 
C = -4*pa*daa + 2*pa*dab -4*dab - 2*daa*dab*pa + 2 * (dab**2)*pb
D = +4*daa*dab*pa

p = [A,B,C,D]
print np.roots(p)
roots = np.roots(p)
for root in roots:
	ca = root
	qa = (daa*pa*(2-ca) + ca *dab*pb) / (ca*(2-ca))
	qb = 1 - (pb*(1-qa))/(2-2*qa-pa)
	print qa,qb

# # Plot it

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

for pa in np.arange(0.1,0.6,0.1):
	pb = 1 - pa
	y = []
	z = []
	for daa in  np.linspace(0, 1,20):
		dab = dba = 1- daa
		A = -1
		B = 2 + 2*pa*daa - pa*dab - dab*pb + 2*dab 
		C = -4*pa*daa + 2*pa*dab -4*dab - 2*daa*dab*pa + 2 * (dab**2)*pb
		D = +4*daa*dab*pa

		p = [A,B,C,D]
		ca = np.roots(p)[1]
		qa = (daa*pa*(2-ca) + ca *dab*pb) / (ca*(2-ca))
		y.append(qa)
		z.append(-(1+(1./qa)))
#func = lambda daa : (daa*pa*(2-ca) + ca *dab*pb) / (ca*(2-ca))
	x = np.linspace(0, 1,20)
	ax1.plot(x,y , label = 'minority fraction = %s'%pa )
	ax2.plot(x,z , label = 'minority fraction = %s'%pa )

ax1.set_xlabel("homophily " r'(d)')
ax1.set_ylabel("dynamic exponent  " r'$\beta(\bar{q})$')
ax1.grid(True)
ax1.legend()
fig1.savefig('minority_dynamic_exponent.eps')

ax2.set_xlabel("homophily " r'(d)')
ax2.set_ylabel("exponent of degree distribution " r'$1 + 1/ \beta(\bar{q})$')
ax2.grid(True)
ax2.legend(loc="lower left")
fig2.savefig('minority_degree_exponent.eps')

