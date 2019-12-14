import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.integrate import trapz

a=0
b=5
n=1000
h=(b-a)/n
r=np.arange(a+h,a+n*h+h,h)
b=-4*r*np.exp(-2*r)

n_s=(1/np.pi)*np.exp(-2*r)
A = -2*np.identity(n)+np.diag(np.ones(n-1),1)+np.diag(np.ones(n-1),-1)
A=A/(h*h)
A1=-A/2-(2/r)*np.identity(n)

eps_old = 1
eps_new = 0

while np.abs(eps_old-eps_new)>1e-5:
	b=-4*np.pi*r*n_s
	U=np.linalg.solve(A,b)
	V_H=(U+r/r[-1])/r
	
	A2=A1+np.diag(V_H)
	D,V = eigh(A2)
	phi=V[:,0]/np.sqrt(trapz(V[:,0]**2,r))
	phi = phi/(np.sqrt(4*np.pi)*r)
	n_s=phi**2
	eps_old = eps_new
	eps_new = D[0]*27.211399 #Convert to eV
	
	


    

plt.plot(r,4*np.pi*r**2*n_s, label = "Data")
Z=2
plt.plot(r,4*Z**3*r**2*np.exp(-2*Z*r), label = "Central Field Approx. (Z = 2)")
Z=27/16
plt.plot(r,4*Z**3*r**2*np.exp(-2*Z*r), label = "Central Field Approx. (Z = 27/16)")
plt.legend()
plt.show()
E=(2*D[0]-4*np.pi*trapz(r**2*V_H*n_s,r))#*27.211399
print('Energy = %f eV' % E)