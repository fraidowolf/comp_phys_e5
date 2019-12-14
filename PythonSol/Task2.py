import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.integrate import trapz


#Initialize
a=0
b=20
n=1000
h=(b-a)/n
r=np.arange(a+h,a+n*h+h,h)

#Solve the kohn-sham equation for hydrogen
A = -2*np.identity(n)+np.diag(np.ones(n-1),1)+np.diag(np.ones(n-1),-1)
A=-A/(2*h*h)-(1/r)*np.identity(n)
D,V = eigh(A)

#Transform to wavefunction
phi=V[:,0]/np.sqrt(trapz(V[:,0]**2,r))
phi = phi/(np.sqrt(4*np.pi)*r)
E=D[0]

print('Energy: %2.2f eV'%(E*27.211396132))

#Plot and compare to analytical expression
plt.plot(r,phi)
plt.plot(r,2/np.sqrt(4*np.pi)*np.exp(-r),'--')
plt.legend(['Numerical solution', 'Analytical expression'])
plt.xlabel('r[$\mathrm{Å}$]')
plt.ylabel(r'$\varphi(r)[Å^{-3/2}]$')
plt.title('Wavefunction for hydrogen')
plt.savefig('wave_hydro.png')
plt.show()
