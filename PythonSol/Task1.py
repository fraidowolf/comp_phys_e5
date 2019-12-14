import numpy as np
import matplotlib.pyplot as plt

#Initialize
a=0
b=100
n=1000
h=(b-a)/n
r=np.arange(a + h,a+n*h + h,h)

# Input electron density
b=-4*r*np.exp(-2*r)

#Solve the poisson equation
A = -2*np.identity(n)+np.diag(np.ones(n-1),1)+np.diag(np.ones(n-1),-1)
A=A/(h*h)
U=np.linalg.solve(A,b)

#Analytical expression
V=1/r-(1+1/r)*np.exp(-2*r)
plt.plot(r, U)
plt.show()

#Plot potentials
plt.plot(r,27.211396132*(U+r/r[-1])/r)
plt.plot(r,27.211396132*V,'--')
plt.legend(['Numerical solution', 'Analytical expression'])
plt.xlabel('r[Å]')
plt.ylabel(r'$V_H(r)[\mathrm{eV}]$')
plt.title('Hartree potential for hydrogen atom')
plt.savefig('hydrogen_pot.png')
plt.show()
