import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("data.csv", delimiter='\t')

dr = data[1,0]-data[0,0]
f = np.zeros_like(data);

for i in range(1,len(f)):
    f[i-1,1]=(1/(dr*i)-(1+1/(dr*i))*np.exp(-2*dr*i))
    f[i-1,0]=data[i,0]

plt.plot(f[:-1,0],f[:-1,1])
plt.plot(data[:,0],data[:,1])
plt.show()
