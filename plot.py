import numpy as np
import matplotlib.pyplot as plt
from sys import exit

rho = np.loadtxt("mc_rhoz1.5.dat")
plt.scatter(rho[:,0], rho[:,1], label="Eq. A=1.5")
dx = rho[1,0] - rho[0,0]
norm = dx * sum(rho[:,1])
print(norm)

rho = np.loadtxt("mc_rhoz1.1.dat")
plt.scatter(rho[:,0], rho[:,1], label="Eq. A=1.1")

dx = rho[1,0] - rho[0,0]
norm = dx * sum(rho[:,1])
print(norm)

time_list = np.loadtxt("time.dat")
x = np.loadtxt("rhoz1.dat")[:,0]

rho  = np.loadtxt("rho1.dat")[:,1]
rho += np.loadtxt("rho2.dat")[:,1]
rho += np.loadtxt("rho3.dat")[:,1]
rho += np.loadtxt("rho4.dat")[:,1]
rho += np.loadtxt("rho5.dat")[:,1]

rho /= 5

plt.scatter(x, rho, label=0)
dx = x[1] - x[0]
norm = dx * sum(rho)
print(norm)

A = 1.5
N = 1000
#x = np.linspace( rho[0,0]*1.1, rho[-1,0]*1.1, N)
#y = norm * np.exp( - A * x * x) * np.sqrt( A / np.pi)
#plt.plot(x,y)

A = 1.1
N = 1000
#x = np.linspace( rho[0,0]*1.1, rho[-1,0]*1.1, N)
#y = norm * np.exp( - A * x * x) * np.sqrt( A / np.pi)
#dx = x[1] - x[0]
#norm = dx * sum(y)
#print(norm)
#plt.plot(x,y)


  
#plt.xlim([-3,3])  
#plt.ylim([0,0.7])  
plt.legend()
plt.show()
