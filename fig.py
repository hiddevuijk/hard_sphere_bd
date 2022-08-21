import numpy as np
import matplotlib.pyplot as plt
from sys import exit

time_list = np.loadtxt("time.dat")
rho = np.loadtxt("rhoz0.dat")
plt.scatter(rho[:,0], rho[:,1], label=0)
dx = rho[1,0] - rho[0,0]
norm = dx * sum(rho[:,1])
print(norm)


for i_t in time_list:
  continue
  i = int(i_t[0])
  if i % 2 == 0: continue
  t = i_t[1]
  rho = np.loadtxt("rhoz{:d}.dat".format(i))
  plt.scatter(rho[:,0], rho[:,1], label=t)
  dx = rho[1,0] - rho[0,0]
  norm = dx * sum(rho[:,1])
  print(norm)
  
#plt.xlim([-3,3])  
#plt.ylim([0,0.7])  
plt.legend()
plt.show()
