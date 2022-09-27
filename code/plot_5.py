import numpy as np
import matplotlib.pyplot as plt

sizes = 7
plt.rc('font', size=sizes)

data = np.loadtxt('prob5.txt', skiprows=41)
N = data[0:-3,0]
log = data[0:-3,2]

plt.plot(N, log, 'x-')
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('N')
plt.ylabel(r'$\log_{10}$(iteration)$-\log_{10}(N)}$')
plt.savefig('prob5.pdf')
plt.show()
