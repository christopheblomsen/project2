#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

sizes = 7
plt.rc('font', size=sizes)


#plot for problem 4, not necessary to include in report unless we want
eig_num = np.loadtxt('eig_vec_num_4.txt')
eig_ana = np.loadtxt('eig_vec_ana_4.txt')
eig_num = eig_num/np.linalg.norm(eig_num)
eig_ana = eig_ana/np.linalg.norm(eig_ana)

#look at plots to find corresponding eigenvectors
fig, [ax_num, ax_ana] = plt.subplots(2, 1)
for i in range(3):
         ax_num.plot(eig_num[:, i], label = f"mode: {i}")
         ax_ana.plot(eig_ana[:, i], label = f"mode: {i} ")
ax_num.set_ylim([-0.5, 0.5])
ax_ana.set_ylim([-0.5, 0.5])
ax_num.set_title('Numeric')
ax_ana.set_title('Analytic')
plt.legend()
#plt.savefig('problem4b.pdf')


#verify that the numeric solutions are eigenvectors,
#this division should give a constant vector
print(eig_num[:, 0]/eig_ana[:, 0])


#plot for problem 6
eig = np.loadtxt('eig_vec_num_6a.txt')
x = np.loadtxt('x_axis_6a.txt')

eig_ana = np.loadtxt("eig_vec_ana_6a.txt")
# We multiply by -1 to get the same eigenvector which is still a solution
# to our system so it is just to get a nicer graph
eig_ana[:, 1] *= -1
eig_ana[:, 2] *= -1


fig = plt.figure()
ax = plt.axes()
for i in range(3):
         ax.plot(x, eig[:, i], label = f"n = {i}")
         ax.plot(x, eig_ana[:, i], 'r--', label = f"analytical n = {i}")
plt.legend()
plt.hlines(0, 0, 1, linestyle='dashed', color='black')
plt.xlabel(r'$\hat{x}$')
plt.ylabel('Height')
plt.savefig('prob_6a.pdf')
plt.show()

eig = np.loadtxt('eig_vec_num_6b.txt')
x = np.loadtxt('x_axis_6b.txt')

eig_ana = np.loadtxt('eig_vec_ana_6b.txt')

fig = plt.figure()
ax = plt.axes()
for i in range(3):
         ax.plot(x, eig[:, i], label = f"n = {i}")
         ax.plot(x, eig_ana[:, i], 'r--', label = f"analytical n = {i}")
plt.hlines(0, 0, 1, linestyle='dashed', color='black')
plt.legend()
plt.xlabel(r'$\hat{x}$')
plt.ylabel('Height')
plt.savefig('prob_6b.pdf')
plt.show()
