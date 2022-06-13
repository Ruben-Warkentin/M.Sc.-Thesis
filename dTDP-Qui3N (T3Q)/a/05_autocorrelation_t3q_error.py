#!/usr/bin/env python
# coding: utf-8

# Importing Packages

# In[33]:


import numpy as np
import csv
import matplotlib.pyplot as plt
import pandas as pd


# Global Variables

# In[34]:


# the list of "a" contains the different vector list file names
a = ['a0_55'] #make a separate f(x) that reads the files
file = []
f = []
frames = 10000
max_tau = 600
n_relaxation_time = 300
x_max = max_tau/100
# the list of "a" contains the different vector list file names

file = "v_list_a0_55.txt"
#open the file
#vector list
sv = []
#mean y value list 
y = []
y_mean = []
y_SE = []
#execute the functions
df_y = pd.DataFrame({})


# Functions

# In[35]:


def read_data(): #note that "f" is defined previously as 
	#open the file
	f = open(file)
	reader = csv.reader(f,delimiter=' ',quoting=csv.QUOTE_NONNUMERIC) # change contents to floats, floats needed to use dot product function later on. 
	next(reader, None) #skip headers
	for row in reader: # each row is a list
		sv.append(row) #write each row to the sv array


# In[38]:


def autoc():#read the file, delimiter ' ' indicates spaces separate columns

	##rotational autocorrelation function is C(tau) = <u(t).u(t+tau)>
	for tau in range (0, max_tau, 1):
		y = []
		for t in range (0, frames-tau, 1):
			C = (np.dot(sv[t], sv[t+tau]))
			y.append(C)
		y_mean.append(np.mean(y)) # the y_mean is the y value used for graphing later on
		n = (frames-tau)/n_relaxation_time # n is an (under-)estimate of the number of independent populations from which C is sampled.
		y_SE.append(np.std(y)/np.sqrt(n))
		### calculate error bar for y values
		### function parameters and function return values


# In[39]:


def autoc_graph():
	#x axis
	# x = list(range(0, max_tau))
	x = np.linspace(0,x_max,max_tau)
	yerr = y_SE
	# plotting the points
	plt.plot(x, y_mean)
	plt.errorbar(x, y_mean, yerr=yerr)

	# # naming the x axis
	# plt.xlabel('tau')
	plt.xlabel('time (ns)')

	# # naming the y axis
	plt.ylabel('autocorrelation')

	# # giving a title to my graph
	plt.title('t3q Rotational Autocorrelation: a0-55')

	# # function to show the plot
	# plt.savefig('200_std_{}.png'.format(i), bbox_inches='tight')
	# plt.savefig('200_std_{}.pdf'.format(i), bbox_inches='tight')
	plt.show()


# Executing Functions

# In[40]:


read_data()
autoc()
autoc_graph()

