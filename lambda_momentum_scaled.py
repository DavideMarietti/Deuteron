#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Reads and prints the Tree structure and what it is containted.

import ROOT


filename_u = '../U1S/lambda_U1S.root'
my_file_u = ROOT.TFile.Open(filename_u)
my_file_u.ls()
my_tree_u = my_file_u.Get('lambda_U1S_tree')
my_tree_u.Print()

filename_c = '../Continuum/lambda_continuum.root'
my_file_c = ROOT.TFile.Open(filename_c)
my_file_c.ls()
my_tree_c = my_file_c.Get('lambda_continuum_tree')
my_tree_c.Print()


# In[3]:


#Translates the Tree into DataFrame, seclects the CMS_p variable and calculates the scaled momentum (CMS_p/E_beam)
#both for continuum and U1S events. Finally gathers these values in 10 bins and calculate the enhancement.

import root_pandas
import numpy as np


df_u = root_pandas.read_root(filename_u, key='lambda_U1S_tree')
p_u = df_u['CMS_p'].to_numpy()
p_u_scaled= p_u / (9.46 * 0.5) #U1S: CMS_E = 9.460 GeV

df_c = root_pandas.read_root(filename_c, key='lambda_continuum_tree')
p_c = df_c['CMS_p'].to_numpy()
p_c_scaled= p_c / (10.58 * 0.5) #Continuum: CMS_E = 10.580 GeV

bins = np.linspace(0, 1, 11) #10 bins
p_c_binned = np.histogram(p_c_scaled, bins)
p_u_binned = np.histogram(p_u_scaled, bins)

p_enhancement = ((p_u_binned[0]) / (p_c_binned[0]) * (98021/100000) - (0.108)) / 0.817
p_enhancement[p_enhancement < 0] = 0


data_enhancement = np.array([2.37096774193548, 3.21505376344086, 2.94623655913978, 2.04838709677419,
                             1.45698924731183, 1.09677419354839, 0.827956989247312, 0, 0, 0])

data_error = np.array([2.9367572894325, 3.53031708306765, 3.0538954263521, 2.12596371224857, 1.54741754101212,
                       1.23240681699612, 0.997985624021796, 0.84418294591618, 0.265694742333778, 0.316083125615906])

error = data_error - data_enhancement

#Visual check
print(data_enhancement)
print(p_enhancement)
print(p_u_binned[0])
print(p_c_binned[0])


# In[4]:


#Plot the momentum-scaled enhancement. Multiple candidates in the same event are considered.

get_ipython().run_line_magic('matplotlib', 'ipympl')
import matplotlib.pyplot as plt


fig, ax = plt.subplots()
new_bins = np.array(bins[:-1]) + 0.5 * (bins[1] - bins[0])
ax.errorbar(new_bins, data_enhancement, yerr=error, fmt='bs', elinewidth=0.7, barsabove=True)
ax.plot(new_bins, p_enhancement, 'g^')


ax.set_title(r'Momentum-Scaled Enhancement for $ggg \rightarrow \Lambda + X$')
ax.set_ylabel(r'$\frac{ggg}{q\bar{q}}$ Enhancement')
ax.set_xlabel(r'Scaled $\Lambda$ Momentum')

ax.yaxis.grid(True)

ax.legend([r'$\Upsilon(1S)$ MC', r'$\Upsilon(1S)$ Data (CLEO)'], fontsize=10)


plt.ylim((-0.5, 4))
plt.hlines(1, 0, 1, color='k')
plt.tight_layout()
plt.savefig('Istogrammi/lambda_momentum_scaled.png')
plt.show()


# In[ ]:




