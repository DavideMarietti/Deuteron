#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Reads and prints the Tree structure and what it is containted.

import ROOT


filename_u = '../U1S/anti-p_U1S.root'
my_file_u = ROOT.TFile.Open(filename_u)
my_file_u.ls()
my_tree_u = my_file_u.Get('anti-p_U1S_tree')
my_tree_u.Print()

filename_c = '../Continuum/anti-p_continuum.root'
my_file_c = ROOT.TFile.Open(filename_c)
my_file_c.ls()
my_tree_c = my_file_c.Get('anti-p_continuum_tree')
my_tree_c.Print()


# In[2]:


#Translates the Tree into DataFrame, seclects the CMS_p variable and calculates the scaled momentum (CMS_p/E_beam)
#both for continuum and U1S events. Finally gathers these values in 10 bins and calculate the enhancement.

import root_pandas
import numpy as np


df_u = root_pandas.read_root(filename_u, key='anti-p_U1S_tree')
p_u = df_u['CMS_p'].to_numpy()
p_u_scaled= p_u / (9.46 * 0.5) #U1S: CMS_E = 9.460 GeV

df_c = root_pandas.read_root(filename_c, key='anti-p_continuum_tree')
p_c = df_c['CMS_p'].to_numpy()
p_c_scaled= p_c / (10.58 * 0.5) #Continuum: CMS_E = 10.580 GeV

bins = np.linspace(0, 1, 11) #10 bins
p_c_binned = np.histogram(p_c_scaled, bins)
p_u_binned = np.histogram(p_u_scaled, bins)

p_enhancement = ((p_u_binned[0]) / (p_c_binned[0]) * (98021/100000) - (0.108)) / 0.817
p_enhancement = np.nan_to_num(p_enhancement)

data_enhancement = np.array([2.07285032458065, 2.58505854729056, 2.03061647304369, 1.50952256717975, 
                             0.94952502263624, 0.611746289958718, 0.368425898927273, 0.191771151455625,
                             0.0817667009407472, 0.25510658215804])

data_error = np.array([2.20441988950276, 2.71270718232044, 2.12707182320442, 1.55801104972376, 1, 
                       0.701657458563536, 0.464088397790056, 0.292817679558011, 0.187845303867404, 
                       0.513812154696133])

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


ax.set_title(r'Momentum-Scaled Enhancement for $ggg \rightarrow \bar{p} + X$')
ax.set_ylabel(r'$\frac{ggg}{q\bar{q}}$ Enhancement')
ax.set_xlabel(r'Scaled $\bar{p}$ Momentum')


ax.yaxis.grid(True)

ax.legend([r'$\Upsilon(1S)$ MC', r'$\Upsilon(1S)$ Data (CLEO)'], fontsize=10)


plt.ylim((-0.5, 4))
plt.hlines(1, 0, 1, color='k')
plt.tight_layout()
#plt.savefig('Istogrammi/anti-p_momentum_scaled.png')


# In[ ]:





# In[ ]:




