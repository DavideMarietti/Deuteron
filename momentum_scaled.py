#!/usr/bin/env python
# coding: utf-8

# In[6]:


#anti-p 

import ROOT
import root_pandas
import numpy as np


filename_u = '../U1S/anti-p_U1S.root'
my_file_u = ROOT.TFile.Open(filename_u)
my_tree_u = my_file_u.Get('anti-p_U1S_tree')

filename_c = '../Continuum/anti-p_continuum.root'
my_file_c = ROOT.TFile.Open(filename_c)
my_tree_c = my_file_c.Get('anti-p_continuum_tree')


df_u = root_pandas.read_root(filename_u, key='anti-p_U1S_tree')
p_u = df_u['CMS_p'].to_numpy()
p_u_scaled = p_u / (9.46 * 0.5) #U1S: CMS_E = 9.460 GeV

df_c = root_pandas.read_root(filename_c, key='anti-p_continuum_tree')
p_c = df_c['CMS_p'].to_numpy()
p_c_scaled = p_c / (10.58 * 0.5) #Continuum: CMS_E = 10.580 GeV

bins = np.linspace(0, 1, 11) #10 bins
p_c_binned = np.histogram(p_c_scaled, bins)
p_u_binned = np.histogram(p_u_scaled, bins)

anti_p_enhance = (((p_u_binned[0]) / (p_c_binned[0])) * (98021/100000) - (0.108)) / 0.817
anti_p_enhance[anti_p_enhance < 0] = 0
anti_p_enhance = np.nan_to_num(anti_p_enhance)


#Visual check
print(anti_p_enhance)
print(p_u_binned[0])
print(p_c_binned[0])


# In[7]:


#p

from numpy import inf

filename_a = '../U1S/p_U1S.root'
my_file_a = ROOT.TFile.Open(filename_a)
my_tree_a = my_file_a.Get('p_U1S_tree')

filename_b = '../Continuum/p_continuum.root'
my_file_b = ROOT.TFile.Open(filename_b)
my_tree_b = my_file_b.Get('p_continuum_tree')


df_a = root_pandas.read_root(filename_a, key='p_U1S_tree')
p_a = df_a['CMS_p'].to_numpy()
p_a_scaled = p_a / (9.46 * 0.5) #U1S: CMS_E = 9.460 GeV

df_b = root_pandas.read_root(filename_b, key='p_continuum_tree')
p_b = df_b['CMS_p'].to_numpy()
p_b_scaled = p_b / (10.58 * 0.5) #Continuum: CMS_E = 10.580 GeV

bins = np.linspace(0, 1, 11) #10 bins
p_b_binned = np.histogram(p_b_scaled, bins)
p_a_binned = np.histogram(p_a_scaled, bins)

p_enhance = (((p_a_binned[0])/(p_b_binned[0])) * (98021/100000) - (0.108)) / 0.817
p_enhance[p_enhance == inf] = 0
p_enhance[p_enhance < 0] = 0
p_enhance = np.nan_to_num(p_enhance)

#Visual check
print(p_enhance)
print(p_a_binned[0])
print(p_b_binned[0])


# In[8]:


#lambda

filename_d = '../U1S/lambda_U1S.root'
my_file_d = ROOT.TFile.Open(filename_d)
my_tree_d = my_file_d.Get('lambda_U1S_tree')

filename_f = '../Continuum/lambda_continuum.root'
my_file_f = ROOT.TFile.Open(filename_f)
my_tree_f = my_file_f.Get('lambda_continuum_tree')

df_d = root_pandas.read_root(filename_d, key='lambda_U1S_tree')
p_d = df_d['CMS_p'].to_numpy()
p_d_scaled = p_d / (9.46 * 0.5) #U1S: CMS_E = 9.460 GeV

df_f = root_pandas.read_root(filename_f, key='lambda_continuum_tree')
p_f = df_f['CMS_p'].to_numpy()
p_f_scaled = p_f / (10.58 * 0.5) #Continuum: CMS_E = 10.580 GeV

bins = np.linspace(0, 1, 11) #10 bins
p_f_binned = np.histogram(p_f_scaled, bins)
p_d_binned = np.histogram(p_d_scaled, bins)

lambda_enhance = (((p_d_binned[0]) / (p_f_binned[0])) * (98021/100000) - (0.108)) / 0.817
lambda_enhance[lambda_enhance < 0] = 0

#Visual check
print(lambda_enhance)
print(p_d_binned[0])
print(p_f_binned[0])


# In[9]:


#phi

filename_q = '../U1S/phi_U1S.root'
my_file_q = ROOT.TFile.Open(filename_q)
my_tree_q = my_file_q.Get('phi_U1S_tree')

filename_w = '../Continuum/phi_continuum.root'
my_file_w = ROOT.TFile.Open(filename_w)
ytree_w = my_file_w.Get('phi_continuum_tree')


df_q = root_pandas.read_root(filename_q, key='phi_U1S_tree')
p_q = df_q['CMS_p'].to_numpy()
p_q_scaled = p_q / (9.46 * 0.5) #U1S: CMS_E = 9.460 GeV

df_w = root_pandas.read_root(filename_w, key='phi_continuum_tree')
p_w = df_w['CMS_p'].to_numpy()
p_w_scaled = p_w / (10.58 * 0.5) #Continuum: CMS_E = 10.580 GeV

bins = np.linspace(0, 1, 11) #10 bins
p_w_binned = np.histogram(p_w_scaled, bins)
p_q_binned = np.histogram(p_q_scaled, bins)

phi_enhance = (((p_q_binned[0])/(p_w_binned[0])) * (98021/100000) - (0.108)) / 0.817
phi_enhance[phi_enhance < 0] = 0.
phi_enhance = np.nan_to_num(phi_enhance)


#Visual check
print(phi_enhance)
print(p_u_binned[0])
print(p_c_binned[0])


# In[10]:


#Plot the momentum-scaled enhancement. Multiple candidates in the same event are considered.

get_ipython().run_line_magic('matplotlib', 'ipympl')
import matplotlib.pyplot as plt


fig, ax = plt.subplots()
new_bins = np.array(bins[:-1]) + 0.5 * (bins[1] - bins[0])
ax.plot(new_bins, lambda_enhance, 'bd')
ax.plot(new_bins, p_enhance, 'gx')
ax.plot(new_bins, anti_p_enhance, 'r.')
ax.plot(new_bins, phi_enhance, 'y*')




ax.set_title(r'Momentum-Scaled Enhancement for $\Upsilon (1S) \rightarrow ggg \rightarrow \Lambda, p, \bar{p}, \Phi$ + X')
ax.set_ylabel(r'$\frac{ggg}{q\bar{q}}$ Enhancement')
ax.set_xlabel(r'Scaled Momentum ($\frac{p_{CMS}}{E_{beam}}$)')

ax.yaxis.grid(True)

ax.legend([r'$\Lambda$', r'$p$', r'$\bar{p}$', r'$\Phi$'], fontsize=10)

plt.ylim((-0.5, 3))
plt.hlines(1, 0, 1, color='k')
plt.tight_layout()
plt.savefig('Istogrammi/momentum_scaled.png')
plt.show()


# In[ ]:




