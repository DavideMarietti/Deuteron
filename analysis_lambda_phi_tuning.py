#!/usr/bin/env python
# coding: utf-8

# In[7]:


#Notebook for TUNING of Pythia.

###################################################################################################################
#I have started generating 100.000 evets both for U1S and continuum (97716 after fragmentation) in order to test 
#the MC.

import ROOT
import root_pandas
import matplotlib.pyplot as plt
import numpy as np
from numpy import inf


filename_lc = 'lambda_continuum_tuning.root'  #Continuum data sample for lambda
filename_lu = 'lambda_U1S_tuning.root'        #U1S data sample sample for lambda

filename_pc = 'phi_continuum_tuning.root'     #Continuum data sample for phi
filename_pu = 'phi_U1S_tuning.root'           #U1S data sample sample for phi


###################################################################################################################
#Reads and prints the Tree structur and what it is containted.

my_file_lc = ROOT.TFile.Open(filename_lc)
my_file_lc.ls()
my_tree_lc = my_file_lc.Get('lambda_continuum_tuning_tree')
my_tree_lc.Print()

my_file_lu = ROOT.TFile.Open(filename_lu)
my_file_lu.ls()
my_tree_lu = my_file_lu.Get('lambda_U1S_tuning_tree')
my_tree_lu.Print()


my_file_pc = ROOT.TFile.Open(filename_pc)
my_file_pc.ls()
my_tree_pc = my_file_pc.Get('phi_continuum_tuning_tree')
my_tree_pc.Print()

my_file_pu = ROOT.TFile.Open(filename_pu)
my_file_pu.ls()
my_tree_pu = my_file_pu.Get('phi_U1S_tuning_tree')
my_tree_pu.Print()


# In[8]:


# LAMBDA CONTINUUM MOMENTUM DISTRIBUTION (the same could be easily done for phi)

###################################################################################################################
#Translates the Tree into DataFrame: each row is a lambda candidate, each column is a candidate's proper variable

df_lc = root_pandas.read_root(filename_lc, key='lambda_continuum_tuning_tree')


###################################################################################################################
#Plot (anti-)p momentum distribution in the CMS.

df_lc.hist('CMS_p', bins=100, range=(0., 5.))

tot_candidates = len(df_lc.index)

plt.xlabel('p [GeV]')
plt.ylabel('Number of candidates')
plt.title(r'$p$ CMS Momentum Distribution (Continuum)')
plt.text(2.4, 165, 'Number evts generated = 97716')
plt.text(2.43, 151, 'Tot. number candidates =')
plt.text(5.2, 151, tot_candidates)
plt.axis([-1, 6, 0, 200])


###################################################################################################################
#Show the result

#plt.savefig('Istogrammi/p_continuum_momentum_dist_tuning.jpg', format='jpg', dpi=900)
#plt.savefig('Istogrammi/p_ccut_continuum_momentum_dist_tuning.jpg', format='jpg', dpi=900)
plt.show()

df_lc.tail(5)


# In[9]:


# Lambda U1S MOMENTUM DISTRIBUTION (the same could be easily done for phi)

###################################################################################################################
#Translates the Tree into DataFrame: each row is a lambda candidate, each column is a candidate's proper variable

df_lu = root_pandas.read_root(filename_lu, key='lambda_U1S_tuning_tree')


###################################################################################################################
#PLOT (anti-)p momentum distribution in the CMS.

df_lu.hist('CMS_p', bins=100, range=(0., 5.))

tot_candidates = len(df_lu.index)

plt.xlabel('p [GeV]')
plt.ylabel('Number of candidates')
plt.title(r'$p$ CMS Momentum Distribution (Continuum)')
plt.text(2.4, 350, 'Number evts generated = 10.000')
plt.text(2.43, 310, 'Tot. number candidates =')
plt.text(5.2, 310, tot_candidates)
plt.axis([-1, 6, 0, 600])


###################################################################################################################
#Show the result

#plt.savefig('Istogrammi/p_continuum_momentum_dist_tuning.jpg', format='jpg', dpi=900)
#plt.savefig('Istogrammi/p_ccut_continuum_momentum_dist_tuning.jpg', format='jpg', dpi=900)
plt.show()

df_lu.tail(5)


# In[10]:


# Lambda ENHANCEMENT - SCALED MOMENTUM

###################################################################################################################
#Calculating the enhancement as function of the SCALED momentum for LAMBDA from MC data

df_lc = root_pandas.read_root(filename_lc, key='lambda_continuum_tuning_tree')

df_lu = root_pandas.read_root(filename_lu, key='lambda_U1S_tuning_tree')


p_lc = df_lc['CMS_p'].to_numpy()
p_lc_scaled= p_lc / (10.58 * 0.5) #Continuum: CMS_E = 10.580 GeV

p_lu = df_lu['CMS_p'].to_numpy()
p_lu_scaled= p_lu / (9.46 * 0.5) #U1S: CMS_E = 9.460 GeV


bins = np.linspace(0, 1, 11) #10 bins
p_lc_binned = np.histogram(p_lc_scaled, bins)
p_lu_binned = np.histogram(p_lu_scaled, bins)

p_enhancement = (((p_lu_binned[0])/(p_lc_binned[0])) * (97716/100000) - (0.108)) / 0.817 #Correction
p_enhancement[p_enhancement == inf] = 0         #Introduced after the above correction for the bins at high scaled 
p_enhancement[p_enhancement < 0] = 0            #momentum, where there could be few particles ad values likes 0/0,
p_enhancement = np.nan_to_num(p_enhancement)    #something/0 or -something can occur.

p_error = p_enhancement * np.sqrt((np.sqrt(p_lu_binned[0])/p_lu_binned[0])*(np.sqrt(p_lu_binned[0])/p_lu_binned[0])
                                  + (np.sqrt(p_lc_binned[0])/p_lc_binned[0])*(np.sqrt(p_lc_binned[0])/p_lc_binned[0])) #error


###################################################################################################################
#Arrays containing enhancement data and errors for protons, extracted from the PAPER's graphics

data_enhancement = np.array([2.37096774193548, 3.21505376344086, 2.94623655913978, 2.04838709677419,
                             1.45698924731183, 1.09677419354839, 0.827956989247312, 0, 0, 0])

data_error = np.array([2.9367572894325, 3.53031708306765, 3.0538954263521, 2.12596371224857, 1.54741754101212,
                       1.23240681699612, 0.997985624021796, 0.84418294591618, 0.265694742333778, 0.316083125615906])

error = data_error - data_enhancement


###################################################################################################################
#Changing and defining variables for the momentum integrated section were I will consider single once multicandidate events

df_lc = df_lc.loc[df_lc['__candidate__'] == 0]
df_lu = df_lu.loc[df_lu['__candidate__'] == 0]

tot_lambda_c = len(df_lc.index)
tot_lambda_u = len(df_lu.index)


###################################################################################################################
#Create a plot of the enhancement as function of the momentum scaled both from data and MC

fig, ax = plt.subplots()
new_bins = np.array(bins[:-1]) + 0.5 * (bins[1] - bins[0])


#Errorbar graphics for data and MC are defined

ax.errorbar(new_bins, data_enhancement, yerr=error, fmt='bs', elinewidth=0.7, barsabove=True)
ax.errorbar(new_bins, p_enhancement, yerr=p_error, fmt='g^', elinewidth=0.6)


#Decorations

ax.set_title(r'Momentum-Scaled Enhancement for $ggg \rightarrow \Lambda + X$')
ax.set_ylabel(r'$\frac{ggg}{q\bar{q}}$ Enhancement')
ax.set_xlabel(r'Scaled $\Lambda$ Momentum')
ax.legend([r'$\Upsilon(1S)$ MC', r'$\Upsilon(1S)$ Data (CLEO)'], fontsize=10)
ax.yaxis.grid(True)


#Plot layout

plt.ylim((-0.5, 4))
plt.hlines(1, 0, 1, color='k')
plt.tight_layout()


###################################################################################################################
#Show the result

#plt.savefig('Istogrammi/tuning_lambda_scaled.jpg', format='jpg', dpi=900)
plt.show()

print('ggg:    ', p_lu_binned[0])
print('q qbar: ', p_lc_binned[0])

print('Enhancement MC :', p_enhancement)
print('Error_enhan MC :', p_error)

print('Enhancement Data :', data_enhancement)


# In[11]:


#PHI ENHANCEMENT - SCALED MOMENTUM

###################################################################################################################
#Calculating the enhancement as function of the SCALED momentum for PHI from MC data

df_pc = root_pandas.read_root(filename_pc, key='phi_continuum_tuning_tree')

df_pu = root_pandas.read_root(filename_pu, key='phi_U1S_tuning_tree')


p_pc = df_pc['CMS_p'].to_numpy()
p_pc_scaled= p_pc / (10.58 * 0.5) #Continuum: CMS_E = 10.580 GeV

p_pu = df_pu['CMS_p'].to_numpy()
p_pu_scaled= p_pu / (9.46 * 0.5) #U1S: CMS_E = 9.460 GeV


bins = np.linspace(0, 1, 11) #10 bins
p_pc_binned = np.histogram(p_pc_scaled, bins)
p_pu_binned = np.histogram(p_pu_scaled, bins)

p_enhancement = (((p_pu_binned[0])/(p_pc_binned[0])) * (97716/100000) - (0.108)) / 0.817 #Correction
p_enhancement[p_enhancement == inf] = 0
p_enhancement[p_enhancement < 0] = 0
p_enhancement = np.nan_to_num(p_enhancement)

p_error = p_enhancement * np.sqrt((np.sqrt(p_pu_binned[0])/p_pu_binned[0])*(np.sqrt(p_pu_binned[0])/p_pu_binned[0])
                                  + (np.sqrt(p_pc_binned[0])/p_pc_binned[0])*(np.sqrt(p_pc_binned[0])/p_pc_binned[0])) #error


###################################################################################################################
#Arrays containing enhancement data and errors for anti-protons, extracted from the PAPER's graphics

data_enhancement = np.array([2.04899481146392, 2.02857308031391, 1.51888251362171, 1.06834703935483,
                             0.816739171279502, 0.764044456648986, 0.662961945917822, 0.674798745501583, 0, 0])

data_error = np.array([2.97849462365591, 2.11827956989247, 1.55913978494624, 1.11827956989247, 0.881720430107527,
                       0.844086021505377, 0.790322580645162, 0.849462365591398, 0.763440860215054,
                       0.333333333333333])

error = data_error - data_enhancement


###################################################################################################################
#Changing and defining variables for the momentum integrated section were I will consider single once multicandidate events

df_pc = df_pc.loc[df_pc['__candidate__'] == 0]
df_pu = df_pu.loc[df_pu['__candidate__'] == 0]

tot_phi_c = len(df_pc.index)
tot_phi_u = len(df_pu.index)


###################################################################################################################
#Create a plot of the enhancement as function of the momentum scaled both from data and MC

fig, ax = plt.subplots()
new_bins = np.array(bins[:-1]) + 0.5 * (bins[1] - bins[0])


#Errorbar graphics for data and MC are defined

ax.errorbar(new_bins, data_enhancement, yerr=error, fmt='bs', elinewidth=0.7, barsabove=True)
ax.errorbar(new_bins, p_enhancement, yerr=p_error, fmt='g^', elinewidth=0.6)


#Decorations

ax.set_title(r'Momentum-Scaled Enhancement for $ggg \rightarrow \Phi + X$')
ax.set_ylabel(r'$\frac{ggg}{q\bar{q}}$ Enhancement')
ax.set_xlabel(r'Scaled $\Phi$ Momentum')
ax.legend([r'$\Upsilon(1S)$ MC', r'$\Upsilon(1S)$ Data (CLEO)'], fontsize=10)
ax.yaxis.grid(True)


#Plot layout

plt.ylim((-0.5, 4))
plt.hlines(1, 0, 1, color='k')
plt.tight_layout()


###################################################################################################################
#Show the result

#plt.savefig('Istogrammi/tuning_phi_scaled.jpg', format='jpg', dpi=900)
plt.show()

print('ggg:    ', p_pu_binned[0])
print('q qbar: ', p_pc_binned[0])

print('Enhancement MC :', p_enhancement)
print('Error_enhan MC :', p_error)

print('Enhancement Data :', data_enhancement)


# In[1]:


# ENHANCEMENT - INTEGRATED

###################################################################################################################
#Calculating the enhancement as function of the INTEGRATED momentum for lamda and phi from MC data using the
#variables earlier defined.

lam = ((tot_lambda_u/tot_lambda_c) * (97716/100000) - 0.108) / 0.817               #Enhancement for protons
phi = ((tot_phi_u/tot_phi_c) * (97716/100000) - 0.108) / 0.817  #Enhancement for anti-protons

MC = np.array([lam, phi])

tlu = np.sqrt(tot_lambda_u)/tot_lambda_u #Defining variable to make the error definition more readeable
tlc = np.sqrt(tot_lambda_c)/tot_lambda_c
tpu = np.sqrt(tot_phi_u)/tot_phi_u
tpc = np.sqrt(tot_phi_c)/tot_phi_c
err_MC = np.array([lam*np.sqrt(tlu*tlu + tlc*tlc), phi*np.sqrt(tpu*tpu + tpc*tpc)]) #Error


#Arrays containing enhancement data and errors for lambda and phi, extracted from the PAPER's graphics

data = np.array([2.668, 1.423])
err_stat_data = np.array([0.027, 0.051])
err_sys_data  = np.array([0.051, 0.065])
err_data = np.sqrt(err_stat_data*err_stat_data + err_sys_data*err_sys_data)


#Errorbar graphics for data and MC are defined

labels = [r'$\Lambda$', r'$\Phi$',]
x_pos = np.arange(len(labels))

fig, ax = plt.subplots()
particles = np.array([0, 1])
ax.errorbar(particles, data, yerr=err_data, fmt='r.', elinewidth=0.6)
ax.errorbar(particles, MC, yerr=err_MC, fmt='g^', elinewidth=0.6)

                 
#Decorations

ax.set_ylabel(r'$\frac{ggg}{q\bar{q}}$ Enhancement')
ax.set_xticks(x_pos)
ax.set_xticklabels(labels, fontsize=15)
ax.set_title('Momentum-Integrated Enhancement for $ggg$ events')
ax.legend([r'$\Upsilon(1S)$ MC', r'$\Upsilon(1S)$ Data (CLEO)'], fontsize=10)
ax.yaxis.grid(True)


#Plot layout

plt.ylim((0, 3))
plt.hlines(1, -1, 2, color='k')
plt.tight_layout()


###################################################################################################################
#Show the result

#plt.savefig('Istogrammi/tuning_integr_enhancement_lp.jpg', format='jpg', dpi=900)
plt.show()

print('MC: [', tot_lambda_u , '/', tot_lambda_c, '=', MC[0], ' ', 
      tot_phi_u , '/', tot_phi_c, '=', MC[1], ']') #Apart the correction
print('MC_err: [', err_MC[0], ' ', err_MC[1], ']')

print('Data:', data)


# 

# In[ ]:




