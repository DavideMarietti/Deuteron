#!/usr/bin/env python
# coding: utf-8

# In[1]:


##############################################################################################################
#I have generated 50.000.000 Upsilon(1S/2S/3S) events and selected the p-n, p-p and n-n pairs that can be
#further analyzed to identify the anti-deuteron candidates with the Cross Section Model.
#NOTE: If not explicity specify, p-n pairs from lambda decay are not considered in this analysis.
#NOTE: An alternative to exclude weak decay products is to apply the cut: r = |r_p - r_n| < 3 fm. (Dal-Raklev)
#NOTE: Taking into account the smallness of the cross section probabilities, the contributions of triple and
#      higher deuteron formation per event have been neglected here.
##############################################################################################################

#Translates the Tree into DataFrame
import ROOT
import root_pandas
import numpy as np


#Show p-n pairs selected events
df_pn = root_pandas.read_root('anti_pn_continuum_recap.root', key='anti_pn_continuum_recap_tree')
df_pn.head(5)


# In[2]:


#Show p-p pairs selected events
df_pp = root_pandas.read_root('anti_pp_continuum_recap.root', key='anti_pp_continuum_recap_tree')
df_pp.head(5)


# In[3]:


#Show n-n pairs selected events
df_nn = root_pandas.read_root('anti_nn_continuum_recap.root', key='anti_nn_continuum_recap_tree')
df_nn.head(5)


# In[4]:


#Concatenate p-n, p-p, n-n DataFrames arrays into a single multi-index DataFrame array containing all possible
#pair combination without duplicates.

import pandas as pd


#Store all pairs information into a single numpy array: each row is a selected p-n, p-p or n-n pair candidate
pairs = np.array(df_pn)
pairs = np.append(pairs, df_pp, axis=0)
pairs = np.append(pairs, df_nn, axis=0)

#Sort rows in crescent value of the '__event__' column elements
pairs = pairs[np.argsort(pairs[:, 2])]

#Define an array that contains the event's values in which have been found at least one pair
a = np.array(df_pn.loc[:,'__event__'])
b = np.array(df_pp.loc[:,'__event__'])
c = np.array(df_nn.loc[:,'__event__'])

event_array = np.unique(np.concatenate((a,b,c),0)) #merge them without duplicates


#Define the second level of the multi-index array
k = [1] #Initialize the array with one single element
for n in range(0, len(event_array)):
    i = len(pairs[np.where(pairs[:,2] == event_array[n])])
    h = np.arange(0,i)
    k = np.append(k,h)

k = np.delete(k, 0) #Delete the first element of the array, defined for the initialization
arrays = [pairs[:,2].astype(int), np.array(k)] #Define the first and second levels of the multi-index array
pairs = np.delete(pairs, [0,1,2,3,4,5,11,12,17,18], axis=1) #Delete negligible columns to improve readability


#Define the multi-index DataFrame array
df = pd.DataFrame(pairs, index=arrays)

df.rename(columns={df.columns[0]: "PDG_N1"}, inplace = True)
df.rename(columns={df.columns[1]: "PDG_N2"}, inplace = True)
df.rename(columns={df.columns[2]:  "px_N1"}, inplace = True)
df.rename(columns={df.columns[3]:  "py_N1"}, inplace = True)
df.rename(columns={df.columns[4]:  "pz_N1"}, inplace = True)
df.rename(columns={df.columns[5]:   "E_N1"}, inplace = True)
df.rename(columns={df.columns[6]:  "px_N2"}, inplace = True)
df.rename(columns={df.columns[7]:  "py_N2"}, inplace = True)
df.rename(columns={df.columns[8]:  "pz_N2"}, inplace = True)
df.rename(columns={df.columns[9]:   "E_N2"}, inplace = True)

df = df.astype({"PDG_N1": int, "PDG_N2": int})


df.tail(25) #Visual check


# In[5]:


#GLOBAL VARIABLES are defined here

#Mass definitions
d_m = 1.87561294257  #+- 5.7e-10 GeV
n_m = 0.93956542052  #+- 5.4e-10 GeV
p_m = 0.93827208816  #+- 2.9e-10 Gev
pi_m = 0.13957039    #+- 1.8e-07 GeV
pi0_m = 0.1349768    #+- 5.0e-07 GeV

#Quadratic mass definitions
d_m2 = pow(d_m, 2)
n_m2 = pow(n_m, 2)
p_m2 = pow(p_m, 2)
pi_m2 = pow(pi_m, 2)
pi0_m2 = pow(pi0_m, 2)


# In[6]:


##############################################################################################################
#Paper Dal-Raklev: definition of the functions that calculate the cross section for the deuteron production
#(starting from p-n, p-p and n-n pairs) as function of k.
##############################################################################################################


# In[7]:


##############################################################################################################
#Function definition: Calculate the cross section of the process n + p -> d + gamma as function of k.
#INPUT:   k = two times the p-n relative momentum in the CMS frame of the pair
#OUTPUT:  sigma = cross section value at the given k
##############################################################################################################

def sez_urto_pn_rad(k): #k = k / 1 GeV
    
    #Best fit values to the parameter used to fit the cross section
    a = np.array([2.30346, -9.366346e01, 2.565390e03, -2.5594101e04, 1.43513109e05, -5.0357289e05,
                  1.14924802e06, -1.72368391e06, 1.67934876e06, -1.01988855e06, 3.4984035e05, -5.1662760e04])
    b = np.array([-5.1885, 2.9196])

    #Calculate the cross section as function of k
    if 0 < k < 1.28:
        
        sigma = 0.0 #Initialize the cross section value in the summation        
        for n in range(0,len(a)):
            sigma = sigma + a[n] * pow(k, n-1)
        
        return sigma
        
    elif k >= 1.28:
        
        sigma = np.exp(-b[0] * k -b[1] * pow(k, 2))
        
        return sigma


##############################################################################################################
#Function definition: Calculate the cross section of the process p + n -> d + pi0 as function of k.
#INPUT:   k = two times the p-n relative momentum in the CMS frame of the pair
#OUTPUT:  sigma = cross section value at the given k
##############################################################################################################

def sez_urto_pn_pi0(k):
    
    #Best fit values to the parameter used to fit the cross section
    pn = np.array([170.0, 1.34, 1.77, 0.38, 0.096])
    
    #Calculate the threshold value of k for this process
    t = d_m2 + pi0_m2 + 2*d_m*pi0_m - p_m2 - n_m2
    k_min = np.sqrt((pow(t, 2) - 4*p_m2*n_m2) / (t + p_m2 + n_m2))

    #Calculate the cross section as function of eta = q / pi_m, q = pion momentum in the CMS frame    
    if k > k_min:
        
        #Translate k (relative momentum among p-n in their CMS) in eta (momentum of the emitted pion in the CMS)
        h = 2*pow(0.5*k, 2) + 2*np.sqrt(p_m2 + pow(0.5*k, 2))*np.sqrt(n_m2 + pow(0.5*k, 2)) - d_m2 - pi0_m2 + p_m2 + n_m2
        q = (0.25*pow(h, 2) - d_m2 * pi0_m2) / (h + d_m2 + pi0_m2)
        eta = np.sqrt(q) / pi0_m
        
        sigma = (pn[0] * pow(eta, pn[1])) / (pow(pn[2] - np.exp(pn[3] * eta), 2) + pn[4])
        
        return 0.5*sigma #Exploit the Isospin Invariance

    else:        
        return 0

    
##############################################################################################################
#Function definition: Calculate the cross section of the process p + n -> d + pi0 + pi0 as function of k.
#INPUT:   k = two times the p-n relative momentum in the CMS frame of the pair
#OUTPUT:  sigma = cross section value at the given k
##############################################################################################################

def sez_urto_pn_pi0pi0(k):
    
    #Best fit values to the parameter used to fit the cross section
    pn = np.array([2.855e06, 1.311e01, 2.961e03, 5.572e00, 1.461e06])
    
    #Calculate the threshold value of k for this process
    t = d_m2 + 4*pi0_m2 + 4*d_m*pi0_m - p_m2 - n_m2 
    k_min = np.sqrt((pow(t, 2) - 4*p_m2*n_m2) / (t + p_m2 + n_m2))
    
    #Calculate the cross section as function of k 
    if k > k_min:
        
        sigma = (pn[0] * pow(k, pn[1])) / (pow(pn[2] - np.exp(pn[3] * k), 2) + pn[4])
       
        return sigma
    
    else:        
        return 0
    

##############################################################################################################
#Function definition: Calculate the cross section of the process p + n -> d + pi+ + pi- as function of k.
#INPUT:   k = two times the p-n relative momentum in the CMS frame of the pair
#OUTPUT:  sigma = cross section value at the given k
##############################################################################################################

def sez_urto_pn_pipi(k):
    
    #Best fit values to the parameter used to fit the cross section
    pn1 = np.array([6.465e06, 1.051e01, 1.979e03, 5.363e00, 6.045e05])
    pn2 = np.array([2.549e15, 1.657e01, 2.330e07, 1.119e01, 2.868e16])
    
    #Calculate the threshold value of k for this process
    t = d_m2 + 4*pi_m2 + 4*d_m*pi_m - p_m2 - n_m2
    k_min = np.sqrt((pow(t, 2) - 4*p_m2*n_m2) / (t + p_m2 + n_m2))

    #Calculate the cross section as function of k 
    if k > k_min:
        
        sigma1 = (pn1[0] * pow(k, pn1[1])) / (pow(pn1[2] - np.exp(pn1[3] * k), 2) + pn1[4])
        sigma2 = (pn2[0] * pow(k, pn2[1])) / (pow(pn2[2] - np.exp(pn2[3] * k), 2) + pn2[4])
    
        sigma = sigma1 + sigma2
    
        return sigma
    
    else:        
        return 0


##############################################################################################################
#Function definition: Calculate the cross section of the process p + p -> d + pi+ as function of k.
#INPUT:   k = two times the p-n relative momentum in the CMS frame of the pair
#OUTPUT:  sigma = cross section value at the given k
##############################################################################################################

def sez_urto_pp_pi(k):
    
    #Best fit values to the parameter used to fit the cross section    
    pp = np.array([170, 1.34, 1.77, 0.38, 0.096])
    
    #Calculate the threshold value of k for this process
    k_min = np.sqrt(d_m2 + pi_m2 + 2*d_m*pi_m -4*p_m2)
    
    #Calculate the cross section as function of eta = q / pi_m, q = pion momentum in the CMS frame
    if k > k_min:
        
        #Translate k (relative momentum among p-n in their CMS) in eta (momentum of the emitted pion in the CMS)
        h = 4*p_m2 - d_m2 - pi_m2 + pow(k, 2)
        q = (0.25*pow(h, 2) - d_m2 * pi0_m2) / (h + d_m2 + pi0_m2)
        eta = np.sqrt(q) / pi_m
        
        sigma = (pp[0] * pow(eta, pp[1])) / (pow(pp[2] - np.exp(pp[3] * eta), 2) + pp[4])
        
        return sigma
    
    else:
        return 0


##############################################################################################################
#Function definition: Calculate the cross section of the process p + p -> d + pi+ + pi0 as function of k.
#INPUT:   k = two times the p-n relative momentum in the CMS frame of the pair
#OUTPUT:  sigma = cross section value at the given k
##############################################################################################################

def sez_urto_pp_pipi0(k):
    
    #Best fit values to the parameter used to fit the cross section
    pp = np.array([5.099e15, 1.656e01, 2.333e07, 1.133e01, 2.868e16])

    #Calculate the threshold value of k for this process
    k_min = np.sqrt(d_m2 + pi0_m2 + pi_m2 + 2*d_m*pi0_m + 2*d_m*pi_m + 2*pi_m*pi0_m - 4*p_m2)
    
    #Calculate the cross section as function of k 
    if k > k_min:
        
        sigma = (pp[0] * pow(k, pp[1])) / (pow(pp[2] - np.exp(pp[3] * k), 2) + pp[4])
        
        return sigma
    
    else:
        return 0

    
##############################################################################################################
#Function definition: Calculate the cross section of the process n + n -> d + pi- as function of k.
#INPUT:   k = two times the p-n relative momentum in the CMS frame of the pair
#OUTPUT:  sigma = cross section value at the given k
##############################################################################################################

def sez_urto_nn_pi(k):
    
    #Best fit values to the parameter used to fit the cross section   
    pp = np.array([170, 1.34, 1.77, 0.38, 0.096])
    
    #Calculate the threshold value of k for this process
    k_min = np.sqrt(d_m2 + pi_m2 + 2*d_m*pi_m -4*n_m2)
    
    #Calculate the cross section as function of eta = q / pi_m, q = pion momentum in the CMS frame
    if k > k_min:
        
        #Translate k (relative momentum among p-n in their CMS) in eta (momentum of the emitted pion in the CMS)
        h = 4*n_m2 - d_m2 - pi_m2 + pow(k, 2)
        q = (0.25*pow(h, 2) - d_m2 * pi0_m2) / (h + d_m2 + pi0_m2)
        eta = np.sqrt(q) / pi_m
        
        sigma = (pp[0] * pow(eta, pp[1])) / (pow(pp[2] - np.exp(pp[3] * eta), 2) + pp[4])
        
        return sigma
    
    else:
        return 0


##############################################################################################################
#Function definition: Calculate the cross section of the process n + n -> d + pi- + pi0 as function of k.
#INPUT:   k = two times the p-n relative momentum in the CMS frame of the pair
#OUTPUT:  sigma = cross section value at the given k
##############################################################################################################

def sez_urto_nn_pipi0(k):
    
    #Best fit values to the parameter used to fit the cross section
    nn = np.array([5.099e15, 1.656e01, 2.333e07, 1.133e01, 2.868e16])

    #Calculate the threshold value of k for this process
    k_min = np.sqrt(d_m2 + pi0_m2 + pi_m2 + 2*d_m*pi0_m + 2*d_m*pi_m + 2*pi_m*pi0_m - 4*n_m2)
    
    #Calculate the cross section as function of k 
    if k > k_min:
        
        sigma = (nn[0] * pow(k, nn[1])) / (pow(nn[2] - np.exp(nn[3] * k), 2) + nn[4])
        
        return sigma
    
    else:
        return 0


# In[8]:


#VISUAL CHECK: Plot the cross section trends respect to k for the considered processes.

import matplotlib.pyplot as plt

y1 = []
y2 = []
y3 = []
y4 = []
z1 = []
z2 = []
z3 = []
z4 = []

x  = np.linspace(0, 4, 10000)

for i in range(0,10000):
    y1.append(sez_urto_pn_rad(x[i]))
for i in range(0,10000):
    y2.append(sez_urto_pn_pi0(x[i]))
for i in range(0,10000):
    y3.append(sez_urto_pn_pi0pi0(x[i]))
for i in range(0,10000):
    y4.append(sez_urto_pn_pipi(x[i]))

for i in range(0,10000):
    z1.append(sez_urto_pp_pi(x[i]))
for i in range(0,10000):
    z2.append(sez_urto_pp_pipi0(x[i]))
for i in range(0,10000):
    z3.append(sez_urto_nn_pi(x[i]))
for i in range(0,10000):
    z4.append(sez_urto_nn_pipi0(x[i]))


plt.figure(figsize=(18, 6))

#Plot SX figure (Fig.4 in the Dal-Raklev paper)
plt.subplot(121)
plt.plot(x,y1, 'k-.')
plt.plot(x,y2, 'g--')
plt.plot(x,y3, 'b:')
plt.plot(x,y4, 'r')

plt.xlabel(r'$k$ [GeV]')
plt.ylabel(r'$\sigma$ $[\mu b]$')
plt.title(r'Cross sections for $n p \rightarrow d X$ processes')
plt.legend([r'$pn\rightarrow d\gamma$', r'$pn\rightarrow d\pi^0$', r'$pn\rightarrow d\pi^0\pi^0$', r'$pn\rightarrow d\pi^+\pi^-$'],
           loc='upper left', fontsize=13)
plt.xlim([0.1,4])
plt.ylim([0.99,5000])
plt.yscale('log')
plt.xscale('log')

#Plot DX figure (Fig.5 in the Dal-Raklev paper)
plt.subplot(122)
plt.plot(x,z1, 'r')
plt.plot(x,z2, 'y')
plt.plot(x,z3, 'b-.')
plt.plot(x,z4, 'k-.')

plt.xlabel(r'$k$ [GeV]')
plt.ylabel(r'$\sigma$ $[\mu b]$')
plt.title(r'Cross sections for $pp \rightarrow dX$ and $n n \rightarrow d X$ processes')
plt.legend([r'$pp\rightarrow d\pi^+$', r'$pp\rightarrow d\pi^+\pi^0$', r'$nn\rightarrow d\pi^-$', r'$nn\rightarrow d\pi^-\pi^0$'],
           loc='upper left', fontsize=13)

plt.xlim([0.1,4])
plt.ylim([0.99,5000])
plt.yscale('log')
plt.xscale('log')


# In[9]:


#BASIC FUNCTIONS are defined here.

import numpy as np


##############################################################################################################  
#Function definition: Extracs tri-vector from a four-vector.
#INPUT:   quadrivector = four-vector
#OUTPUT:  tri-vector = tri-vector in form of numpy array
##############################################################################################################

def trivector(quadrivector):
    
    quadrivector = np.array(quadrivector)
    tri_vector =   np.array([quadrivector[1], quadrivector[2], quadrivector[3]])

    return tri_vector


##############################################################################################################
#Function definition: Calculates the QUADRATIC norm in minkowsky space, metric = diag(+,-,-,-).
#INPUT:   quadrivector = four-vector
#OUTPUT:  minkowsky_norm = the quadratic norm of the four-vector in minkowsky space, i.e. a number
##############################################################################################################

def normM(quadrivector):
    
    quadrivector = np.array(quadrivector)
    minkowsky_norm = 0
    
    for i in quadrivector:
        minkowsky_norm -= i*i  
        
    minkowsky_norm += 2*quadrivector[0] * quadrivector[0]
        
    return minkowsky_norm


# In[10]:


import numpy as np
from numpy import linalg as la


##############################################################################################################
#Function deinition: Given a certain four-vector, it operates the LT along a certain direction defined by the
#                    Lorentz's beta and returns the trasformed four-vector.
#INPUT:   quadrivector = four-vector in phase space, i.e. quadrivector = [E, px, py, pz]
#         bata = Lorentz's beta parameter, i.e. b = [bx, by, bz]
#OUTPUT:  quadrivector_transformed = LT trasformed four-vector, i.e. quadrivector_tranf = [E', px', py', pz']
##############################################################################################################

def boost(quadrivector, beta):
    
    quadrivector = np.array(quadrivector)
    beta = np.array(beta)

    #Calculates factors for the LT
    bx, by, bz = beta                    #Extracts beta's components
    beta_norm =  la.norm(beta)           #Lorentz's beta
    beta_norm2 = pow(beta_norm, 2)
        
    gamma = pow(1 - beta_norm2, -0.5)    #Lorentz's gamma
        
    h = (gamma - 1) / beta_norm2         #Factor employed in the trasformation
    
    
    #Matrix of the generic LT
    boost_matrix = np.array([-gamma * np.insert(beta, 0, -1),
                            [-gamma * bx, 1 + h * pow(bx,2), h * bx * by, h * bx * bz],
                            [-gamma * by, h * bx * by, 1 + h * pow(by,2), h * by * bz],
                            [-gamma * bz, h * bx * bz, h * bz * by, 1 + h * pow(bz,2)]])
    
    
    #Calculates the transformed quadrivector
    quadrivector_transformed = boost_matrix.dot(quadrivector)
    
    return quadrivector_transformed


# In[11]:


import numpy as np
from numpy.random import default_rng

rng = default_rng()


##############################################################################################################
#Function deinition: Given the E_cms and the mass of the pions (deuteron mass is already known), it calculates 
#                    the deuteron's momentum according to the theory of three body decay (DALITZ PLOT).
#                    We can immagine a general process: p-n -> 1 + 2 + 3, where 3 is deuteron, 1 and 2 pions.
#INPUT:   M = E_cms
#         m_1 = mass of thee firs pion
#         m_2 = mass of the second pion, m_3 is implicitly the deuteron's mass 
#OUTPUT:  deuteron_momentum = momentum of the formed deuteron
##############################################################################################################

def dalitz(M, m_1, m_2):
    
    m2_1 = pow(m_1, 2)
    m2_2 = pow(m_2, 2)
    
    #Drawn random invariant masses m_12 and m_23 uniformly within the kinematically allowed region
    m_12 = rng.uniform(pow(m_1 + m_2, 2), pow(M - d_m, 2)) #m_12
    
    E_2 = (m_12 - m2_1 + m2_2) / (2*np.sqrt(m_12))
    E_3 = (pow(M, 2) - m_12 - d_m2) / (2*np.sqrt(m_12))
    m_23_min = pow(E_2 + E_3, 2) - pow(np.sqrt(pow(E_2, 2) - m2_2) + np.sqrt(pow(E_3, 2) - d_m2), 2)
    m_23_max = pow(E_2 + E_3, 2) - pow(np.sqrt(pow(E_2, 2) - m2_2) - np.sqrt(pow(E_3, 2) - d_m2), 2)

    condition = False
    while(condition == False):
        
        m_23 = rng.uniform(pow(m_2 + d_m, 2), pow(M - m_1, 2)) #m_12

        if m_23_min < m_23 and m_23 < m_23_max:
            condition = True
            
        else:
            m_12 = rng.uniform(pow(m_1 + m_2, 2), pow(M - d_m, 2)) #m_23
    
            E_2 = (m_12 - m2_1 + m2_2) / (2*np.sqrt(m_12))
            E_3 = (pow(M, 2) - m_12 - d_m2) / (2*np.sqrt(m_12))
            m_23_min = pow(E_2 + E_3, 2) - pow(np.sqrt(pow(E_2, 2) - m2_2) + np.sqrt(pow(E_3, 2) - d_m2), 2)
            m_23_max = pow(E_2 + E_3, 2) - pow(np.sqrt(pow(E_2, 2) - m2_2) - np.sqrt(pow(E_3, 2) - d_m2), 2)

    #The deuteron's momentum is infered from the kinimatics of the three body decay
    p = (pow(M, 2) - pow(np.sqrt(m_12) + d_m, 2))*(pow(M, 2) - pow(np.sqrt(m_12) - d_m, 2))
    deuteron_momentum = (np.sqrt(p) / M) * 0.5

    return deuteron_momentum


# In[12]:


import numpy as np
from numpy import linalg as la


##############################################################################################################
#Function definition: Calculates the resulting deuteron four-momentum in the n + p -> d + gamma process 
#                     assuming p and n four-momenta are known, while gamma is completly unknown.
#                     NOTE: Deuteron and gamma are emitted back-to-back in the CoM system, in a random
#                     direction drawn from an isotropic distribution.
#INPUT:   proton =  four-momentum of the proton given in the CMS of the p-n
#         neutron = four-momentum of the neutron given in the CMS of the p-n
#OUTPUT:  p_deuteron = momentum's magnitude of the deuteron given in the CMS of the colliding e+e-
##############################################################################################################

def deuteron_momentum_rad(proton, neutron, beta):
    
    #Uniform emission throughout the entire solid angle 
    gamma_theta = np.arccos(rng.uniform(-1, 1))
    gamma_phi =   rng.uniform(0, 2 * np.pi)
    
    #Taking advantage of the relativistic kinematical invariants, calculates the deuteron's energy
    p = (la.norm(trivector(neutron)) + la.norm(trivector(proton))) * 0.5
    deuteron_energy = (p_m2 + n_m2 + d_m2) + 2 * (proton[0] * neutron[0] + pow(p, 2))  
    deuteron_energy = deuteron_energy / (2 * (proton[0] + neutron[0]))
    
    #From deuteron's energy infers the momentum magnitude
    deuteron_momentum = pow(pow(deuteron_energy, 2) - d_m2, 0.5)
    
    #Generate a random emission of gamma in the p-n CMS    
    gamma = np.array([np.sin(gamma_theta) * np.cos(gamma_phi), np.sin(gamma_theta) * np.sin(gamma_phi),
                      np.cos(gamma_theta)])
    
    #Obtain deuteron's tri-momentum which in the p-n CMS points in opposite versus respect gamma 
    deuteron_phase_space = gamma * (-deuteron_momentum)
    
    #Apply inverse LT to obtain the final deuteron's four-momentum in the CMS of e+e-
    deuteron = np.insert(deuteron_phase_space, 0, deuteron_energy)
    deuteron = boost(deuteron, -beta)
    
    p_deuteron = la.norm(trivector(deuteron))

    return p_deuteron


##############################################################################################################
#Function definition: Calculates the resulting deuteron four-momentum in the N1 + N2 -> d + pi process 
#                     assuming N1 and N2 four-momenta are known, while pion is completly unknown.
#                     NOTE: Deuteron and pion are emitted back-to-back in the CoM system, in a random
#                     direction drawn from an isotropic distribution.
#INPUT:   N1 = four-momentum of the nucleon given in the CMS of the pair
#         N2 = four-momentum of the nucleon given in the CMS of the pair
#         m =  mass of the pion, necessary to properly take into account the phace-space.
#OUTPUT:  p_deuteron = momentum's magnitude of the deuteron given in the CMS of the colliding e+e-
##############################################################################################################

def deuteron_momentum_pi(proton, neutron, beta, m):

    #Uniform emission throughout the entire solid angle 
    pi_theta = np.arccos(rng.uniform(-1, 1))
    pi_phi =   rng.uniform(0, 2 * np.pi)
    
    #Taking advantage of the relativistic kinematical invariants, calculates the deuteron's energy
    p = (la.norm(trivector(neutron)) + la.norm(trivector(proton))) * 0.5
    deuteron_energy = (p_m2 + n_m2 + d_m2 + pow(m, 2)) + 2 * (proton[0] * neutron[0] + pow(p, 2))  
    deuteron_energy = deuteron_energy / (2 * (proton[0] + neutron[0]))
    
    #From deuteron's energy infers the momentum magnitude
    deuteron_momentum = pow(pow(deuteron_energy, 2) - d_m2, 0.5)

    #Generate a random emission of the pion in the N1-N2 CMS    
    pi = np.array([np.sin(pi_theta) * np.cos(pi_phi), np.sin(pi_theta) * np.sin(pi_phi), np.cos(pi_theta)])
    
    #Obtain deuteron's tri-momentum which in the N1-N2 CMS points in opposite versus respect the pion 
    deuteron_phase_space = pi * (-deuteron_momentum)
    
    #Apply inverse LT to obtain the final deuteron's four-momentum in the CMS of e+e-
    deuteron = np.insert(deuteron_phase_space, 0, deuteron_energy)
    deuteron = boost(deuteron, -beta)
    
    p_deuteron = la.norm(trivector(deuteron))

    return p_deuteron


##############################################################################################################
#Function definition: Calculates the resulting deuteron four-momentum in the N1 + N2 -> d + pi + pi process 
#                     assuming N1 and N2 four-momenta are known, while pions are completly unknown.
#                     NOTE: Deuteron and pions are emitted back-to-back in the CoM system, in a random
#                     direction drawn from an isotropic distribution.
#INPUT:   proton =   four-momentum of the proton given in the CMS of the p-n
#         neutron =  four-momentum of the neutron given in the CMS of the p-n
#         m_1, m_2 = mass of the pions, necessary to properly take into account the phace-space
#OUTPUT:  p_deuteron = momentum's magnitude of the deuteron given in the CMS of the colliding e+e-
##############################################################################################################

def deuteron_momentum_pipi(proton, neutron, beta, m_1, m_2):
        
    M = proton[0] + neutron[0] #E_cms
    
    #Calculate the momentum of the deuteron in the CoM frame according to the paper's suggestion: DALITZ PLOT
    deuteron_momentum = dalitz(M, m_1, m_2)

    #From deuteron's momentum infers the energy value
    deuteron_energy = np.sqrt(pow(deuteron_momentum, 2) + d_m2)
    
    #Drawn the deuteron's direction of emission from an isotropic distribution in the N1-N2 CoM frame
    d_theta = np.arccos(rng.uniform(-1, 1))
    d_phi =   rng.uniform(0, 2 * np.pi)
    
    d = np.array([np.sin(d_theta) * np.cos(d_phi), np.sin(d_theta) * np.sin(d_phi), np.cos(d_theta)])
    
    #Obtain deuteron's tri-momentum in the N1-N2 CoM frame 
    deuteron_phase_space = d * (deuteron_momentum)
    
    #Apply inverse LT to obtain the final deuteron's four-momentum in the CMS of e+e-
    deuteron = np.insert(deuteron_phase_space, 0, deuteron_energy)
    deuteron = boost(deuteron, -beta)
    
    p_deuteron = la.norm(trivector(deuteron))
        
    return p_deuteron


# In[13]:


import numpy as np
from numpy import linalg as la


##############################################################################################################
#Function definition: Calculate k, the relative momentum of the nucleons in their CMS, and return theif four-
#                     momentum in the CMS of the pair, k and the Lorentz's beta used to boost.
#INPUT:   N_1 = four-momentum of the proton in the e+e- CMS
#         N_2 = four-momentum of the neutron in the e+e- CMS
#OUTPUT:  k =        magnitude of two times the relative momentum of the nucleons in their CMS
#         beta_CMS = Lorentz's beta used to boost in the CMS of the pair
#         N_1, N_2 = four-momentum of the nucleons in their CMS system of reference
##############################################################################################################

def CoM_system(N_1, N_2):

    #Calculate the four-momentum and tri-momentum of the p-n pair
    pair = N_1 + N_2
    pair_energy = pair[0]
    pair_phase_space = trivector(pair)

    beta_CMS = pair_phase_space / pair_energy #beta to boost in the CMS of the p-n pair

    #Apply the LT in order to obtain four-vectors in the p-n CMS system of reference
    N_1 = boost(N_1, beta_CMS)
    N_2 = boost(N_2, beta_CMS)
    
    #Calculate the relative momentum in the p-n CMS 
    k = la.norm(trivector(N_1) - trivector(N_2))
    
    return k, beta_CMS, N_1, N_2


# In[14]:


#Central part of the notebook

import numpy as np
from numpy import linalg as la


sigma_0_invers = 0.89 * 1e-6  #ub^-1, invers value of the Cross Section fenomenological parameter

p_deuter = []  #Definition of a list in which store the momentum of all deuteron candidates
event = []

#LOOP over each EVENT in which appear a pair previously saved, avoiding double counting
for i in range(0, len(event_array)):
    
    n = len(df.loc[event_array[i]]) #Number of anti-nucleon pairs in that event
    
    momentum_pair = [] #Temporarily store the deuteron's momentum from the pair combinations
    prob_pair = []     #Temporarily store the deuteron's formation probability from the pair combinations
    ID = []            #Temporarily store the identity of the nucleon that produced the deuteron
    
    
    #LOOP over each PAIR in the event
    for j in range(0,n):
        
        #Extract the four-momenta of the nucleon from the pair
        pair = df.loc[event_array[i],j]
        N_1 = np.array([pair.loc['E_N1'], pair.loc['px_N1'], pair.loc['py_N1'], pair.loc['pz_N1']]) #nucleon 1
        N_2 = np.array([pair.loc['E_N2'], pair.loc['px_N2'], pair.loc['py_N2'], pair.loc['pz_N2']]) #nucleon 2
        
        momentum_process = [] #Temporarily store the deuteron's momentum in that particular process
        prob_process = []     #Temporarily store the deuteron's formation probability in that process
        
        k, beta, N_1, N_2 = CoM_system(N_1, N_2) #boost pairs in the CoM and calculate the relative momentum 
        
        r = rng.uniform(0,1,4) #Array containing four random number uniformly distributed

        #If it is a p-n pair
        ######################################################################################################
        if pair.loc['PDG_N1'] == 2212 and pair.loc['PDG_N2'] == 2112:

            prob_0 = sigma_0_invers * sez_urto_pn_rad(k)
            if r[0] < prob_0:
                momentum_process.append(deuteron_momentum_rad(N_1, N_2, beta))
                prob_process.append(prob_0)
            prob_1 = sigma_0_invers * sez_urto_pn_pi0(k)
            if r[1] < prob_1:
                momentum_process.append(deuteron_momentum_pi(N_1, N_2, beta, pi0_m))
                prob_process.append(prob_1)
            prob_2 = sigma_0_invers * sez_urto_pn_pipi(k)
            if r[2] < prob_2:
                momentum_process.append(deuteron_momentum_pipi(N_1, N_2, beta, pi_m, pi_m))
                prob_process.append(prob_2)
            prob_3 = sigma_0_invers * sez_urto_pn_pi0pi0(k)
            if r[3] < prob_3:
                momentum_process.append(deuteron_momentum_pipi(N_1, N_2, beta, pi0_m, pi0_m))
                prob_process.append(prob_3)
            
            #Check if one deuteron is formed and temporarily store it in momentum_pair[]
            if len(momentum_process) == 1:
                momentum_pair.append(momentum_process[0])
                prob_pair.append(prob_process[0])
                ID.append(N_1[0])
                ID.append(N_2[0])
                
            #If two deuterons are formed choose randomly one of them according to the relative cross sections
            if len(momentum_process) > 1:
                prob_process = np.sort(prob_process)
                v = rng.uniform(0,1)
                if v < (prob_process[0] / (prob_process[0] + prob_process[1])):
                    momentum_pair.append(momentum_process[0])
                    prob_pair.append(prob_process[0])
                    ID.append(N_1[0])
                    ID.append(N_2[0])
                else:
                    momentum_pair.append(momentum_process[1])
                    prob_pair.append(prob_process[1])
                    ID.append(N_1[0])
                    ID.append(N_2[0])

        #If it is a p-p pair
        ######################################################################################################
        elif pair.loc['PDG_N1'] == 2212 and pair.loc['PDG_N2'] == 2212:

            prob_0 = sigma_0_invers * sez_urto_pp_pi(k)
            if r[0] < prob_0:
                momentum_process.append(deuteron_momentum_pi(N_1, N_2, beta, pi_m))
                prob_process.append(prob_0)
            prob_1 = sigma_0_invers * sez_urto_pp_pipi0(k)
            if r[1] < prob_1:
                momentum_process.append(deuteron_momentum_pipi(N_1, N_2, beta, pi_m, pi0_m))
                prob_process.append(prob_1)
            
            #Check if one deuteron is formed and temporarily store it in momentum[]
            if len(momentum_process) == 1:
                    momentum_pair.append(momentum_process[0])
                    prob_pair.append(prob_process[0])
                    ID.append(N_1[0])
                    ID.append(N_2[0])

            #If two deuterons are formed choose randomly one of them according to the relative cross sections
            if len(momentum_process) == 2:
                prob_process = np.sort(prob_process)
                v = rng.uniform(0,1)
                if v < (prob_process[0] / (prob_process[0] + prob_process[1])):
                    momentum_pair.append(momentum_process[0])
                    prob_pair.append(prob_process[0])
                    ID.append(N_1[0])
                    ID.append(N_2[0])   
                else:
                    momentum_pair.append(momentum_process[1])
                    prob_pair.append(prob_process[1])
                    ID.append(N_1[0])
                    ID.append(N_2[0])    

        #If it is a n-n pair
        ######################################################################################################
        elif pair.loc['PDG_N1'] == 2112 and pair.loc['PDG_N2'] == 2112:

            prob_0 = sigma_0_invers * sez_urto_nn_pi(k)
            if r[0] < prob_0:
                momentum_process.append(deuteron_momentum_pi(N_1, N_2, beta, pi_m))
                prob_process.append(prob_0)
            prob_1 = sigma_0_invers * sez_urto_nn_pipi0(k)
            if r[1] < prob_1:
                momentum_process.append(deuteron_momentum_pipi(N_1, N_2, beta, pi_m, pi0_m))
                prob_process.append(prob_1)

            #Check if one deuteron is formed and temporarily store it in momentum[]
            if len(momentum_process) == 1:
                momentum_pair.append(momentum_process[0])
                prob_pair.append(prob_process[0])
                ID.append(N_1[0])
                ID.append(N_2[0])
                
            #If two deuterons are formed choose randomly one of them according to the relative cross sections
            if len(momentum_process) > 1:
                prob_process = np.sort(prob_process)
                v = rng.uniform(0,1)
                if v < (prob_process[0] / (prob_process[0] + prob_process[1])):
                    momentum_pair.append(momentum_process[0])
                    prob_pair.append(prob_process[0])
                    ID.append(N_1[0])
                    ID.append(N_2[0])
                else:
                    momentum_pair.append(momentum_process[1])
                    prob_pair.append(prob_process[1])
                    ID.append(N_1[0])
                    ID.append(N_2[0])
    
    
    #If a deuteron is formed, exclude the involved nucleons from been used in the formation of other deuterons
    ##########################################################################################################
    if len(momentum_pair) == 1:
        p_deuter.append(momentum_pair[0])
        event.append(event_array[i])
            
    elif len(momentum_pair) > 1:
        event.append(event_array[i])
        
        if ID[0] == ID[2] or ID[0] == ID[3] or ID[1] == ID[2] or ID[1] == ID[3]:
            prob_pair = np.sort(prob_pair)
            v = rng.uniform(0,1,1)
            if v < (prob_pair[0] / (prob_pair[0] + prob_pair[1])):
                p_deuter.append(momentum_pair[0])
            else:
                p_deuter.append(momentum_pair[1])
        else:
            p_deuter.append(momentum_pair[0])
            p_deuter.append(momentum_pair[1])


#Translate lists into numpy array
p_deuter = np.array(p_deuter)
event = np.array(event)

#Print results
print('# deuterons produced:', len(p_deuter)) #Print the number of anti-deuterons formed
print(p_deuter)
print(event)
print(len(event))


# In[16]:


#Save results in files .dat
np.savetxt('prod_ass_anti_sez_urto_continuum_babar.dat', event)
np.savetxt('anti_event_urto_continuum_babar.dat', event_array)


# In[38]:


#Plot the momentum distribution of the produced deuterons

import matplotlib.pyplot as plt


N_ev = 5e07 #Total number of events generated

#Calculate the appropriate weights
weight = np.ones_like(p_deuter) / N_ev * 5.0 #Weight

#Plot histogram and corrispondent decorations 
plt.hist(p_deuter, bins=10, range=(0,4), weights=weight)
plt.xlabel('p [GeV]')
plt.ylabel(r'$\frac{1}{N_{ev}}$ $\frac{dn}{dp}$ $[Gev^{-1}]$')
plt.title(r'p distribution of the selected $d$ from $\Upsilon(1S)$')
plt.xlim((-0.1, 4.1))

#Show results
#plt.savefig('Sezione_urto/Istogrammi/p_deuter_sez_U1S.jpg', format='jpg', dpi=900)

plt.show()


# In[26]:


#Visual check of the kinematical allowed region of the three body decay, populated uniformally.


M = 2.5 #GeV, E_cms

list1 = []
list2 = []
m_1 = m_2 = pi_m
m2_1 = m2_2 = pow(pi_m, 2)

for i in range (0,10000):
    
    #Drawn random invariant masses m_12 and m_23 uniformly within the kinematically allowed region
    m_12 = rng.uniform(pow(m_1 + m_2, 2), pow(M - d_m, 2)) #m_12
    
    E_2 = (m_12 - m2_1 + m2_2) / (2*np.sqrt(m_12))
    E_3 = (pow(M, 2) - m_12 - d_m2) / (2*np.sqrt(m_12))
    m_23_min = pow(E_2 + E_3, 2) - pow(np.sqrt(pow(E_2, 2) - m2_2) + np.sqrt(pow(E_3, 2) - d_m2), 2)
    m_23_max = pow(E_2 + E_3, 2) - pow(np.sqrt(pow(E_2, 2) - m2_2) - np.sqrt(pow(E_3, 2) - d_m2), 2)

    condition = False
    while(condition == False):
        
        m_23 = rng.uniform(pow(m_2 + d_m, 2), pow(M - m_1, 2)) #m_12

        if m_23_min < m_23 and m_23 < m_23_max:
            condition = True
            
        else:
            m_12 = rng.uniform(pow(m_1 + m_2, 2), pow(M - d_m, 2)) #m_23
    
            E_2 = (m_12 - m2_1 + m2_2) / (2*np.sqrt(m_12))
            E_3 = (pow(M, 2) - m_12 - d_m2) / (2*np.sqrt(m_12))
            m_23_min = pow(E_2 + E_3, 2) - pow(np.sqrt(pow(E_2, 2) - m2_2) + np.sqrt(pow(E_3, 2) - d_m2), 2)
            m_23_max = pow(E_2 + E_3, 2) - pow(np.sqrt(pow(E_2, 2) - m2_2) - np.sqrt(pow(E_3, 2) - d_m2), 2)
    
    list1.append(m_23)
    list2.append(m_12)

    
#Plot
plt.scatter(list1, list2, s=0.5)
plt.xlabel(r'$m_{d \pi}^2$')
plt.ylabel(r'$m_{\pi\pi}^2$')
plt.title(r'Dalitz Plot: $p n \rightarrow d \pi \pi$')

plt.show()


# In[27]:


#Calculate the Branching Ratio
#CLEO result for anti-d production is BR = (3.36 +- 0.23 +- 0.25) e-05

#BR(N1 + N2 -> d + X)
d = len(p_deuter)
BR_d = d / N_ev
err_d = BR_d * np.sqrt(pow(np.sqrt(d)/d, 2) + pow(np.sqrt(N_ev)/N_ev, 2))


#Print Results
print('BR(N1 + N2 -> d + X) =', BR_d, '+-', err_d)


# In[ ]:




