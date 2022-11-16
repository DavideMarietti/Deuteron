#!/usr/bin/env python
# coding: utf-8

# In[1]:


##############################################################################################################
#I have generated 5.000.000 Upsilon(1S) events at 9.46 GeV and selected the p-n pairs that can be further
#analyzed to identify deuteron candidates with the Advanced Coalescence Model.
#NOTE: If not explicity specify, p-n pairs from lambda decay are not considered in this analysis.
#NOTE: An alternative to exclude weak decay products is to apply the cut: r = |r_p - r_n| < 3 fm. (Dal-Raklev)
#NOTE: Taking into account the smallness of the cross section probabilities, the contributions of triple and
#      higher deuteron formation per event have been neglected here.
##############################################################################################################

#Translates the Tree into DataFrame
import ROOT
import root_pandas
import numpy as np


#Show p-n pairs selected
df_pn = root_pandas.read_root('pn_U3S_recap.root', key='pn_U3S_recap_tree')
df_pn.tail(5)


# In[2]:


#Create a multi-index DataFrame array containing all possible pair combination in each relevant event: first 
#index is the "__event__" value, while second index is the "__cantidate__" value.

import pandas as pd


#First level of multi-index: array that contains the "__event__" in which have been found at least one pair
pairs = np.array(df_pn)
event_array = np.unique(np.array(df_pn.loc[:,'__event__']))

#Second level of the multi-index array: the "__cantidate__" value of the pairs related to that event
k = [1] #Initialize the array with one single element
for n in range(0, len(event_array)):
    i = len(pairs[np.where(pairs[:,2] == event_array[n])])
    h = np.arange(0,i)
    k = np.append(k,h)

k = np.delete(k, 0) #Delete the first element of the array, defined for the initialization

arrays = [pairs[:,2].astype(int), np.array(k)] #Define the first and second levels of the multi-index array
pairs = np.delete(pairs, [0,1,2,3,4,5,9,10,15,16], axis=1) #Delete negligible columns to improve readability

#Define the multi-index DataFrame array
df = pd.DataFrame(pairs, index=arrays)

df.rename(columns={df.columns[0]: "px_p+"}, inplace = True)
df.rename(columns={df.columns[1]: "py_p+"}, inplace = True)
df.rename(columns={df.columns[2]: "pz_p+"}, inplace = True)
df.rename(columns={df.columns[3]: "E_p+" }, inplace = True)
df.rename(columns={df.columns[4]: "px_n0"}, inplace = True)
df.rename(columns={df.columns[5]: "py_n0"}, inplace = True)
df.rename(columns={df.columns[6]: "pz_n0"}, inplace = True)
df.rename(columns={df.columns[7]: "E_n0" }, inplace = True)


df.head(10) #Visual check


# In[3]:


#GLOBAL VARIABLES are defined here.

#Mass definitions
d_m = 1.87561294257  #+- 5.7e-10 GeV
n_m = 0.93956542052  #+- 5.4e-10 GeV
p_m = 0.93827208816  #+- 2.9e-10 Gev

#Quadratic mass definitions
d_m2 = pow(d_m, 2)
n_m2 = pow(n_m, 2)
p_m2 = pow(p_m, 2)


# In[4]:


#BASIC FUNCTIONS are definde here.

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


# In[5]:


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


# In[6]:


##############################################################################################################
#Function definition: Calculate the probability to form a deuteron according to the different wave functions.
#INPUT:  q =     equivalent to k, it is half the relative momentum among the nucleons in their CoM system
#        sigma = fenomenologocal parameter rapresenting the size of the formation region
#        theta = angle between the direction of motion of the pair in their CoM and the z-axis
#        beta =  Lorentz's beta to boost in the CoM frame of the pair
#OUTPUT: w = probability to form deuteron
##############################################################################################################

def gauss_one(q, sigma):
    
    d =  3.2 #fm
    
    d2 = pow(d, 2)
    q2 = pow(q, 2)
    sigma2 = pow(sigma, 2)
    
    C = pow(d2 / (d2 + 4*sigma2), 3.0/2.0)
    
    #w = pow(np.pi, -3.0/2.0) * (1/pow(d, 3)) * np.exp(-q2/d2) #Module squared wave function
    w = 3 * C * np.exp(-q2*d2) #Probability to form deuteron
    
    return w


def gauss_two_phi0(q, sigma):
    
    #Parameter obtained from fits
    delta = 0.581
    d_1 =   3.979 #fm
    d_2 =   0.890 #fm
    
    d2_1 =   pow(d_1, 2)
    d2_2 =   pow(d_2, 2)
    q2 =     pow(q, 2)
    sigma2 = pow(sigma, 2)
    
    #Calculate the probability to form deuteron
    C_1 = pow(d2_1 / (d2_1 + 4*sigma2), 3.0/2.0)
    C_2 = pow(d2_2 / (d2_2 + 4*sigma2), 3.0/2.0)

    #w = pow(np.pi, -3.0/2.0)*((delta/pow(d_1, 3))*np.exp(-q2/d2_1) + ((1-delta)/pow(d_2, 3))*np.exp(-q2/d2_2))
    w = 3 * (C_1*delta*np.exp(-q2*d2_1) + C_2*(1-delta)*np.exp(-q2*d2_2))  #Probability to form deuteron

    return w


def gauss_two_r3(q, sigma):

    #Parameter obtained from fits
    delta = 0.247
    d_1 =   5.343 #fm
    d_2 =   1.810 #fm
    
    d2_1 =   pow(d_1, 2)
    d2_2 =   pow(d_2, 2)
    q2 =     pow(q, 2)
    sigma2 = pow(sigma, 2)
    
    #Calculate the probability to form deuteron
    C_1 = pow(d2_1 / (d2_1 + 4*sigma2), 3.0/2.0)
    C_2 = pow(d2_2 / (d2_2 + 4*sigma2), 3.0/2.0)

    #w = pow(np.pi, -3.0/2.0)*((delta/pow(d_1, 3))*np.exp(-q2/d2_1) + ((1-delta)/pow(d_2, 3))*np.exp(-q2/d2_2))
    w = 3 * (C_1*delta*np.exp(-q2*d2_1) + C_2*(1-delta)*np.exp(-q2*d2_2))  #Probability to form deuteron
    
    return w


# In[7]:


#VISUAL CHECK: Plot the cross section trends respect to k for the considered processes.

import matplotlib.pyplot as plt

y1 = []
y2 = []
y3 = []

x = np.linspace(0, 0.5, 1000)
#x = np.linspace(0, 10, 1000)

for i in range(0,len(x)):
    y1.append(gauss_one(x[i], 7.0))
for i in range(0,len(x)):
    y2.append(gauss_two_phi0(x[i], 7.0))
for i in range(0,len(x)):
    y3.append(gauss_two_r3(x[i], 7.0))

plt.plot(pow(x, 1), y1, 'k-.')
plt.plot(pow(x, 1) ,y2, 'g--')
plt.plot(pow(x, 1) ,y3, 'b:')

plt.xlabel(r'$q$ $[GeV]$')
plt.title(r'Probability Distribution')
plt.legend([r'one-gauss', r'two-gauss $\phi_0$-fit', r'two-gauss $r_3$-fit'], loc='upper right', fontsize=10)
#plt.yscale('log')

plt.show()


# In[8]:


import numpy as np
from numpy import linalg as la
from numpy.random import default_rng

rng = default_rng()

##############################################################################################################
#Function definition: calculates the resulting deuteron momentum in the n + p -> d + gamma process assuming
#                     p and n four-momenta are known, while gamma is completly unknown and its direction of
#                     emission is randomly generated.
#INPUT:   proton =  four-momentum of the proton given in the CMS of the p-n
#         neutron = four-momentum of the neutron given in the CMS of the p-n
#OUTPUT:  p_deuteron = momentum of the deuteron given in the CMS of the colliding e+e-.
##############################################################################################################

def deuteron_momentum(proton, neutron, beta, k):
    
    #Uniform emission throughout the entire solid angle 
    gamma_theta = np.arccos(rng.uniform(-1, 1))
    gamma_phi =   rng.uniform(0, 2 * np.pi)
    
    #Taking advantage of the relativistic kinematical invariants, calculates the deuteron's energy
    deuteron_energy = (p_m2 + n_m2 + d_m2) + 2*proton[0]*neutron[0] + 2*pow(k, 2)
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


# In[9]:


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

def CoM_system(proton, neutron):

    #Calculate the four-momentum and tri-momentum of the p-n pair
    pair = proton + neutron
    pair_energy = pair[0]
    pair_phase_space = trivector(pair)
    
    beta_CMS = pair_phase_space / pair_energy #beta to boost in the CMS of the p-n pair

    #Apply the LT in order to obtain four-vectors in the p-n CMS system of reference
    proton =  boost(proton,  beta_CMS)
    neutron = boost(neutron, beta_CMS)
    
    #Calculate the relative momentum in the p-n CMS 
    k = la.norm(trivector(proton) - trivector(neutron)) * 0.5
    
    return k, beta_CMS, proton, neutron


# In[ ]:


#Central part of the notebook

import numpy as np
from numpy import linalg as la


sigma = 8.81 #GeV^-1   Sigma: fenomenologocal parameter, rapresenting a sort of interaction volume

p_deuter = [] #Definition of a list in which store the momentum of all selected deuteron

#LOOP over each EVENT in which appear a pair previously saved, avoiding double counting
##############################################################################################################
for i in range(0, len(event_array)):
    
    n = len(df.loc[event_array[i]]) #Number of nucleon pairs in this even

    momentum = [] #Temporarily store the deuteron's momentum from the pair combinations
    prob = []     #Temporarily store the deuteron's formation probability from the pair combinations
    ID = []       #Temporarily store the identity of the nucleons that produced the deuteron
    
    #LOOP over each PAIR in the event
    ##########################################################################################################
    for j in range(0,n):
        
        #Extract the four-momenta of the nucleons from the selected pair
        pair = df.loc[event_array[i],j]
        proton =  np.array([pair.loc['E_p+'], pair.loc['px_p+'], pair.loc['py_p+'], pair.loc['pz_p+']])
        neutron = np.array([pair.loc['E_n0'], pair.loc['px_n0'], pair.loc['py_n0'], pair.loc['pz_n0']])
        
        q, beta, proton, neutron = CoM_system(proton, neutron) #boost pairs in the CoM and calculate the relative momentum 

        #SWITCH: select the correct condition depending if the ONE or TWO gaussian model is choosen 
        #if q < 0.25: #For the gauss_one case:  we take into consideration only pairs satisfing this condition
        if q < 0.5:   #For the gauss_two cases: we take into consideration only pairs satisfing this condition

            random_value = rng.uniform(0,1) #Random value uniformaly distributed between 0 and 1
            w = gauss_two_phi0(q, sigma)    #Probabilty that this pair can form a deuteron
            
            #Select deuteron candidates according to the advanced coalescence model
            if random_value < w:
                momentum.append(deuteron_momentum(proton, neutron, beta, q))
                prob.append(w)
                ID.append(proton[0])
                ID.append(neutron[0])


    #If a deuteron is formed, exclude the involved nucleons from been used in the formation of other deuterons
    if len(momentum) == 1:
        p_deuter.append(momentum[0])
         
    elif len(momentum) > 1:
        if ID[0] == ID[2] or ID[0] == ID[3] or ID[1] == ID[2] or ID[1] == ID[3]:
            prob = np.sort(prob)
            v = rng.uniform(0,1)
            if v < (prob[0] / (prob[0] + prob[1])):
                p_deuter.append(momentum[0])
            else:
                p_deuter.append(momentum[1])
        else:
            p_deuter.append(momentum[0])
            p_deuter.append(momentum[1])


#Translete lists into numpy array
p_deuter = np.array(p_deuter)

#Print results
print('# deuterons produced:', len(p_deuter))
print(p_deuter)

#Save results in files .dat
#np.savetxt('p_advance_U3S_babar.dat', p_deuter)


# In[12]:


#np.savetxt('p_advance_U3S_babar.dat', p_deuter)


# In[57]:


#Plot the momentum distribution of the produced deuterons

import matplotlib.pyplot as plt


#Calculate the appropriate weights
N_ev = 5e07 #Total number of events generated
weight = np.ones_like(p_deuter) / N_ev * 5.0 #Weight

#Plot histogram and corrispondent decorations 
plt.hist(p_deuter, bins=20, range=(0, 4), weights=weight)
plt.xlabel('p [GeV]')
plt.ylabel(r'$\frac{1}{N_{ev}}$ $\frac{dn}{dp}$ $[Gev^{-1}]$')
plt.title(r'p distribution of the selected $d$ from $\Upsilon(1S)$')
plt.xlim((-0.1, 4.1))

#Show results
#plt.savefig('Sezione_urto/Istogrammi/p_deuter_sez_U1S.jpg', format='jpg', dpi=900)

plt.show()


# In[ ]:


# Calculate the Branching Ratio
#CLEO result for anti-d production is BR = (3.36 +- 0.23 +- 0.25) e-05

#BR(p + n -> d + X)
d = len(p_deuter)
BR_d = d / N_ev
err_d = BR_d * np.sqrt(pow(np.sqrt(d)/d, 2) + pow(np.sqrt(N_ev)/N_ev, 2))


#Print Results
print('BR(p + n -> d + X) =', BR_d, '+-', err_d)


# In[ ]:




