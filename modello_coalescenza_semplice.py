#!/usr/bin/env python
# coding: utf-8

# In[153]:


##############################################################################################################
#I have generated 5.000.000 Upsilon(1S/2S/3S) events and selected the (anti-)p-n pairs that can be  
#further analyzed to identify the (anti-)deuteron candidates with the Simple Coalescence Model.
#NOTE: if not explicity specify, p-n pairs from lambda decay are not considered in this analysis.
#NOTE: an alternative to exclude weak decay products is to apply the cut: r = |r_p - r_n| < 3 fm. (Dal-Raklev)
##############################################################################################################

#Translates the Tree into DataFrame: each row is a (anti-)deuterium candidate and in its columns are stored
#the four-momenta of the (anti-)p-n mother particles and the deuterium itself.

import ROOT
import root_pandas


df = root_pandas.read_root('deuterons_continuum_recap.root', key='deuterons_continuum_recap_tree')
df.head(5)


# In[154]:


#Basic function definitions are reported here, at the beginning of the notebook.

import numpy as np


##############################################################################################################  
#Function definition: extracs tri-vector from a four-vector.
#INPUT:   quadrivector = four-vector
#OUTPUT:  tri-vector = tri-vector in form of numpy array
##############################################################################################################

def trivector(quadrivector):
    
    quadrivector = np.array(quadrivector)
    tri_vector =   np.array([quadrivector[1], quadrivector[2], quadrivector[3]])

    
    return tri_vector


##############################################################################################################
#Function definition: calculates the QUADRATIC norm in minkowsky space, metric = diag(+,-,-,-).
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


# In[155]:


import numpy as np
from numpy import linalg as la


##############################################################################################################
#Function deinition: given a certain four-vector, it operates the LT along a certain direction defined by the
#                    Lorentz's beta and returns the trasformed four-vector.
#
#INPUT:   quadrivector = four-vector in phase space, i.e. quadrivector = [E, px, py, pz]
#         bata = Lorentz's beta parameter, i.e. b = [bx, by, bz]
#OUTPUT:  quadrivector_transformed = LT trasformed four-vector, i.e. quadrivector_tranf = [E', px', py', pz']
##############################################################################################################

def boost(quadrivector, beta):

    #Calculates factors for the LT
    beta = np.array(beta)
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


# In[2]:


#Mass definitions
deuterium_mass = 1.87561294257  #(+- 5.7e-10 GeV)
neutron_mass =   0.93956542052  #(+- 5.4e-10 GeV)
proton_mass =    0.93827208816  #(+- 2.9e-10 Gev)

#Quadratic mass definitions
deuterium_mass2 = pow(deuterium_mass, 2)
neutron_mass2 =   pow(neutron_mass, 2)
proton_mass2 =    pow(proton_mass, 2)


# In[157]:


import numpy as np
from numpy import linalg as la
from numpy.random import default_rng


##############################################################################################################
#Function definition: calculates the resulting deuterium four-momentum in the n + p -> d or n + p -> d + gamma
#                     process assuming p and n four-momenta are known, while gamma is completly unknown and
#                     its direction of emission is randomly generated. Several options are avaible below.
#INPUT:   proton =  four-momentum of the proton given in the CMS of the p-n
#         neutron = four-momentum of the neutron given in the CMS of the p-n
#OUTPUT:  deuterium = four-momentum of the deuter given in the CMS of the colliding e+e-.
##############################################################################################################


##############################################################################################################
#OPTION 1: the gamma emission is considered. Radioactive capture is the dominating deuterium formation process
#          at low CMS momentum differences required by the coalescence model, four-momentum is preserved.
#          Process: p + n -> d + gamma
##############################################################################################################

def deuterium_quadrivector_opt1(proton, neutron, beta):
    
    #Taking advantage of the relativistic kinematical invariants, calculates the deuterium's energy
    p = (la.norm(trivector(neutron)) + la.norm(trivector(proton))) * 0.5
    deuterium_energy = (proton_mass2 + neutron_mass2 + deuterium_mass2) + 2 * (proton[0] * neutron[0] + pow(p, 2))  
    deuterium_energy = deuterium_energy / (2 * (proton[0] + neutron[0]))
    
    #From deuterium's energy infers the momentum magnitude
    deuterium_momentum = pow(pow(deuterium_energy, 2) - deuterium_mass2, 0.5)

    #Uniform emission throughout the entire solid angle 
    gamma_theta = np.arccos(default_rng().uniform(-1, 1))
    gamma_phi =   default_rng().uniform(0, 2 * np.pi)
    
    #Generate a random emission of gamma in the p-n CMS    
    gamma = np.array([np.sin(gamma_theta) * np.cos(gamma_phi), np.sin(gamma_theta) * np.sin(gamma_phi),
                      np.cos(gamma_theta)])
    
    #Obtain deuterium's tri-momentum which in the p-n CMS points in opposite versus respect gamma 
    deuterium_phase_space = gamma * (-deuterium_momentum)
    
    #Apply inverse LT to obtain the final deuterium's four-momentum in the CMS of e+e-
    deuterium = np.insert(deuterium_phase_space, 0, deuterium_energy)
    deuterium = boost(deuterium, -beta)
    
    return deuterium
        

##############################################################################################################
#OPTION 2: I'm considering p-n pairs as a unique particle, i.e. the deuterium. Gamma emission is ignored and
#          four-momentum is not preserved.
#          Process: p + n -> d
##############################################################################################################

def deuterium_quadrivector_opt2(proton, neutron, beta):
    
    beta_norm2 = la.norm(beta)           #Lorentz'beta
    gamma = pow(1 - beta_norm2, -0.5)    #Lorentz'gamma
    
    deuterium_E = proton[0] + neutron[0]
    deuterium_E = deuterium_E * gamma
    deuterium_p = deuterium_E * beta
    
    deuterium = np.insert(deuterium_p, 0, deuterium_E)
    
    return deuterium


##############################################################################################################
#OPTION 3: as option 2, but deuterium's for-momentum is simply obtained from Modular Analyis variable values.
#INPUT:    n is the row index of the array containing all deuterium candidates, i.e a integer number
##############################################################################################################

def deuterium_quadrivector_opt3(n):
    
    deuterium_E = df.iloc[n,12]
    deuterium_p = np.array([df.iloc[n,7], df.iloc[n,8], df.iloc[n,9]])
    
    deuterium = np.insert(deuterium_p, 0, deuterium_E)
    
    return deuterium


# In[158]:


import numpy as np
from numpy import linalg as la


##############################################################################################################
#Function definition: selects (anti-)deuterium candidates applying the SIMPLE COALESCENCE condition and 
#                     returns their four-momentum.
#INPUT:   proton =  four-momentum of the proton in the e+e- CMS
#         neutron = four-momentum of the neutron in the e+e- CMS
#         k_cut =   coalescence momentum paramenter, i.e a number
#OUTPUT:  deuterium = four-momentum of the (anti-)deuterium in the e+e- CMS
##############################################################################################################


def deuterium (proton, neutron, k_cut):

    proton =  np.array(proton)
    neutron = np.array(neutron)
    
    #Calculate total the four-momentum and tri-momentum of the p-n pair
    pair = proton + neutron
    pair_phase_space = trivector(pair)
    pair_energy = pair[0]
    
    #Calculate the Lorentz's beta to boost in the CMS of the p-n pair
    beta_CMS = pair_phase_space / pair_energy
    beta_CMS = np.array(beta_CMS)

    #Apply the LT in order to obtain four-vectors in the p-n CMS
    proton =  boost(proton, beta_CMS)
    neutron = boost(neutron, beta_CMS)
    
    #Calculate the relative momentum in the p-n CMS
    k = (trivector(proton) - trivector(neutron)) * 0.5
    k = la.norm(k)
    k_min = abs((2*deuterium_mass*(proton[0] + neutron[0]) - deuterium_mass2 - proton_mass2 - neutron_mass2 - 2*proton[0]*neutron[0])) 
    k_min = np.sqrt(k_min * 0.5)
    
    #Apply the SIMPLE COALESCENCE condition in order to select deuteron candidates
    if k_min < k and k < k_cut:
        return deuterium_quadrivector_opt1(proton, neutron, beta_CMS) #Return the (anti-)d four-momentum
    else:
        return np.nan


# In[164]:


#Central part of the notebook

import numpy as np
from numpy import linalg as la


max_index = len(df.index)   #Total number of (anti-)p-n pairs saved
k_cut = 0.05941#GeV         #Coalescence momentum (fenomenologocal parameter)

p_anti_d = []          #Definition of a list in which store the momentum of all anti-d candidates
p_d = []               #Definition of a list in which store the momentum of all d candidates

counter=0
#Loop over each (anti-)p-n pair previously saved as a potential (anti-)deuterium candidate
for n in range(0,max_index):
    
    event = df.iloc[n,2] #Event considered

    #If PDG < 0
    if df.iloc[n,6] < 0: 
        anti_proton =  np.array([df.iloc[n,18], df.iloc[n,13], df.iloc[n,14], df.iloc[n,15]]) #four-momentum
        anti_neutron = np.array([df.iloc[n,24], df.iloc[n,19], df.iloc[n,20], df.iloc[n,21]]) #[E,p_x,p_y,p_z]
        
        #Return anti-deuterium four-momentum and event if it satisfies the SIMPLE COALESCENCE condition
        anti_deuter = deuterium(anti_proton, anti_neutron, k_cut)

        if anti_deuter is not np.nan:
            p_norm = la.norm(trivector(anti_deuter))  #Momentum of the selected anti-deututerium
            p_anti_d.append(p_norm)
            
            #print('dbar =', anti_proton)
            counter+=1
            if counter%100 == 0:
                print(counter)
        
    else:
        proton =  np.array([df.iloc[n,18], df.iloc[n,13], df.iloc[n,14], df.iloc[n,15]])
        neutron = np.array([df.iloc[n,24], df.iloc[n,19], df.iloc[n,20], df.iloc[n,21]])

        deuter = deuterium(proton, neutron, k_cut)
                                                        
        if deuter is not np.nan:
            p_norm = la.norm(trivector(deuter))  #Momentum of the selected deututerium
            p_d.append(p_norm)
            
            #print('d    =', deuter)
        

#Translete lists into numpy array
p_anti_d = np.array(p_anti_d)
p_d = np.array(p_d)


# In[165]:


print("# anti-deuterons produced:", len(p_anti_d))
print(p_anti_d)


# In[169]:


np.savetxt('p_simple_continuum_babar.dat', p_anti_d)


# In[170]:


##############################################################################################################
#Paper Gustafson-Hakkinen (deuteron)
##############################################################################################################


# In[171]:


#Momentum distribution of the selected deuterons (Fig.5)

import matplotlib.pyplot as plt


#Calculate the appropriate weights
N_ev = 5e07 #Total number of events generated
weight = np.ones_like(p_d) / N_ev * 5.0 #Weight

#Plot histogram and corrispoondent decorations 
plt.hist(p_d, bins=20, range=(0, 4), weights=weight)
plt.xlabel('p [GeV]')
plt.ylabel(r'$\frac{1}{N_{ev}}$ $\frac{dn}{dp}$ [c\Gev]')
plt.title(r'p distribution of the selected $d$ from $\Upsilon(1S)$')
plt.xlim((-0.1, 4.1))

#Show results
#plt.savefig('Coalescenza_semplice/Istogrammi/p_deuterium_distr.jpg', format='jpg', dpi=900)

plt.show()


# In[86]:


#Calculate the Branching Ratio
#CLEO result for anti-d production is BR = (3.36 +- 0.23 +- 0.25) e-05
#CLEO result for d anti-d associate production is BR = (0.9 +- 0.6) e-03

#BR(p + n -> dbar + gamma)
anti_d = len(p_anti_d)
BR_anti_d = anti_d / N_ev
err_anti_d = BR_anti_d * np.sqrt(pow(np.sqrt(anti_d)/anti_d, 2) + pow(np.sqrt(N_ev)/N_ev, 2))

#BR(p + n -> d + gamma)
d = len(p_d)
BR_d = d / N_ev
err_d = BR_d * np.sqrt(pow(np.sqrt(d)/d, 2) + pow(np.sqrt(N_ev)/N_ev, 2))

#Print Results
print('BR(p + n -> dbar + gamma) =', BR_anti_d, '+-', err_anti_d)
print('BR(p + n -> d + gamma) =   ', BR_d, '+-', err_d)


# In[39]:


#Momentum distribution of the p-n pairs. (Fig.2)

import matplotlib.legend as lg

#Load data from continuum distribution
k_d_continuum = np.loadtxt('k_pairs_distr_continuum.dat')

#Calculate the appropriate weights
weight = np.ones_like(k_d) / (N_ev * k_d * k_d) #Weight
weight_continuum = np.ones_like(k_d_continuum) / (9774466 * k_d_continuum * k_d_continuum)

#Plot histogram and corrispondent decorations 
plt.hist(k_d, bins=100, range=(0.05,0.5), weights=weight)
plt.hist(k_d_continuum, bins=100, range=(0.05,0.5), weights=weight_continuum, color='r')
plt.xlabel('k [GeV]')
plt.ylabel(r'$\frac{1}{N_{ev}}$ $\frac{dn}{d^3k}$ $[Gev^{-3}]$')
plt.title(r'Momentum distribution of the $p-n$ pairs')
plt.legend([r'$\Upsilon(1S)$ (9.46 GeV)', r'$q\bar{q}$ (10.58 GeV)'])
plt.xlim((0, 0.5))


#Show results
#plt.savefig('Coalescenza_semplice/Istogrammi/k_pairs_distribution_U1S.jpg', format='jpg', dpi=900)

plt.show()


# In[14]:


##############################################################################################################
#Paper Artoisenet-Braaten (anti-deuterium)
##############################################################################################################


# In[15]:


#Distribution of (anti-)p-n pairs as a function of the relative momentum k (considering also pairs from lambda)

import matplotlib.pyplot as plt

k_anti_d = np.loadtxt('k_anti_pairs_distr_U1S.dat')

#Calculate the appropriate weights
N_ev = 5e07
weight = np.ones_like(k_anti_d) / N_ev #Weight

#Plot histogram and corrispoondent decorations 
plt.hist(k_anti_d, bins=40, range=(0, 2), weights=weight)
plt.xlabel('k [GeV]')
plt.ylabel('Event fraction / 50 MeV bin')
plt.title(r'Momentum distribution of the $\bar{p}-\bar{n}$ pairs (included $\Lambda$ decays)')
plt.ylim((1e-6, 0.01))
plt.yscale('log')

#Show the k^2 trend at k < 200 Mev of the event fraction of (anti-)p-n pairs
x = np.linspace(0,0.5,100)
y = 0.2*0.2*x*x
plt.plot(x,y, 'r')

#Show results
#plt.savefig('Coalescenza_semplice/Istogrammi/k_pairs_distribution_U1S.jpg', format='jpg', dpi=900)

plt.show()


# In[7]:


import numpy as np

BR_teo = 2.81*pow(10,-5)
err_BR_teo = 0.57*pow(10,-5)

BR_exp = 6.0*pow(10,-8)
err_BR_exp = 2.8*pow(10,-8)

BR= BR_exp / BR_teo

err_BR = BR * np.sqrt(pow(err_BR_exp/BR_exp, 2) + pow(err_BR_teo/BR_teo, 2))

print(err_BR)


# In[ ]:




