#!/usr/bin/env python
# coding: utf-8

# In[6]:


#Perform the parbolic fit in order to find the chi-squared's minimum, i.e. the best fit fenomelogical parameter.

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq


#Data to be fitted
X = np.array([1.5, 1.4, 1.3, 1.2, 1.1, 1])
Y = np.array([3.07806, 1.795392, 0.2383639, 0.22571651, 2.00085, 3.468625])


#Standard form of quadratic function
def func(params, x):
    a, b, c = params
    return a * x * x + b * x + c


# Error function, that is, the difference between the value obtained by fitting curve and the actual value
def error(params, x, y):
    return func(params, x) - y


# Solving parameters
def sloveParam():
    
    p0 = [0.05, -6, 200]
    param = leastsq(error, p0, args=(X, Y))
    
    return param


#Display the final result
def solution():
    
    param = sloveParam()
    a, b, c = param[0]
    
    #Definition of the x's value vertice; i.e. the best fenomenological fit parameter parameter
    v = - b / (2*a)
    
    #Print useful infos
    print("a =", a, "b =", b, "c =", c)
    #print("cost:" + str(param[1]))
    print("The curve of solution is:")
    print("y = "+str(round(a,2))+"x*x+"+str(round(b,2))+"x+"+str(c))
    
    #Drawing the data sample
    plt.figure(figsize=(8,6))
    plt.scatter(X, Y, color="green", label="data sample", linewidth=2)

    # Drawing fitted lines
    x = np.linspace (min(X), max(X), 100) #Draw 100 continuous points directly from the x's minimum to x's maximum
    y = a * x * x + b * x + c #Function
    plt.plot(x,y,color="red",label="solution line",linewidth=2)
    plt.legend()# Draw Legend
    plt.show()
    
    #Display the best fenomenological fit parameter parameter and the associated error
    error1 = ((-b + np.sqrt(pow(b, 2) - 4*a*(c-1))) / a) * 0.5
    error2 = ((-b - np.sqrt(pow(b, 2) - 4*a*(c-1))) / a) * 0.5
    error = abs(error1 - error2) * 0.5
    print(error1)
    print(error2)
    print("Best fit fenomenological parameter:", v, "+-", error)


#Call the function
solution()


# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec  # for unequal plot boxes
import scipy.optimize
    

d_m = 1.87561294257  #+- 5.7e-10 GeV
d_m2 = pow(d_m, 2)

#Define the fitting function
def maxwell_distr(x, a, b):
    return a * pow(x,2) * np.exp(-b * np.sqrt(pow(x,2) + d_m2))

# read in spectrum from data file
# x=p_cm, y=dBR/dp, dy=y uncertainty
x, y, ys = np.loadtxt("Spectrum.txt", unpack=True)

# initial guesses for fitting parameters
a = 1.0e-02
b = 4.0

# fit data using SciPy's Levenberg-Marquart method
nlfit, nlpcov = scipy.optimize.curve_fit(maxwell_distr, x, y, p0=[a, b], sigma=dy)

# unpack fitting parameters
a, b = nlfit

# unpack uncertainties in fitting parameters from diagonal of covariance matrix
da, db = \ [np.sqrt(nlpcov[j,j]) for j in range(nlfit.size)]

# create fitting function from fitted parameters
f_fit = np.linspace(0.0, 3.0, 100)
s_fit = maxwell_distr(f_fit, a, b)

# Calculate residuals and reduced chi squared
resids = y - maxwell_distr(x, a, b)
redchisqr = ((resids/dy)**2).sum()/float(x.size-6)

# Create figure window to plot data
fig = plt.figure(1, figsize=(8,8))
gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

# Top plot: data and fit
ax1 = fig.add_subplot(gs[0])
ax1.plot(f_fit, s_fit)
ax1.errorbar(f, s, yerr=ds, fmt='or', ecolor='black')
ax1.set_xlabel('p')
ax1.set_ylabel('dBR/dp')
ax1.text(0.7, 0.95, 'a = {0:0.1x}$\pm${1:0.1x}'
         .format(a, da), transform = ax1.transAxes)
ax1.text(0.7, 0.90, 'b = {0:0.2x}$\pm${1:0.2x}'
         .format(b, db), transform = ax1.transAxes)


# Bottom plot: residuals
ax2 = fig.add_subplot(gs[1])
ax2.errorbar(x, resids, yerr = dy, ecolor="black", fmt="ro")
ax2.axhline(color="gray", zorder=-1)
ax2.set_xlabel('p')
ax2.set_ylabel('dBR/dp')
ax2.set_ylim(-20, 20)
ax2.set_yticks((-20, 0, 20))

plt.show()


# In[101]:


from numpy.polynomial import polynomial as P
import numpy as np
import matplotlib.pyplot as plt


#Data to be fitted
#-----------------------------------------------------------------------------------#
X = np.array([8,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8])
Y = np.array([4.12680425155675,2.21912079118501,
              2.32423857386899,3.26771080011906,5.02279753675983,5.99302203232199,7.90378630295337,
              11.4963818392275,13.0362019322854])

err = np.array([1.03140646011561,0.613958544351195,
                0.552733831781154,0.662216651857646,0.852279692036009,0.950672318088458,1.12463572210897,
                1.38985772689866,1.49083509072659])  


#--------------------------------------#
c = P.polyfit(X,Y,2,w=1/err)

a = c[2]
b = c[1]
c = c[0]
v = - b / (2*a)

#Print useful info
#---------------------------------------------------#
print("a =", a, "b =", b, "c =", c)
#print("cost:" + str(param[1]))
#print("The curve of solution is:")
#print("y = "+str(round(a,2))+"x*x+"+str(round(b,2))+"x+"+str(c))

#Drawing the data sample
plt.figure(figsize=(8,6))
plt.errorbar(X, Y, yerr= err, color="green", label="data sample", linestyle='', marker='o')

# Drawing fitted lines
x = np.linspace (min(X), max(X), 100) #Draw 100 continuous points directly from the x's minimum to x's maximum
y = a * x * x + b * x + c #Function
plt.plot(x,y,color="red",label="solution line",linewidth=2)
plt.legend()# Draw Legend
plt.show()

#Display the best fenomenological fit parameter parameter and the associated error
#error1 = ((-b + np.sqrt(pow(b, 2) - 4*a*(c-1))) / a) * 0.5
#error2 = ((-b - np.sqrt(pow(b, 2) - 4*a*(c-1))) / a) * 0.5
#error = abs(error1 - error2) * 0.5

print("Best fit fenomenological parameter:", v)


# In[102]:


err_1 = (-b + 2*np.sqrt(a))/(2*a)
print(err_1)

print('errore sup/inf:', err_1 - v)


# In[ ]:




