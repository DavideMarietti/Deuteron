#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Perform the parabolic fit in order to find the chi-squared's minimum, i.e. the best fit fenomelogical parameter.

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq


#Data to be fitted
X = np.array([74, 72, 70, 68, 66, 64, 62, 60, 58, 56])
Y = np.array([37.6659884851318, 20.3862279838524, 9.24766707861548, 2.69072558645573, 0.534580457662901,
              1.32873938915425, 4.95480826433039, 11.1006907662353, 18.0569226633534, 28.4992253327862])
err_Y = np.array([3.66170988260421, 2.59739558007358, 1.70696888552154, 0.902465661369881, 0.383726060438828,
                  0.56226402585359, 1.07466349043977, 1.54763577948001, 1.89506139563689, 2.26556604504505])

#Standard form of quadratic function
def func(params, x):
    a, b, c = params
    return a * x * x + b * x + c


#Error function, that is, the difference between the value obtained by fitting curve and the actual value
def error(params, x, y):
    return func(params, x) - y


#Solving parameters
def sloveParam():
    
    p0 = [0.05, -6, 200]
    param = leastsq(error, p0, args=(X, Y))
    
    return param

#Parameters
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

#Drawing fitted lines
x = np.linspace (min(X), max(X), 100) #Draw 100 continuous points directly from the x's minimum to x's maximum
y = a * x * x + b * x + c #Function
plt.plot(x,y,color="red",label="fit line",linewidth=2)
plt.legend()# Draw Legend
plt.show()
    
#Display the best fenomenological fit parameter parameter and the associated error
error1 = ((-b + np.sqrt(pow(b, 2) - 4*a*(c-1))) / a) * 0.5
error2 = ((-b - np.sqrt(pow(b, 2) - 4*a*(c-1))) / a) * 0.5
error = abs(error1 - error2) * 0.5
print(error1)
print(error2)
print("Best fit fenomenological parameter:", v, "+-", error)


# In[2]:


###################################################################################################################
#Create a minimization plot that contains error bars

fig, ax = plt.subplots()

#Errorbar graphics for data and MC are defined
ax.errorbar(X, Y, yerr=err_Y, fmt='bs', elinewidth=0.7, barsabove=True)

#Decorations
ax.set_title(r'Simple Coalescence - $\Upsilon(1S)+\Upsilon(2S)+\Upsilon(3S)$')
ax.set_ylabel(r'$\chi^2$')
ax.set_xlabel(r'k [MeV/c]')
ax.yaxis.grid(True, linestyle='--')

#Plot layout
plt.tight_layout()
plt.xticks(np.arange(min(X)-2, max(X)+2, 2))


###################################################################################################################
#Show the result

plt.savefig('Istogrammi_chi2/simple_chi2_U1S_U2S_U3S_babar.jpg', format='jpg', dpi=900)
plt.show()


# In[3]:


#Experimental data from BaBar, CLEO and Argus are reported here.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


def binning(x):
    
    x_new = []
    n = len(x)
    for i in range(0,n-1):
        x_new.append(x[i] + (x[i+1]-x[i])*0.5)
        
    return np.array(x_new)


def width(a):
    
    a_new = []
    n = len(a)
    for i in range(0,n-1):
        a_new.append(a[i+1]-a[i])
        
    return np.array(a_new)


#Babar U1S
xbins_U1S_babar = np.array([0.35, 0.65, 0.85, 1.0, 1.20, 2.25])
xerr_U1S_babar =  0.5*width(xbins_U1S_babar)
xbins_U1S_babar = binning(xbins_U1S_babar)
y_U1S_babar =     np.array([12.032, 27.723, 34.373, 28.149, 6.103])
y_U1S_babar =     y_U1S_babar * 1.0e-06
yerr_U1S_babar =  np.array([7.120, 8.700, 13.047, 9.569, 2.876])
yerr_U1S_babar =  yerr_U1S_babar * 1.0e-06

#Babar U2S
xbins_U2S_babar = np.array([0.35, 0.55, 0.7, 0.8, 0.9, 1.0, 1.15, 1.30, 1.55, 2.25])
xerr_U2S_babar =  0.5*width(xbins_U2S_babar)
xbins_U2S_babar = binning(xbins_U2S_babar)
y_U2S_babar =     np.array([15.775, 17.515, 23.924, 25.115, 31.250, 19.896, 17.149, 12.479, 3.322])
y_U2S_babar =     y_U2S_babar * 1.0e-06
yerr_U2S_babar =  np.array([2.106, 2.289, 3.388, 3.663, 4.578, 3.480, 3.480, 2.839, 1.465])
yerr_U2S_babar =  yerr_U2S_babar * 1.0e-06

#Babar U3S
xbins_U3S_babar = np.array([0.35, 0.6, 0.7, 0.8, 0.9, 1.0, 1.10, 1.25, 1.45, 2.25])
xerr_U3S_babar =  0.5*width(xbins_U3S_babar)
xbins_U3S_babar = binning(xbins_U3S_babar)
y_U3S_babar =     np.array([12.658, 19.098, 20.157, 13.683, 23.972, 19.182, 13.464, 10.052, 5.434])
y_U3S_babar =     y_U3S_babar * 1.0e-06
yerr_U3S_babar =  np.array([7.713, 3.084, 3.393, 5.384, 5.075, 5.387, 5.848, 4.150, 1.441])
yerr_U3S_babar =  yerr_U3S_babar * 1.0e-06

#Babar Continuum
xbins_continuum_babar = np.array([0.35, 0.60, 0.75, 0.85, 0.95, 1.05, 1.25, 1.4, 1.65, 2.25])
xerr_continuum_babar =  0.5*width(xbins_continuum_babar)
xbins_continuum_babar = binning(xbins_continuum_babar)
y_continuum_babar =     np.array([2.080, 1.691, 2.261, 2.249, 2.507, 2.016, 2.077, 1.478, 0.593])
y_continuum_babar =     y_continuum_babar * 1.0e-06
yerr_continuum_babar =  np.array([0.563, 0.527, 0.713, 0.701, 0.773, 0.575, 0.611, 0.425, 0.366])
yerr_continuum_babar =  yerr_continuum_babar * 1.0e-06

#CLEO U1S
xbins_cleo = np.array([0.45, 0.65, 0.85, 1.05, 1.25, 1.45])
xerr_cleo =  0.5*width(xbins_cleo)
xbins_cleo = binning(xbins_cleo)
y_cleo =     np.array([2.2, 3.0, 2.9, 2.5, 1.5])
y_cleo =     y_cleo * 1.0e-05
yerr_cleo =  np.array([0.54, 0.45, 0.50, 0.45, 0.45])
yerr_cleo =  yerr_cleo * 1.0e-05

#Argus U1S+U2S
xbins_argus = np.array([0.45, 0.7, 0.95, 1.2, 1.45, 1.7])
xerr_argus =  0.5*width(xbins_argus)
xbins_argus = binning(xbins_argus)
y_argus =     np.array([3.7, 6.1, 5.7, 3.7, 1.0])
y_argus =     y_argus * 1.0e-05
u_err_argus = np.array([3.61, 4.23, 4.23, 4.01, 3.82])
u_err_argus = u_err_argus * 1.0e-05
l_err_argus = np.array([2.02, 2.75, 2.75, 2.42, 1.02 ])
l_err_argus = l_err_argus * 1.0e-05


# In[48]:


#Show experimental data's plot
###################################################################################################################

fig, ax = plt.subplots()


#Errorbar plots
ax.errorbar(xbins_U1S_babar, y_U1S_babar, yerr=yerr_U1S_babar, fmt='bs', elinewidth=0.7, barsabove=True)
ax.errorbar(xbins_cleo, y_cleo,  yerr=yerr_cleo, fmt='r^', elinewidth=0.7)
ax.errorbar(xbins_argus, y_argus,  yerr=[l_err_argus, u_err_argus], fmt='gd', elinewidth=0.7)

#Decorations
ax.set_title(r'Experimental Momentum Distributions - $\Upsilon(1S)$')
ax.set_ylabel(r'$\frac{dBR}{dp}$ $[Gev^{-1}c]$')
ax.set_xlabel(r'p [GeV/c]')
ax.legend([r'BaBar', r'CLEO', r'Argus ($\Upsilon(1S)+\Upsilon(2S))$'], fontsize=10)
ax.yaxis.grid(True, linestyle='--')

#Plot layout
#plt.ylim((0, 10e-05))
plt.xlim((0, 2.3))
plt.tight_layout()


###################################################################################################################
#Show the result

#plt.savefig('Istogrammi_spettro/exp_distr_U1S.jpg', format='jpg', dpi=900)
plt.show()


# In[38]:


#Spectra in momentum from model output once set the best fit parameter value.

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


N_ev = 5.0e+07 #Number of events generated


def make_error_boxes(ax, xdata, ydata, xerror, yerror, fc='xkcd:powder blue', ec='None', alpha=0.5):

    # Create list for all the error patches
    errorboxes = []

    # Loop over data points; create box from errors at each point
    for xc, yc, xe, ye in zip(xdata, ydata, xerror.T, yerror.T):
        rect = Rectangle((xc-xe, yc-ye), xe*2, ye*2)
        errorboxes.append(rect)

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=fc, alpha=alpha, edgecolor=ec)

    # Add collection to axes
    ax.add_collection(pc)

    # Plot errorbars
    ax.errorbar(xdata, ydata, xerr=xerror, fmt='None', ecolor='xkcd:bluegrey', elinewidth=0.7, label='Experiment (BaBar)')

    
#Model data
##################################################################################################################

#Load data
y_simple =  np.loadtxt("./Coalescenza_semplice/p_simple_U1S_babar.dat", unpack=True)
y_sez =     np.loadtxt("./Sezione_urto/p_sez_urto_U1S_babar.dat", unpack=True)
y_advance = np.loadtxt("./Coalescenza_avanzata/p_advance_U1S_babar.dat", unpack=True)
y_boost =   np.loadtxt("./Coalescenza_avanzata/p_advance_U1S_boost_babar.dat", unpack=True)

#xbins = np.array([0.35, 0.60, 0.75, 0.85, 0.95, 1.05, 1.25, 1.4, 1.65, 2.25])
x_simple =      np.linspace(0.16, 2.46, 12)
bins_x_simple = np.array(x_simple[:-1]) + 0.5 * (x_simple[1] - x_simple[0])
x_sez =         np.linspace(0.2, 2.5, 12)
bins_x_sez =    np.array(x_sez[:-1]) + 0.5 * (x_sez[1] - x_sez[0])
x_boost =       np.linspace(0.24, 2.54, 12)
bins_x_boost =  np.array(x_boost[:-1]) + 0.5 * (x_boost[1] - x_boost[0])


#Apply hist selection for simpple coalescence
y_simple =  np.histogram(y_simple, x_simple)
y_sez =     np.histogram(y_sez, x_sez)
#y_advance = np.histogram(y_advance, xbins)
y_boost =   np.histogram(y_boost, x_boost)

#x = binning(xbins)

y_simple_err =  (np.sqrt(y_simple[0]) / N_ev) / (x_simple[1]-x_simple[0])
y_sez_err =     (np.sqrt(y_sez[0]) / N_ev) / (x_sez[1]-x_sez[0])
#y_advance_err = (np.sqrt(y_advance[0]) / N_ev) / width(xbins)
y_boost_err =   (np.sqrt(y_boost[0]) / N_ev) / (x_boost[1]-x_boost[0])

weight = (np.ones_like(bins_x_simple) / N_ev) / (x_simple[1]-x_simple[0])


#Plotting
##################################################################################################################

fig, ax = plt.subplots()

make_error_boxes(ax, xbins_U1S_babar, y_U1S_babar, xerr_U1S_babar, yerr_U1S_babar)

#ax.errorbar(xbins_cleo, y_cleo,  yerr=yerr_cleo, fmt='ko', elinewidth=0.7,  label='Experiment (CLEO)')

ax.errorbar(bins_x_simple, y_simple[0]*weight,  yerr=y_simple_err,  fmt='ys', elinewidth=0.7, label='Simple Coalescence')
ax.errorbar(bins_x_sez, y_sez[0]*weight,        yerr=y_sez_err,     fmt='bd', elinewidth=0.7, label='Cross Section Model')
#ax.errorbar(x, y_advance[0]*weight, yerr=y_advance_err, fmt='gv', elinewidth=0.7, label='Advanced Coalescence')
ax.errorbar(bins_x_boost, y_boost[0]*weight,    yerr=y_boost_err,   fmt='r.', elinewidth=0.7, label='Advanced Coalescence')

ax.legend(fontsize=8)

#Decorations
ax.set_title(r'Antideuteron Spectra in $\Upsilon(1S) \rightarrow \bar{d} + X$')
ax.set_ylabel(r'$\frac{dBR}{dp}$ $[Gev^{-1}c]$')
ax.set_xlabel(r'p [GeV/c]')
ax.yaxis.grid(True, linestyle='--')

#ax.set_ylim([0, 4.2e-06])
#ax.set_xlim([0, 2.3])

plt.savefig('Istogrammi_spettro/spettro_U1S_babar_2.jpg', format='jpg', dpi=900)


# In[45]:


#Fit the spectrum in momentum with Maxwell's distribution

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec  # for unequal plot boxes
import scipy.optimize
from scipy import stats

    
#Define deuteron mass
d_m = 1.87561294257  #+- 5.7e-10 GeV
d_m2 = pow(d_m, 2)

#Define the fitting function
def maxwell_distr(x, a, b):
    return a * pow(x,2) * np.exp(-b * np.sqrt(pow(x,2) + d_m2))

#Read spectrum from data (p=p_cm, y=dBR/dp, dy=y uncertainty)
y =  y_simple[0]*weight
dy = y_simple_err

#Initial guesses for fitting parameters
a = 1.0e-02
b = 4.0

#Fit data using SciPy's Levenberg-Marquart Method
nlfit, nlpcov = scipy.optimize.curve_fit(maxwell_distr, bins_x_simple, y, p0=[a, b], sigma=dy)

#Unpack fitting parameters
a, b = nlfit

#Unpack uncertainties in fitting parameters from diagonal of covariance matrix
da, db = [np.sqrt(nlpcov[j,j]) for j in range(nlfit.size)]

#Create fitting function from fitted parameters
x_fit = np.linspace(0.0, 3.0, 100)
y_fit = maxwell_distr(x_fit, a, b)

#Calculate residuals, reduced chi squared and probability
resids = y - maxwell_distr(bins_x_simple, a, b)
chisq =  ((resids/dy)**2).sum()
df =     float(x.size-2)
chisqR = chisq/df
prob =   1 - stats.chi2.cdf(chisq, df)

#Create figure window to plot data
fig = plt.figure(1, figsize=(8,8))
gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

#Top plot: data and fit
ax1 = fig.add_subplot(gs[0])
ax1.plot(x_fit, y_fit, 'xkcd:dark sky blue')
ax1.errorbar(bins_x_simple, y, yerr=dy, fmt='or', ecolor='black', elinewidth=0.7)
ax1.set_title(r'Advanced Coalescence in $\Upsilon(1S) \rightarrow \bar{d} + X$', fontsize=15)
ax1.set_xlabel(r'p [GeV/c]')
ax1.set_ylabel(r'$\frac{dBR}{dp}$ $[Gev^{-1}c]$')

a * pow(x,2) * np.exp(-b * np.sqrt(pow(x,2) + d_m2))

ax1.text(0.66, 0.93, r'$\frac{dBR}{dp} = a \cdot p^2 \cdot exp(-b \sqrt{p^2 +m_d^2})$', transform = ax1.transAxes)
ax1.text(0.66, 0.86, f'a = ({round(a,2)} $\pm$ {round(da,2)})'r' $GeV^{-3}$', transform = ax1.transAxes)
ax1.text(0.66, 0.80, f'b = ({round(b,2)} $\pm$ {round(db,2)})'r' $GeV^{-1}$', transform = ax1.transAxes)
ax1.text(0.66, 0.73, f'$\chi^2_R$ = {round(chisqR,2)}   ($P$ = {"{0:.2%}".format(prob)})', transform = ax1.transAxes)

#Bottom plot: residuals
ax2 = fig.add_subplot(gs[1])
ax2.errorbar(bins_x_simple, resids, yerr=dy, ecolor="black", fmt="ro", elinewidth=0.7)
ax2.axhline(color='xkcd:dark sky blue', zorder=-1)
ax2.set_xlabel(r'p [GeV]')
ax2.set_ylabel(r'$\frac{dBR}{dp}$ $[Gev^{-1}]$')
#ax2.set_ylim(-6e-06, 6e-06)
#ax2.set_yticks((-2, 0, 2))


#Show result
plt.show()

#fig.savefig('Istogrammi_spettro/fit_simple_U1S_babar_2.jpg', format='jpg', dpi=900)


# In[ ]:




