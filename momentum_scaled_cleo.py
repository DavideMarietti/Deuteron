#!/usr/bin/env python
# coding: utf-8

# In[13]:


import numpy as np
import matplotlib.pyplot as plt

bins = np.linspace(0, 1, 11) #10 bins

##################################################################################################################

p_data_enhancement = np.array([3.19021739130435, 2.73913043478261, 2.19021739130435, 1.41847826086956,
                             0.956521739130434, 0.766304347826086, 0.592391304347825, 0.559782608695651, 0, 0])

p_data_error = np.array([3.79891304347826, 2.80978260869565, 2.27717391304348, 1.47282608695652, 
                             1.01086956521739, 0.869565217391303, 0.744565217391303, 0.815217391304347, 
                             0.152173913043477, 0.195652173913043])

p_error = p_data_error - p_data_enhancement

##################################################################################################################

anti_data_enhancement = np.array([2.07285032458065, 2.58505854729056, 2.03061647304369, 1.50952256717975, 
                             0.94952502263624, 0.611746289958718, 0.368425898927273, 0.191771151455625,
                             0.0817667009407472, 0.25510658215804])

anti_data_error = np.array([2.20441988950276, 2.71270718232044, 2.12707182320442, 1.55801104972376, 1, 
                       0.701657458563536, 0.464088397790056, 0.292817679558011, 0.187845303867404, 
                       0.513812154696133])

anti_error = anti_data_error - anti_data_enhancement

##################################################################################################################

lambda_data_enhancement = np.array([2.37096774193548, 3.21505376344086, 2.94623655913978, 2.04838709677419,
                             1.45698924731183, 1.09677419354839, 0.827956989247312, 0, 0, 0])

lambda_data_error = np.array([2.9367572894325, 3.53031708306765, 3.0538954263521, 2.12596371224857, 1.54741754101212,
                       1.23240681699612, 0.997985624021796, 0.84418294591618, 0.265694742333778, 0.316083125615906])

lambda_error = lambda_data_error - lambda_data_enhancement

##################################################################################################################

phi_data_enhancement = np.array([2.04899481146392, 2.02857308031391, 1.51888251362171, 1.06834703935483,
                             0.816739171279502, 0.764044456648986, 0.662961945917822, 0.674798745501583, 0, 0])

phi_data_error = np.array([2.97849462365591, 2.11827956989247, 1.55913978494624, 1.11827956989247, 0.881720430107527,
                       0.844086021505377, 0.790322580645162, 0.849462365591398, 0.763440860215054,
                       0.333333333333333])

phi_error = phi_data_error - phi_data_enhancement

##################################################################################################################

fig, ax = plt.subplots()
new_bins = np.array(bins[:-1]) + 0.5 * (bins[1] - bins[0])

ax.errorbar(new_bins, p_data_enhancement, yerr=p_error, fmt='bs', elinewidth=0.7, barsabove=True, label=r"$ggg \rightarrow p + X$")
ax.errorbar(new_bins, anti_data_enhancement, yerr=anti_error, fmt='g^', elinewidth=0.7, barsabove=True, label=r"$ggg \rightarrow \bar{p} + X$")
ax.errorbar(new_bins, lambda_data_enhancement, yerr=lambda_error, fmt='yd', elinewidth=0.7, barsabove=True, label=r"$ggg \rightarrow \Lambda + X$")
ax.errorbar(new_bins, phi_data_enhancement, yerr=phi_error, fmt='r<', elinewidth=0.7, barsabove=True, label=r"$ggg \rightarrow \Phi + X$")

ax.set_title('CLEO Experimental Data')
ax.set_ylabel(r'$\frac{ggg}{q\bar{q}}$ Enhancement')
ax.set_xlabel(r'Scaled Momentum ($p / E_{beam}$)')

ax.yaxis.grid(True)
ax.legend(fontsize=8)


plt.ylim((-0.5, 4))
ax.hlines(1, 0, 1, color='k')
plt.tight_layout()

plt.savefig('Istogrammi/enhancement_scaled_cleo.jpg', format='jpg', dpi=900)
plt.show()


# In[ ]:




