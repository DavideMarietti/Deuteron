#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np


#lambda
l = ((7466/4934) * (98021/100000) - (0.108)) / 0.817

#p
p = ((21861/13465) * (98021/100000) - (0.108)) / 0.817

#anti-p
pbar = ((21768/13389) * (98021/100000) - (0.108)) / 0.817

#phi
phi = ((3315/4434) * (98021/100000) - (0.108)) / 0.817

MC = np.array([l, p, pbar, phi])

print(MC)


# In[2]:


get_ipython().run_line_magic('matplotlib', 'ipympl')
import numpy as np
import matplotlib.pyplot as plt


particles = np.array([0, 1, 2, 3])

data = np.array([2.668, 1.623, 1.634, 1.423])
err_stat_data = np.array([0.027, 0.014, 0.014, 0.051])
err_sys_data = np.array([0.051, 0.116, 0.111, 0.065])
err_data = np.sqrt(err_stat_data*err_stat_data + err_sys_data*err_sys_data)


labels = ['$\Lambda$', '$p$', r'$\bar{p}$', '$\Phi$']
x_pos = np.arange(len(labels))

fig, ax = plt.subplots()
ax.errorbar(particles, data, yerr=err_data, fmt='r.', elinewidth=0.5)
ax.plot(particles, MC, 'g^')

ax.set_ylabel(r'$\frac{ggg}{q\bar{q}}$ Enhancement')
ax.set_xticks(x_pos)
ax.set_xticklabels(labels, fontsize=15)
ax.set_title('Momentum-Integrated Enhancement for $ggg$ events')
ax.yaxis.grid(True)

ax.legend([r'$\Upsilon(1S)$ MC', r'$\Upsilon(1S)$ Data (CLEO)'], fontsize=10)


plt.ylim((0, 3.5))
plt.hlines(1, -1, 4, color='k')
plt.tight_layout()
plt.savefig('Istogrammi/momentum_integr_enhancement.png')
plt.show()


# In[ ]:





# In[ ]:




