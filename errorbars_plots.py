#!/usr/bin/env python
# coding: utf-8

# In[1]:


## import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D

simple = [75.1, 74.3, 71.6, 73.4, 63.8]
upper_simple = [2.4, 2.8, 2.0, 1.3, 3.2]
lower_simple = [2.3, 2.6, 2.9, 1.5, 3.0]

sez = [1.87, 1.74, 1.57, 1.72, 1.13]
upper_sez = [0.22, 0.19, 0.16, 0.14, 0.14]
lower_sez = [0.17, 0.16, 0.21, 0.13, 0.16]

advance = [1.58, 1.62, 1.68,1.62,1.89]
upper_advance = [0.06, 0.07, 0.09, 0.04, 0.1]
lower_advance = [0.05, 0.06, 0.08, 0.03, 0.09]

x = [9.46, 10.02, 10.35, 10.58]

#------------------------------------------------#
elements = 5
x_1 = [9, 9.5, 10, 10.5, 11]
y_1 = advance[3]-lower_advance[3]
y_2 = advance[3]+upper_advance[3]

y_3 = [y_1] * elements
y_4 = [y_2] * elements

area = plt.fill_between(x_1, y_3, y_4, color='xkcd:bluegrey', alpha=0.15,)
line = plt.hlines(1.62, 9, 11, color='grey', linestyle='--', linewidth=0.8)
Y1 = plt.errorbar(x[0], advance[0], yerr=[[lower_advance[0]],[upper_advance[0]]], fmt='ys', ecolor='black', elinewidth=1)
Y2 = plt.errorbar(x[1], advance[1], yerr=[[lower_advance[1]],[upper_advance[1]]], fmt='bs', fillstyle ='none', ecolor='black', elinewidth=1)
Y3 = plt.errorbar(x[2], advance[2], yerr=[[lower_advance[2]],[upper_advance[2]]], fmt='go', ecolor='black', elinewidth=1)
#plt.errorbar(x[3], simple[3], yerr=[[lower_simple[3]],[upper_simple[3]]], fmt='o', color='xkcd:bluegrey', capsize=3, ecolor='black', elinewidth=1, label=r'$\Upsilon (1S,2S,3S) \rightarrow\bar{d}X$')
continuo = plt.errorbar(x[3], advance[4], yerr=[[lower_advance[4]],[upper_advance[4]]], fmt='ro', fillstyle ='none', ecolor='black', elinewidth=1)


#plt.xticks([])
#plt.title('Advanced Coalescence Model')
plt.xlabel('$\sqrt{s}$ [GeV]')
plt.ylabel(r'$\sigma$ [$b$]')
plt.xlim(9,11)
#plt.grid(False, linestyle='--')
#plt.legend([(area, line), Y1, Y2, Y3, continuo], [r'$\Upsilon (1,2,3S) \rightarrow\bar{d}X$', r'$\Upsilon (1S) \rightarrow\bar{d}X$', r'$\Upsilon (2S) \rightarrow\bar{d}X$', r'$\Upsilon (3S) \rightarrow\bar{d}X$', r'$e^+ e^- \rightarrow\bar{d}X$'])

plt.savefig('./best_fit_advance_E.pdf', format='pdf', dpi=500)
plt.show()


# In[2]:


# sigma = 3.36 fm

import numpy as np
import matplotlib.pyplot as plt

C = 0.71

x = np.array([-0.25, -0.5, -1, -2, -3, -4, -5, -6, -8, -10, -14])

y_1S =  C * np.array([236, 337, 482, 670, 807, 917, 1026, 1120, 1274, 1414, 1676]) * 1e-09
err_1S = np.sqrt(C * np.array([236, 337, 482, 670, 807, 917, 1026, 1120, 1274, 1414, 1676])) * 1e-09
y_2S =  C * np.array([195, 267, 396, 526, 638, 757, 837, 920, 1060, 1162, 1376]) * 1e-09
err_2S = np.sqrt(C * np.array([195, 267, 396, 526, 638, 757, 837, 920, 1060, 1162, 1376])) * 1e-09
y_3S =  C * np.array([171, 258, 365, 515, 609, 706, 776, 850, 976, 1100, 1301]) * 1e-09
err_3S = np.sqrt(C * np.array([171, 258, 365, 515, 609, 706, 776, 850, 976, 1100, 1301])) * 1e-09
y_c =   np.array([65, 94, 136, 187, 222, 263, 291, 331, 374, 411, 485]) / (3.69*1e09)
err_c =  np.sqrt(np.array([65, 94, 136, 187, 222, 263, 291, 331, 374, 411, 485])) / (3.69*1e09)


#Sum of Y(1S) and Y(2S)
summ = C * ((102/260)*np.array([236, 337, 482, 670, 807, 917, 1026, 1120, 1274, 1414, 1676]) + (158/260)*np.array([195, 267, 396, 526, 638, 757, 837, 920, 1060, 1162, 1376]))
y = summ * 1e-09
err_stat_y = np.sqrt(summ) * 1e-09
        

sigma_upper_1S = np.array([18, 22, 30, 32, 44, 49, 51, 64, 72, 67, 87])
sigma_upper_2S = np.array([11, 12, 28, 22, 33, 52, 40, 52, 53, 63, 72])
sigma_upper    = ((102/260)*sigma_upper_1S + (158/260)*sigma_upper_2S) * 1e-09
sigma_lower_1S = np.array([13, 18, 26, 36, 39, 59, 53, 53, 74, 80, 89])
sigma_lower_2S = np.array([13, 19, 8, 23, 48, 40, 50, 46, 40, 75, 87])
sigma_lower    = ((102/260)*sigma_lower_1S + (158/260)*sigma_lower_2S) * 1e-09

err_stat_sigma_upper = np.sqrt(pow(sigma_upper, 2) + pow(err_stat_y, 2))
err_stat_sigma_lower = np.sqrt(pow(sigma_lower, 2) + pow(err_stat_y, 2))


#Experimental limit
x_exp = np.array([-2, -6, -10, -14])
y_exp = np.array([1.5, 0.97, 0.71, 0.63]) * 1e-06
x_exp_err = np.array([2, 2, 2, 2])
y_exp_err = np.array([1.2e-07, 8.271e-08, 2.995e-07, 8.04e-08])


fig = plt.figure(figsize=(9, 5))

plt.errorbar(x_exp, y_exp, xerr=x_exp_err, yerr=y_exp_err, fmt='k.', uplims=True, label='Belle U.L.', elinewidth=2)
plt.errorbar(x, y, yerr=([err_stat_sigma_lower,err_stat_sigma_upper]), fmt='bo', fillstyle='none', ecolor='k', elinewidth=1, label=r'Prediction')

plt.ylabel(r'$\mathcal{B}$($\Upsilon(1S,2S) \rightarrow H X$)$\mathcal{B}$($H \rightarrow \Lambda p \pi^-$)')
plt.xlabel(r'$m_H - 2m_{\Lambda}$ [MeV]')
plt.ylim(2e-08,5e-05)
plt.xlim(-16,2.5)
plt.xticks(np.arange(-16, 4, 2))
#plt.grid(True, linestyle='--')
plt.yscale('log')
plt.tight_layout()
plt.xticks(np.array([2, 0, -2, -4, -6, -8, -10, -12, -14, -16]))


#legend
fontP = FontProperties()
fontP.set_size('medium')
plt.legend( bbox_to_anchor=(1.0, 1), loc='upper left', prop=fontP).get_frame().set_alpha(None)
fig.subplots_adjust(right=0.75)

#lines and shadows
plt.hlines(2.81 * 1e-05, -19, 3, color='k', linestyle='--', label='dds')
plt.text(-12.5, 3 * 1e-05, r'$\mathcal{B}$($\Upsilon(1S) \rightarrow \bar{d} X$)', ha='right', va='bottom')

plt.vlines(0, 1e-08, 1e-03, color='r', linestyle='--', label='dds')
plt.text(-0.2, 2.5 * 1e-05, r'$H \rightarrow \Lambda p \pi^-$', color='r', ha='right', va='top')
plt.text( 1.8, 2.5 * 1e-05, r'$H \rightarrow \Lambda \Lambda$', color='r', ha='right', va='top')

plt.vlines(-7.13, 1e-08, 1e-03, color='c', linestyle='--', label='dds')
plt.text(-7, 3.3 * 1e-05, r'"Nagara event"', color='c', ha='left', va='bottom')

#areas
plt.axvspan(0, 3, alpha=0.3, color='k', hatch='///', fill=False)
plt.axvspan(-7.13, -17, alpha=0.6, color='c', hatch='///', fill=False)
plt.axvspan(-18, -3.1, alpha=0.2, color='grey')

#plt.savefig('./BR_H_12S_nagara_2.pdf', format='pdf', dpi=500)
plt.show()


# In[3]:


y_1S =  C * np.array([236, 337, 482, 670, 807, 917, 1026, 1120, 1274, 1414, 1676]) * 1e-09
err_1S = np.sqrt(C * np.array([236, 337, 482, 670, 807, 917, 1026, 1120, 1274, 1414, 1676])) * 1e-09

#plt.errorbar(x_exp, y_exp, xerr=x_exp_err, yerr=y_exp_err, fmt='k.', uplims=True, label='Belle U.L.', elinewidth=2)
plt.errorbar(x, y_1S, yerr=err_1S, fmt='bo', fillstyle='none', ecolor='k', elinewidth=1, label=r'Prediction')

plt.ylabel(r'$\mathcal{B}$($\Upsilon(1S) \rightarrow H X$)$\mathcal{B}$($H \rightarrow \Lambda p \pi^-$)')
plt.xlabel(r'$m_H - 2m_{\Lambda}$ [MeV]')
plt.ylim(1e-08,6e-05)
plt.xlim(-16,2.5)
plt.xticks(np.arange(-16, 4, 2))
#plt.grid(True, linestyle='--')
plt.yscale('log')
plt.tight_layout()
plt.xticks(np.array([2, 0, -2, -4, -6, -8, -10, -12, -14, -16]))


#lines and shadows
plt.hlines(2.81 * 1e-05, -19, 3, color='k', linestyle='--', label='dds')
plt.text(-15.6, 3.1 * 1e-05, r'$\mathcal{B}$($\Upsilon(1S) \rightarrow \bar{d} X$)', ha='left', va='bottom')

plt.vlines(0, 1e-08, 1e-03, color='r', linestyle='--', label='dds')
plt.text(-0.2, 2.5 * 1e-05, r'$H \rightarrow \Lambda p \pi^-$', color='r', ha='right', va='top')
plt.text( 0.3, 2.5 * 1e-05, r'$H \rightarrow \Lambda \Lambda$', color='r', ha='left', va='top')

plt.vlines(-7.13, 1e-08, 1e-03, color='c', linestyle='--', label='dds')
plt.text(-7, 3.35 * 1e-05, r'"Nagara event"', color='c', ha='left', va='bottom')

#areas
plt.axvspan(0, 3, alpha=0.3, color='k', hatch='///', fill=False)
plt.axvspan(-7.13, -17, alpha=0.6, color='c', hatch='///', fill=False)
plt.axvspan(-18, -3.1, alpha=0.2, color='grey')

plt.savefig('./BR_H_1S_nagara_2.pdf', format='pdf', dpi=500)
plt.show()


# In[4]:


y_2S =  C * np.array([195, 267, 396, 526, 638, 757, 837, 920, 1060, 1162, 1376]) * 1e-09
err_2S = np.sqrt(C * np.array([195, 267, 396, 526, 638, 757, 837, 920, 1060, 1162, 1376])) * 1e-09

plt.errorbar(x, y_2S, yerr=err_2S, fmt='bo', fillstyle='none', ecolor='k', elinewidth=1, label=r'Prediction')

plt.ylabel(r'$\mathcal{B}$($\Upsilon(2S) \rightarrow H X$)$\mathcal{B}$($H \rightarrow \Lambda p \pi^-$)')
plt.xlabel(r'$m_H - 2m_{\Lambda}$ [MeV]')
plt.ylim(1e-08,6e-05)
plt.xlim(-16,2.5)
plt.xticks(np.arange(-16, 4, 2))
#plt.grid(True, linestyle='--')
plt.yscale('log')
plt.tight_layout()
plt.xticks(np.array([2, 0, -2, -4, -6, -8, -10, -12, -14, -16]))


#lines and shadows
plt.hlines(2.64 * 1e-05, -19, 3, color='k', linestyle='--', label='dds')
plt.text(-15.6, 2.95 * 1e-05, r'$\mathcal{B}$($\Upsilon(2S) \rightarrow \bar{d} X$)', ha='left', va='bottom')

plt.vlines(0, 1e-08, 1e-03, color='r', linestyle='--', label='dds')
plt.text(-0.2, 2.3 * 1e-05, r'$H \rightarrow \Lambda p \pi^-$', color='r', ha='right', va='top')
plt.text( 0.3, 2.3 * 1e-05, r'$H \rightarrow \Lambda \Lambda$', color='r', ha='left', va='top')

plt.vlines(-7.13, 1e-08, 1e-03, color='c', linestyle='--', label='dds')
plt.text(-7, 3.35 * 1e-05, r'"Nagara event"', color='c', ha='left', va='bottom')

#areas
plt.axvspan(0, 3, alpha=0.3, color='k', hatch='///', fill=False)
plt.axvspan(-7.13, -17, alpha=0.6, color='c', hatch='///', fill=False)
plt.axvspan(-18, -3.1, alpha=0.2, color='grey')

#plt.savefig('./BR_H_2S_nagara_2.pdf', format='pdf', dpi=500)
plt.show()


# In[5]:


y_3S =  C * np.array([171, 258, 365, 515, 609, 706, 776, 850, 976, 1100, 1301]) * 1e-09
err_3S = np.sqrt(C * np.array([171, 258, 365, 515, 609, 706, 776, 850, 976, 1100, 1301])) * 1e-09

plt.errorbar(x, y_3S, yerr=err_3S, fmt='bo', fillstyle='none', ecolor='k', elinewidth=1, label=r'Prediction')

plt.ylabel(r'$\mathcal{B}$($\Upsilon(3S) \rightarrow H X$)$\mathcal{B}$($H \rightarrow \Lambda p \pi^-$)')
plt.xlabel(r'$m_H - 2m_{\Lambda}$ [MeV]')
plt.ylim(1e-08,6e-05)
plt.xlim(-16,2.5)
plt.xticks(np.arange(-16, 4, 2))
#plt.grid(True, linestyle='--')
plt.yscale('log')
plt.tight_layout()
plt.xticks(np.array([2, 0, -2, -4, -6, -8, -10, -12, -14, -16]))


#lines and shadows
plt.hlines(2.33 * 1e-05, -19, 3, color='k', linestyle='--', label='dds')
plt.text(-15.6, 2.55 * 1e-05, r'$\mathcal{B}$($\Upsilon(3S) \rightarrow \bar{d} X$)', ha='left', va='bottom')

plt.vlines(0, 1e-08, 1e-03, color='r', linestyle='--', label='dds')
plt.text(-0.2, 1.95 * 1e-05, r'$H \rightarrow \Lambda p \pi^-$', color='r', ha='right', va='top')
plt.text( 0.3, 1.95 * 1e-05, r'$H \rightarrow \Lambda \Lambda$', color='r', ha='left', va='top')

plt.vlines(-7.13, 1e-08, 1e-03, color='c', linestyle='--', label='dds')
plt.text(-7, 3.35 * 1e-05, r'"Nagara event"', color='c', ha='left', va='bottom')

#areas
plt.axvspan(0, 3, alpha=0.3, color='k', hatch='///', fill=False)
plt.axvspan(-7.13, -17, alpha=0.6, color='c', hatch='///', fill=False)
plt.axvspan(-18, -3.1, alpha=0.2, color='grey')

#plt.savefig('./BR_H_3_nagara_2.pdf', format='pdf', dpi=500)
plt.show()


# In[6]:


y_c =   np.array([65, 94, 136, 187, 222, 263, 291, 331, 374, 411, 485]) / (3.69*1e09)
err_c =  np.sqrt(np.array([65, 94, 136, 187, 222, 263, 291, 331, 374, 411, 485])) / (3.69*1e09)

plt.errorbar(x, y_c, yerr=err_c, fmt='bo', fillstyle='none', ecolor='k', elinewidth=1, label=r'Prediction')

plt.ylabel(r'$\frac{\sigma(e^+e^- \rightarrow \bar{d} X)}{\sigma(e^+e^- \rightarrow hadrons)}$ $\mathcal{B}$($H \rightarrow \Lambda p \pi^-$)')
plt.xlabel(r'$m_H - 2m_{\Lambda}$ [MeV]')
plt.ylim(1e-08,6e-05)
plt.xlim(-16,2.5)
plt.xticks(np.arange(-16, 4, 2))
#plt.grid(True, linestyle='--')
plt.yscale('log')
plt.tight_layout()
plt.xticks(np.array([2, 0, -2, -4, -6, -8, -10, -12, -14, -16]))


#lines and shadows
plt.hlines(3.01 * 1e-06, -19, 3, color='k', linestyle='--', label='dds')
plt.text(-15.6, 3.1* 1e-06, r'$\frac{\sigma(e^+e^- \rightarrow \bar{d} X)}{\sigma(e^+e^- \rightarrow hadrons)}$', ha='left', va='bottom')

plt.vlines(0, 1e-08, 1e-03, color='r', linestyle='--', label='dds')
plt.text(-0.2, 2.6 * 1e-06, r'$H \rightarrow \Lambda p \pi^-$', color='r', ha='right', va='top')
plt.text( 0.3, 2.6 * 1e-06, r'$H \rightarrow \Lambda \Lambda$', color='r', ha='left', va='top')

plt.vlines(-7.13, 1e-08, 1e-03, color='c', linestyle='--', label='dds')
plt.text(-7, 3.35 * 1e-05, r'"Nagara event"', color='c', ha='left', va='bottom')

#areas
plt.axvspan(0, 3, alpha=0.3, color='k', hatch='///', fill=False)
plt.axvspan(-7.13, -17, alpha=0.6, color='c', hatch='///', fill=False)
plt.axvspan(-18, -3.1, alpha=0.2, color='grey')

#plt.savefig('./BR_H_c_nagara_2.pdf', format='pdf', dpi=500)
plt.show()


# In[7]:


#Original
# function of sigma

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

C = 0.71

x = np.array([-0.25, -0.5, -1, -2, -3, -4, -5, -6, -8, -10, -14])
x_new = np.array([-0.25, -0.5, -1, -2, -3, -4, -5, -6, -8, -10, -14, -16])

# sigma = 2
##################################################################################################################
y_1S_1 =   np.array([670, 917, 1274, 1778, 2121, 2425, 2723, 2952, 3436, 3845, 4503, 4791]) * 1e-09
err_1S_1 = np.sqrt(np.array([670, 917, 1274, 1778, 2121, 2425, 2723, 2952, 3436, 3845, 4503, 4791])) * 1e-09
y_2S_1 =   np.array([525, 756, 1060, 1477, 1841, 2155, 2412, 2608, 2981, 3337, 3928, 4215]) * 1e-09
err_2S_1 = np.sqrt(np.array([525, 756, 1060, 1477, 1841, 2155, 2412, 2608, 2981, 3337 ,3928, 4215])) * 1e-09

summ_1_new = C * ((102/260)*np.array([670, 917, 1274, 1778, 2121, 2425, 2723, 2952, 3436, 3845, 4503, 4791]) + (158/260)*np.array([525, 756, 1060, 1477, 1841, 2155, 2412, 2608, 2981, 3337, 3928, 4215]))
summ_1 = C * ((102/260)*np.array([670, 917, 1274, 1778, 2121, 2425, 2723, 2952, 3436, 3845, 4503]) + (158/260)*np.array([525, 756, 1060, 1477, 1841, 2155, 2412, 2608, 2981, 3337, 3928]))
y_1 = summ_1 * 1e-09
y_1_new = summ_1_new * 1e-09
err_y_1 = np.sqrt(summ_1) * 1e-09


# sigma = 2.5
##################################################################################################################
y_1S_2 =   np.array([463, 606, 838, 1161, 1400, 1619, 1810, 1937, 2221, 2453, 2896]) * 1e-09
err_1S_2 = np.sqrt(np.array([463, 606, 838, 1161, 1400, 1619, 1810, 1937, 2221, 2453, 2896])) * 1e-09
y_2S_2 =   np.array([355, 484, 678, 958, 1149, 1333, 1494, 1657, 1930, 2176, 2548]) * 1e-09
err_2S_2 = np.sqrt(np.array([355, 484, 678, 958, 1149, 1333, 1494, 1657, 1930, 2176, 2548])) * 1e-09

summ_2 = C * ((102/260)*np.array([463, 606, 838, 1161, 1400, 1619, 1810, 1937, 2221, 2453, 2896]) + (158/260)*np.array([355, 484, 678, 958, 1149, 1333, 1494, 1657, 1930, 2176, 2548]))
y_2 = summ_2 * 1e-09
err_y_2 = np.sqrt(summ_2) * 1e-09


# sigma = 3.36
##################################################################################################################
y_1S_3 =   np.array([236, 337, 482, 670, 807, 917, 1026, 1120, 1274, 1414, 1676]) * 1e-09
err_1S_3 = np.sqrt(np.array([236, 337, 482, 670, 807, 917, 1026, 1120, 1274, 1414, 1676])) * 1e-09
y_2S_3 =   np.array([195, 267, 396, 526, 638, 757, 837, 920, 1060, 1162, 1376]) * 1e-09
err_2S_3 = np.sqrt(np.array([195, 267, 396, 526, 638, 757, 837, 920, 1060, 1162, 1376])) * 1e-09

summ_3 = C * ((102/260)*np.array([236, 337, 482, 670, 807, 917, 1026, 1120, 1274, 1414, 1676]) + (158/260)*np.array([195, 267, 396, 526, 638, 757, 837, 920, 1060, 1162, 1376]))
y_3 = summ_3 * 1e-09
err_y_3 = np.sqrt(summ_3) * 1e-09


# sigma = 5
##################################################################################################################
y_1S_5 =   np.array([120, 166, 213, 300, 378, 436, 484, 532, 606, 673, 791, 838]) * 1e-09
err_1S_5 = np.sqrt(np.array([120, 166, 213, 300, 378, 436, 484, 532, 606, 673, 791, 838])) * 1e-09
y_2S_5 =   np.array([100, 137, 177, 247, 307, 355, 398, 422, 484, 529, 626, 678]) * 1e-09
err_2S_5 = np.sqrt(np.array([100, 137, 177, 247, 307, 355, 398, 422, 484, 529, 626, 678])) * 1e-09

summ_5_new = C * ((102/260)*np.array([120, 166, 213, 300, 378, 436, 484, 532, 606, 673, 791, 838]) + (158/260)*np.array([100, 137, 177, 247, 307, 355, 398, 422, 484, 529, 626, 678]))
summ_5 = C * ((102/260)*np.array([120, 166, 213, 300, 378, 436, 484, 532, 606, 673, 791]) + (158/260)*np.array([100, 137, 177, 247, 307, 355, 398, 422, 484, 529, 626]))
y_5_new = summ_5_new * 1e-09
y_5 = summ_5 * 1e-09
err_y_5 = np.sqrt(summ_5) * 1e-09
        

# sigma = 8
##################################################################################################################
y_1S_8 =   np.array([52, 70, 99, 133, 161, 175, 191, 203, 234, 267, 312]) * 1e-09
err_1S_8 = np.sqrt(np.array([52, 70, 99, 133, 161, 175, 191, 203, 234, 267, 312])) * 1e-09
y_2S_8 =   np.array([39, 55, 79, 110, 132, 148, 163, 170, 194, 217, 252]) * 1e-09
err_2S_8 = np.sqrt(np.array([39, 55, 79, 110, 132, 148, 163, 170, 194, 217, 252])) * 1e-09

summ_8 = C * ((102/260)*np.array([52, 70, 99, 133, 161, 175, 191, 203, 234, 267, 312]) + (158/260)*np.array([39, 55, 79, 110, 132, 148, 163, 170, 194, 217, 252]))
y_8 = summ_8 * 1e-09
err_y_8 = np.sqrt(summ_8) * 1e-09


#Experimental limit
##################################################################################################################
x_exp = np.array([-2, -6, -10, -14])
y_exp = np.array([1.5, 0.97, 0.71, 0.63]) * 1e-06
x_exp_err = np.array([2, 2, 2, 2])
y_exp_err = np.array([1.2e-07, 8.271e-08, 2.995e-07, 8.04e-08])


# Forbidden area
##################################################################################################################
x_forb = np.array([-16, -14, -8.7, -7.2,  -5.57, -4.44, -3.87, -3.1, -2.17, -1.4, -0.97, -0.71, -0.54, -0.25])
y_forb = C * np.array([7933, 6584, 3294.4, 2511.9, 1706.8, 1213.5, 1004.5, 717.5, 439.5, 225.9, 146.2, 93.1, 61.9, 16]) * 1e-09


fig = plt.figure(figsize=(9, 5))

plt.errorbar(x_exp, y_exp, xerr=x_exp_err, yerr=y_exp_err, fmt='k.', uplims=True, label='Belle U.L.', elinewidth=2)
plt.errorbar(x, y_1, yerr=err_y_1,  ls='none',  marker='o', color='tab:orange', ecolor='black', elinewidth=1, label=r'$\sigma$ = 2.0 fm')
#plt.errorbar(x, y_2, yerr=err_y_2,  ls='none', marker='v', color='tab:green', capsize=5, ecolor='black', elinewidth=0.5, label=r'$\sigma$ = 2.5 fm')
plt.errorbar(x, y_3, yerr=err_y_3, fmt='bo', fillstyle='none', ecolor='black', elinewidth=1, label=r'$\sigma$ = 3.36 fm')
plt.errorbar(x, y_5, yerr=err_y_5,  ls='none', marker='s', color='gold', ecolor='black', elinewidth=1, label=r'$\sigma$ = 5.0 fm')
#plt.errorbar(x, y_8, yerr=err_y_8,  ls='none', marker='>',  color='tab:cyan', capsize=5, ecolor='black', elinewidth=0.5, label=r'$\sigma$ = 8.0 fm')


plt.ylabel(r'$\mathcal{B}$($\Upsilon(1S,2S) \rightarrow H X$)$\mathcal{B}$($H \rightarrow \Lambda p \pi^-$)')
plt.xlabel(r'$m_H - 2m_{\Lambda}$ [MeV]')
plt.ylim(2e-08,5e-05)
plt.xlim(-16,2.5)
plt.xticks(np.arange(-16, 4, 2))
#plt.grid(True, linestyle='--')
plt.yscale('log')
plt.tight_layout()

#legend
fontP = FontProperties()
fontP.set_size('medium')
plt.legend( bbox_to_anchor=(1.0, 1), loc='upper left', prop=fontP).get_frame().set_alpha(None)
fig.subplots_adjust(right=0.75)

#lines and shadows
plt.hlines(2.81 * 1e-05, -19, 3, color='k', linestyle='--', label='dds')
plt.text(-12.5, 3 * 1e-05, r'$\mathcal{B}$($\Upsilon(1S) \rightarrow \bar{d} X$)', ha='right', va='bottom')

plt.vlines(0, 1e-08, 1e-03, color='r', linestyle='--', label='dds')
plt.text(-0.2, 2.5 * 1e-05, r'$H \rightarrow \Lambda p \pi^-$', color='r', ha='right', va='top')
plt.text( 1.8, 2.5 * 1e-05, r'$H \rightarrow \Lambda \Lambda$', color='r', ha='right', va='top')

plt.vlines(-7.13, 1e-08, 1e-03, color='c', linestyle='--', label='dds')
plt.text(-7, 3.3 * 1e-05, r'"Nagara event"', color='c', ha='left', va='bottom')

plt.axvspan(0, 3, alpha=0.3, color='k', hatch='///', fill=False)
plt.axvspan(-7.13, -17, alpha=0.6, color='c', hatch='///', fill=False)
#plt.axvspan(-18, -3.1, alpha=0.2, color='grey')
plt.fill_between(x_forb, 0, y_forb, alpha=0.2, color='grey')

plt.fill_between(x_new, y_1_new, y_5_new, alpha=0.15, color='tab:green')


fig.savefig('./BR_sigma_nagara_2.pdf', format='pdf', dpi=300)
plt.show()


# In[16]:


# function of sigma

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

C = 0.71

x = np.array([-0.25, -0.5, -1, -2, -3, -4, -5, -6, -8, -10, -14])
x_new = np.array([-0.25, -0.5, -1, -2, -3, -4, -5, -6, -8, -10, -14, -16])
x_limit1 = np.array([0, -4])
x_limit2 = np.array([-4, -8])
x_limit3 = np.array([-8, -12])
x_limit4 = np.array([-12, -16])

y_limit = np.array([1, 1])

# sigma = 2
##################################################################################################################
y_1S_1 =   np.array([670, 917, 1274, 1778, 2121, 2425, 2723, 2952, 3436, 3845, 4503, 4791]) * 1e-09
err_1S_1 = np.sqrt(np.array([670, 917, 1274, 1778, 2121, 2425, 2723, 2952, 3436, 3845, 4503, 4791])) * 1e-09
y_2S_1 =   np.array([525, 756, 1060, 1477, 1841, 2155, 2412, 2608, 2981, 3337, 3928, 4215]) * 1e-09
err_2S_1 = np.sqrt(np.array([525, 756, 1060, 1477, 1841, 2155, 2412, 2608, 2981, 3337 ,3928, 4215])) * 1e-09

summ_1_new = C * ((102/260)*np.array([670, 917, 1274, 1778, 2121, 2425, 2723, 2952, 3436, 3845, 4503, 4791]) + (158/260)*np.array([525, 756, 1060, 1477, 1841, 2155, 2412, 2608, 2981, 3337, 3928, 4215]))
summ_1 = C * ((102/260)*np.array([670, 917, 1274, 1778, 2121, 2425, 2723, 2952, 3436, 3845, 4503]) + (158/260)*np.array([525, 756, 1060, 1477, 1841, 2155, 2412, 2608, 2981, 3337, 3928]))
y_1 = summ_1 * 1e-09
y_1_new = summ_1_new * 1e-09
err_y_1 = np.sqrt(summ_1) * 1e-09


# sigma = 2.5
##################################################################################################################
y_1S_2 =   np.array([463, 606, 838, 1161, 1400, 1619, 1810, 1937, 2221, 2453, 2896]) * 1e-09
err_1S_2 = np.sqrt(np.array([463, 606, 838, 1161, 1400, 1619, 1810, 1937, 2221, 2453, 2896])) * 1e-09
y_2S_2 =   np.array([355, 484, 678, 958, 1149, 1333, 1494, 1657, 1930, 2176, 2548]) * 1e-09
err_2S_2 = np.sqrt(np.array([355, 484, 678, 958, 1149, 1333, 1494, 1657, 1930, 2176, 2548])) * 1e-09

summ_2 = C * ((102/260)*np.array([463, 606, 838, 1161, 1400, 1619, 1810, 1937, 2221, 2453, 2896]) + (158/260)*np.array([355, 484, 678, 958, 1149, 1333, 1494, 1657, 1930, 2176, 2548]))
y_2 = summ_2 * 1e-09
err_y_2 = np.sqrt(summ_2) * 1e-09


# sigma = 3.36
##################################################################################################################
y_1S_3 =   np.array([236, 337, 482, 670, 807, 917, 1026, 1120, 1274, 1414, 1676]) * 1e-09
err_1S_3 = np.sqrt(np.array([236, 337, 482, 670, 807, 917, 1026, 1120, 1274, 1414, 1676])) * 1e-09
y_2S_3 =   np.array([195, 267, 396, 526, 638, 757, 837, 920, 1060, 1162, 1376]) * 1e-09
err_2S_3 = np.sqrt(np.array([195, 267, 396, 526, 638, 757, 837, 920, 1060, 1162, 1376])) * 1e-09

summ_3 = C * ((102/260)*np.array([236, 337, 482, 670, 807, 917, 1026, 1120, 1274, 1414, 1676]) + (158/260)*np.array([195, 267, 396, 526, 638, 757, 837, 920, 1060, 1162, 1376]))
y_3 = summ_3 * 1e-09
err_y_3 = np.sqrt(summ_3) * 1e-09


# sigma = 5
##################################################################################################################
y_1S_5 =   np.array([120, 166, 213, 300, 378, 436, 484, 532, 606, 673, 791, 838]) * 1e-09
err_1S_5 = np.sqrt(np.array([120, 166, 213, 300, 378, 436, 484, 532, 606, 673, 791, 838])) * 1e-09
y_2S_5 =   np.array([100, 137, 177, 247, 307, 355, 398, 422, 484, 529, 626, 678]) * 1e-09
err_2S_5 = np.sqrt(np.array([100, 137, 177, 247, 307, 355, 398, 422, 484, 529, 626, 678])) * 1e-09

summ_5_new = C * ((102/260)*np.array([120, 166, 213, 300, 378, 436, 484, 532, 606, 673, 791, 838]) + (158/260)*np.array([100, 137, 177, 247, 307, 355, 398, 422, 484, 529, 626, 678]))
summ_5 = C * ((102/260)*np.array([120, 166, 213, 300, 378, 436, 484, 532, 606, 673, 791]) + (158/260)*np.array([100, 137, 177, 247, 307, 355, 398, 422, 484, 529, 626]))
y_5_new = summ_5_new * 1e-09
y_5 = summ_5 * 1e-09
err_y_5 = np.sqrt(summ_5) * 1e-09
        

# sigma = 8
##################################################################################################################
y_1S_8 =   np.array([52, 70, 99, 133, 161, 175, 191, 203, 234, 267, 312]) * 1e-09
err_1S_8 = np.sqrt(np.array([52, 70, 99, 133, 161, 175, 191, 203, 234, 267, 312])) * 1e-09
y_2S_8 =   np.array([39, 55, 79, 110, 132, 148, 163, 170, 194, 217, 252]) * 1e-09
err_2S_8 = np.sqrt(np.array([39, 55, 79, 110, 132, 148, 163, 170, 194, 217, 252])) * 1e-09

summ_8 = C * ((102/260)*np.array([52, 70, 99, 133, 161, 175, 191, 203, 234, 267, 312]) + (158/260)*np.array([39, 55, 79, 110, 132, 148, 163, 170, 194, 217, 252]))
y_8 = summ_8 * 1e-09
err_y_8 = np.sqrt(summ_8) * 1e-09


#Experimental limit
##################################################################################################################
x_exp = np.array([-2, -6, -10, -14])
y_exp = np.array([1.5, 0.97, 0.71, 0.63]) * 1e-06
x_exp_err = np.array([2, 2, 2, 2])
y_exp_err = np.array([1.2e-07, 8.271e-08, 2.995e-07, 8.04e-08])


# Forbidden area
##################################################################################################################
x_forb = np.array([-16, -14, -8.7, -7.2,  -5.57, -4.44, -3.87, -3.1, -2.17, -1.4, -0.97, -0.71, -0.54, -0.25])
y_forb = C * np.array([7933, 6584, 3294.4, 2511.9, 1706.8, 1213.5, 1004.5, 717.5, 439.5, 225.9, 146.2, 93.1, 61.9, 16]) * 1e-09


#fig = plt.figure(figsize=(9, 5))

plt.errorbar(x_exp, y_exp, xerr=x_exp_err, fmt='k.', label='Belle U.L.', elinewidth=2)
plt.plot(x, y_1, color='tab:green', label=r'$\sigma$ range limits')
#plt.errorbar(x, y_2, yerr=err_y_2,  ls='none', marker='v', color='tab:green', capsize=5, ecolor='black', elinewidth=0.5, label=r'$\sigma$ = 2.5 fm')
plt.errorbar(x, y_3, yerr=err_y_3, fmt='bo', fillstyle='none', ecolor='blue', elinewidth=1, label=r'$\sigma$ = 3.36 fm')
plt.plot(x, y_5, color='tab:green')
#plt.errorbar(x, y_8, yerr=err_y_8,  ls='none', marker='>',  color='tab:cyan', capsize=5, ecolor='black', elinewidth=0.5, label=r'$\sigma$ = 8.0 fm')


plt.ylabel(r'$\mathcal{B}$($\Upsilon(1S,2S) \rightarrow H X$)$\mathcal{B}$($H \rightarrow \Lambda p \pi^-$)')
plt.xlabel(r'$m_H - 2m_{\Lambda}$ [MeV]')
plt.ylim(2e-08,5e-05)
plt.xlim(-16,0)
plt.xticks(np.arange(-16, 2, 2))
#plt.grid(True, linestyle='--')
plt.yscale('log')
plt.tight_layout()

#legend
fontP = FontProperties()
fontP.set_size('small')
plt.legend(loc='upper right', prop=fontP).get_frame().set_alpha(None)

#lines and shadows
plt.hlines(2.81 * 1e-05, -19, 3, color='k', linestyle='--', label='dds')
plt.text(-12.5, 3.1 * 1e-05, r'$\mathcal{B}$($\Upsilon(1S) \rightarrow \bar{d} X$)', ha='right', va='bottom')

plt.vlines(-7.13, 1e-08, 1e-03, color='c', linestyle='--', label='dds')
plt.text(-7, 2.4 * 1e-08, r'"Nagara event"', color='c', ha='left', va='bottom')

#plt.axvspan(0, 3, alpha=0.3, color='k', hatch='///', fill=False)
plt.axvspan(-7.13, -17, alpha=0.5, color='c', hatch='//', fill=False)
#plt.axvspan(-18, -3.1, alpha=0.2, color='grey')
#plt.fill_between(x_forb, 0, y_forb, alpha=0.2, color='grey')
#plt.fill_between(x_new, y_1_new, y_5_new, alpha=0.15, color='tab:green')
plt.fill_between(x_limit1, [y_exp[0],y_exp[0]], y_limit, alpha=0.2, hatch='\\\\', facecolor="none", edgecolor="lightgrey", linewidth=0.0)
plt.fill_between(x_limit2, [y_exp[1],y_exp[1]], y_limit, alpha=0.2, hatch='\\\\', facecolor="none", edgecolor="lightgrey", linewidth=0.0)
plt.fill_between(x_limit3, [y_exp[2],y_exp[2]], y_limit, alpha=0.2, hatch='\\\\', facecolor="none", edgecolor="lightgrey", linewidth=0.0)
plt.fill_between(x_limit4, [y_exp[3],y_exp[3]], y_limit, alpha=0.2, hatch='\\\\', facecolor="none", edgecolor="lightgrey", linewidth=0.0)


plt.savefig('./BR_sigma_nagara_2.pdf', format='pdf', dpi=500)
plt.show()
print( y_3)


# In[9]:


import numpy as np
import matplotlib.pyplot as plt

sigma = np.array([0.2,0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 11])
B = 1 / (1115 * pow(sigma/ 197, 2))

plt.plot(sigma,B,'bo')

plt.ylabel(r'$B_H,_{max}$ [MeV]')
plt.xlabel(r'$\sigma$ [fm]')
#plt.yscale('log')
plt.grid(axis='y', linestyle='--')
plt.xticks(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]))
plt.yticks([0, 1, 2, 3, 4, 6, 8, 10])

plt.ylim(0,10)
plt.xlim(0,10.5)
plt.fill_between(sigma, B, 100, alpha=0.3, color='grey')

#plt.savefig('./B_max.jpg', format='jpg', dpi=300)
plt.show()
print(B)


# In[10]:


import numpy as np
import matplotlib.pyplot as plt

B = np.array([0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])

sigma = (1 / np.sqrt(1115 * B)) * 197
r = (pow((4 * B * pow(1115,2)) / (1115*2), -0.5)) * 197

plt.plot(B, sigma, 'bo', label=r'$\sigma_{max}$ $(B_H)$      ')
plt.xlabel(r'$B_H$ [MeV]')
plt.ylabel(r'$\sigma_{max}$ [fm]')
plt.ylim(0,10)
plt.xlim(0, 10.5)
plt.legend(loc='upper center')
plt.fill_between(B, sigma, 100, alpha=0.3, color='grey')
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
plt.yticks([0,1, 2, 3, 4, 6, 8, 10])
plt.grid(axis='y', linestyle='--')

plt.twinx()

plt.plot(B, r, 'rs', label=r'$r_{RMS}$ $(B_H)$     ')
plt.xlabel(r'$B_H$ [MeV]')
plt.ylabel(r'$r_{RMS}$ [fm]')
plt.ylim(0,10)
plt.xlim(0, 10.5)
plt.yticks([1, 2, 3, 4, 6, 8, 10])
plt.legend()

print(sigma)
#plt.savefig('./sigma_max.jpg', format='jpg', dpi=300plt.savefig('./sigma_max.jpg', format='jpg', dpi=300#)
plt.show()


# In[11]:


simple = np.array([16.4, 11.0, 11.6, 2.7]) * 1e-04
err_simple = np.array([2.4, 2.0, 2.2, 1.6]) * 1e-04
x_sez = np.array([0, 1, 2, 3])

sez = np.array([11.0, 10.6, 8.6, 3.6]) * 1e-04
err_sez = np.array([2.0, 2.0, 1.9, 1.8]) * 1e-04
x_simple = np.array([-0.1, 0.9, 1.9, 2.9])

advance = np.array([13.2, 12.1, 10.3, 5.4]) * 1e-04
err_advance = np.array([2.2, 2.1, 2.1, 2.2]) * 1e-04
x_advance = np.array([0.1, 1.1, 2.1, 3.1])

label = [r'$\frac{\mathcal{B}[\Upsilon(1S)\rightarrow d\bar{d} X]}{\mathcal{B}[\Upsilon(1S)\rightarrow \bar{d} X]}$', r'$\frac{\mathcal{B}[\Upsilon(2S)\rightarrow d\bar{d} X]}{\mathcal{B}[\Upsilon(2S)\rightarrow \bar{d} X]}$', r'$\frac{\mathcal{B}[\Upsilon(3S)\rightarrow d\bar{d} X]}{\mathcal{B}[\Upsilon(3S)\rightarrow \bar{d} X]}$', r'$\frac{\sigma(e^+e^-\rightarrow d\bar{d} X)\, \sigma(e^+e^- \rightarrow {\rm hadrons})}{\sigma(e^+e^- \rightarrow \bar{d} X)}$']

fig = plt.figure(figsize=(9, 5))

plt.errorbar(x_simple, simple, yerr=err_simple, fmt='bs',  fillstyle='none', ecolor='black', elinewidth=1, label='Simple Coalescence')
plt.errorbar(x_sez, sez, yerr=err_sez, fmt='go', ecolor='black', elinewidth=1, label='Cross Section Model')
plt.errorbar(x_advance, advance, yerr=err_advance, fmt='ro', fillstyle='none', ecolor='black', elinewidth=1, label='Advanced Coalescence')

plt.xticks(x_simple, label, fontsize='large')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#plt.title('Advanced Coalescence Model')
#plt.xlabel('$\sqrt{s}$ [GeV]')
plt.ylabel('Rates ratio')
plt.xlim(-0.5,3.5)
plt.ylim(0,2.25e-3)
#plt.grid(axis = 'y', linestyle='--')
plt.legend()

plt.savefig('./ass_prod.pdf', format='pdf', dpi=00)
plt.show()


# In[12]:


import numpy as np
import matplotlib.pyplot as plt

my = np.array([75.1, 63.8])
err_up_my = np.array([2.4, 3.2])
err_low_my = np.array([2.3, 3.0])
x_my = np.array([0.2, 1.05])

paper0 = np.array([65, 1])
x_paper0 = np.array([-0.2, 0.8])

paper1 = np.array([133.1,145.1])*0.5
err_up_paper1 = np.array([3.3, 4.0])*0.5
err_low_paper1 = np.array([3.3, 4.3])*0.5
x_paper1 = np.array([0.95, 0.1])

paper2 = np.array([140, 1])*0.5
err_up_paper2 = np.array([163-140,1])*0.5
err_low_paper2 = np.array([140-105, 1])*0.5
x_paper2 = np.array([-0.1, 0.9])

paper3 = np.array([133, 1])*0.5
err_up_paper3 = np.array([163-140,1])*0.5
err_low_paper3 = np.array([140-105, 1])*0.5
x_paper3 = np.array([0., 1.])

x = np.array([0, 1])

label = [r'$\mathcal{B}[\Upsilon(1S)\rightarrow \bar{d} X]$',  r'$\frac{\sigma(e^+e^- \rightarrow \bar{d} X)}{\sigma(e^+e^- \rightarrow {\rm hadrons})}$']

#fig = plt.figure(figsize=(9, 5))
plt.errorbar(x_paper1, paper1, yerr=[err_low_paper1,err_up_paper1], fmt='bs', ecolor='black', elinewidth=1, label='Dal and Raklev')
plt.errorbar(x_paper2, paper2, yerr=[err_low_paper2,err_up_paper2], fmt='yo', ecolor='black', elinewidth=1, label='Artoisenet and Braaten')
plt.errorbar(x_my, my, yerr=[err_low_my,err_up_my], fmt='go', fillstyle='none', ecolor='black', elinewidth=1, label='Our Results')
plt.plot(x_paper0, paper0, 'rs', fillstyle='none', label='Gustafson and Häkkinen')
plt.plot(x_paper3, paper3, 'k<', fillstyle='none', label='Wild and Ibarra')

plt.xticks(x, label, fontsize='large')
#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.xlim(-0.5,1.5)
plt.ylim(50, 90)
#plt.title('Simple Coalescence Model')
plt.ylabel(r'$k_{cut}$ [MeV]')
#plt.grid(axis = 'y', linestyle='--')
plt.legend()

plt.savefig('./comparison_simple.pdf', format='pdf', dpi=500)
plt.show()


# In[13]:


import numpy as np
import matplotlib.pyplot as plt

my = np.array([1.87, 1.13])
err_up_my = np.array([0.22, 0.14])
err_low_my = np.array([0.17, 0.16])
x_my = np.array([0.05, 1.05])

paper1 = np.array([1.25, 1.13])
x_paper1 = np.array([-0.05, 0.95])

x = np.array([0, 1])

label = [r'$\mathcal{B}[\Upsilon(1S)\rightarrow \bar{d} X]$',  r'$\frac{\sigma(e^+e^- \rightarrow \bar{d} X)}{\sigma(e^+e^- \rightarrow {\rm hadrons})}$']

#fig = plt.figure(figsize=(9, 5))
plt.plot(x_paper1, paper1, 'bs', label='Dal and Raklev')
plt.errorbar(x_my, my, yerr=[err_low_my,err_up_my], fmt='go', fillstyle='none', ecolor='black', elinewidth=1, label='Our Results')


plt.xticks(x, label, fontsize='large')
#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.xlim(-0.5,1.5)
#plt.ylim(50, 90)
#plt.title('Cross Section based Model')
plt.ylabel(r'$1/\sigma_0$ [b$^{-1}$]')
#plt.grid(axis = 'y', linestyle='--')
#plt.legend()

plt.savefig('./comparison_sez.pdf', format='pdf', dpi=500)
plt.show()


# In[14]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

figure(figsize=(10, 5), dpi=300)

x = np.array([133.1, 135.1, 152.1, 193.3, 226.6, 208.9])
upper_x = np.array([4.0, 3.3, 6.3, 18.6, 13.0, 5.3])
lower_x = np.array([4.3, 3.3, 7.3, 23.6, 14.6, 5.7])

y = ['9.46 GeV', '10.58 GeV', '53 GeV', '91.2 GeV', '318 GeV', '7 TeV']
labels = [r'CLEO ($\Upsilon$ decay)', 'BaBar ($e^+e^-$)', 'ISR (pp)', 'ALEPH (Z decay)', 'ZEUS ($e^-$p)', 'ALICE (pp)']


plt.errorbar(x, labels, xerr=[lower_x,upper_x], fmt='ro', capsize=3, ecolor='black', elinewidth=1)

plt.ylim(-1,6)
plt.xlim(100,260)
plt.grid(axis='x', linestyle='--')
plt.xlabel(r'Coalescence momentum  $2 k_{cut}$ [MeV]')


plt.twinx()
plt.errorbar(x, y, xerr=[lower_x,upper_x], fmt='ro', capsize=3, ecolor='black', elinewidth=1)

plt.ylim(-1,6)
plt.xlim(100,270)
plt.ylabel('$\sqrt{s}$ [GeV]')
plt.grid(axis='x', linestyle='--')


#plt.savefig('./best_fit_comparison.jpg', format='jpg', dpi=300)
plt.show()


# In[15]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

figure(figsize=(10, 5), dpi=300)

x = np.array([133.1, 135.1, 152.1, 193.3, 226.6, 208.9])
upper_x = np.array([4.0, 3.3, 6.3, 18.6, 13.0, 5.3])
lower_x = np.array([4.3, 3.3, 7.3, 23.6, 14.6, 5.7])

y = np.array(['9.46 GeV', '10.58 GeV', '53 GeV', '91.2 GeV', '318 GeV', '7 TeV'])
labels = np.array([r'CLEO ($\Upsilon(1S)$)', 'BaBar ($e^+e^-$)', 'ISR (pp)', 'ALEPH (Z decay)', 'ZEUSS ($e^-$p)', 'ALICE (pp)'])

Y = np.array([9.46, 10.58, 53, 91.2, 318, 7000])

#plt.errorbar(x, labels, xerr=[lower_x,upper_x], fmt='ro', capsize=3, ecolor='black', elinewidth=1)
plt.errorbar(x, Y, xerr=[lower_x,upper_x], fmt='ro', capsize=3, ecolor='black', elinewidth=1)

plt.grid(axis='x', linestyle='--')
plt.xlabel(r'Coalescence momentum $p_0$ [MeV]')
plt.xlim(100,270)
plt.ylim(3,10000)

plt.ylabel('$\sqrt{s}$ [GeV]')
plt.grid(axis='x', linestyle='--')
plt.yscale('log')

for label, a, b in zip(labels[1:4], x[1:4], Y[1:4]):
    plt.annotate(label, xy=(a, b), xytext=(-20, 20),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='grey', alpha=0.1),
        arrowprops=dict(arrowstyle = '->', connectionstyle='angle,angleA=180,angleB=90,rad=0', lw=0.5))
    
labels_new = np.append(labels[0],labels[4])
labels_new = np.append(labels_new,labels[5])
x_new = np.append(x[0],x[4])
x_new = np.append(x_new,x[5])
Y_new = np.append(Y[0],Y[4])
Y_new = np.append(Y_new,Y[5])

for label, a, b in zip(labels_new, x_new, Y_new):
    plt.annotate(label, xy=(a, b), xytext=(20, -20),
        textcoords='offset points', ha='left', va='top',
        bbox=dict(boxstyle='round,pad=0.5', fc='grey', alpha=0.1),
        arrowprops=dict(arrowstyle = '->', connectionstyle='angle,angleA=180,angleB=90,rad=3', lw=0.5))
    
#plt.savefig('./best_fit_comparison_log.jpg', format='jpg', dpi=300)
plt.show()

