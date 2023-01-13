import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator

rc('mathtext', default='regular')
#from matplotlib.ticker import ScalarFormatter
#data
time = (482.7932,
        841.7043,
        1112.0819)
#success_rate_d = (0.71,0.63,0.59,0.48)
ord = (1.720791529,
       1.791327333,
       1.842339066)
seq_ord = (1.7133,
           1.7893,
           1.8444)
apd = (1.319918648,
       1.344098549,
       1.372871571)
ekf = (1.7133,
       1.7893,
       1.8444)
#put labels to all data points

# Create a subplots and twin axes
f, ax1 = plt.subplots()
#f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
# zoom-in / limit the view to different portions of the data
#ax2.set_ylim(0.5,0.97)  # outliers only
ax1.set_ylim(1.2, 1.9)  # most of the data
# plot the same data on both axes
#ax1.spines['bottom'].set_visible(False)
#ax2.spines['top'].set_visible(False)
#ax1.plot(time,ekf)
ax1.plot(time, ord, color="black")
#ax1.plot(time, seq)
ax1.plot(time, apd,linestyle='dashed',color="black")
ax1.plot(time, seq_ord,linestyle='dashed',color="black")
#ax1.plot(noc, ekf)

#plt.xscale('log', basex=10)
#limits

#ax1.set_xlim(0.1,81)

ax1.grid(color='gray', alpha=0.5, linestyle='dashed', linewidth=0.5)
#ax2.grid(color='gray', alpha=0.5, linestyle='dashed', linewidth=0.5)
#ax2.grid(color='gray', alpha=0.5, linestyle='dashed', linewidth=0.5)
#saving and showing fig

#ax1.set_xticks(noc)
#ax1.set_xticklabels([1,5,10,30,50,80], fontsize=12)
#ax1.set_yticks(arrayY)
#ax1.set_yticklabels([0.3,0.39,0.47,0.5,0.6,0.65,0.72,0.76,1],fontsize=12)

#legend
lns1 = ax1.plot(time, ord, color="black", label='Only recent data (15k UEs)')
lns2 = ax1.plot(time, seq_ord, color="black", marker="o", markersize=7, label='Sequential (only recent data) (15k UEs)')
#lns3 = ax1.plot(time, seq, color="red", marker="o", markersize=7, label='Sequential (only recent data) (5000 UEs)')
lns4 = ax1.plot(time, apd, color="black", marker="s", markersize=7, label='All previous data (15k UEs)',linestyle='dashed')
#lns5 = ax1.plot(time, ekf, color="black", marker="*", markersize=7, label='EKF (10k UEs)')
#lnsv = plt.axvline(x=time[0],label='localization update time',linestyle ="dashed", color ="orange", alpha=0.9,linewidth=0.5)

lns = lns1+lns2+lns4

labs = [l.get_label() for l in lns]

for i in range(len(time)):
    my_selected_date = time[i]
    lnv1 = ax1.vlines(my_selected_date, 0, seq_ord[i]+0.03, linestyles ="dashed", colors ="brown", alpha=0.5,linewidth=0.5)   # changed
    #lnv2 = ax1.vlines(my_selected_date, 2.51, ekf[i]+0.005, linestyles ="dashed", colors ="orange",alpha=0.9,linewidth=0.5)   # changed
ax1.axvline(x=time[0],ymin=0.0, ymax=apd[0]+0.03 ,label='localization update time',linestyle ="dashed", color ="brown", alpha=0.2,linewidth=0.5)
#ax1.legend(lns2, labs, loc=1, fontsize=10)
#ax1.legend(loc=4, fontsize=10)
ax1.legend(loc=7, fontsize=10)
ax1.set_title("Outdoor LOS (Total emergency devices localized: 1017)",fontsize=14)
#ax2.set_xlabel("time (sec)",fontsize=14)
plt.ylabel(r"RMSE (m)",fontsize=14)
#ax2.set_ylabel(r"average messages per node (Dr)",fontsize=14)
d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them


# single vline with specific ymin and ymax
#ax2.vlines(x=39.25, ymin=25, ymax=150, colors='green', ls=':', lw=2, label='vline_single - partial height')

# place legend outside
#ax1.legend(loc=1)
ax1.tick_params(axis='both', which='major', labelsize=11)
plt.show()