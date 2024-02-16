#%%
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
from plot_utils import mm, color_by_cellname,colors
from BBPresults import load_postres
from listcells import *
from scipy.stats import linregress, spearmanr
## Link : https://raphaelvallat.com/correlation.html

def annotate_reg(fig, ax, x, y, lr, s=""):

    xtext = (x.min() + x.max()) / 2
    ytext = lr.slope*xtext + lr.intercept

    # find run and rise of (xtext, ytext) vs (0, n)
    dx, dy = xtext, ytext - lr.intercept

    # normalize to ax and fig scales
    xfig, yfig = fig.get_size_inches()
    xnorm = dx * xfig / (x.max() - x.min())
    ynorm = dy * yfig / (y.max() - y.min())

    # find normalized annotation angle in radians
    rotation = np.rad2deg(np.arctan2(ynorm, xnorm))
    # rotation = np.rad2deg(np.arctan2(lr.slope, 1))

    if s== "":
        s = f'$R^2={lr.rvalue**2:.3f}$'
    ax.annotate(
        s,
        (xtext, ytext), xycoords='data',
        ha='center', va='bottom',
        rotation=rotation, rotation_mode='anchor',
    )

    return ax


#%% Linear regression between polarisation length at one freq and cell size in y dir
###############################################################################
savefig=1

with open("results/cell_bounds.pkl", "rb") as f:
    bounds= pkl.load(f)

with open("results/AllCell_polarisationLength.pkl", "rb") as f:
    Lp = pkl.load(f)


freq = 0
i_freq = np.where(Lp["freqs"]==freq)[0][0]

L, y,c = [], [], []
for cell in PC + VIP + SST+ PV:
    y.append( np.diff(bounds[cell][0:2,1])[0]/1e3 )
    L.append(Lp[cell][i_freq])
    c.append(color_by_cellname(cell))
L,  y = np.array(L),  np.array(y)
c = np.array(c, dtype=object)

## for correlation between cell size // to EF & polarization Length 
a,b = y,L 
lr = linregress(a, b)
xn = np.linspace(a.min(), a.max(),100)
yn = xn*lr.slope + lr.intercept


## PLOT  FIGURE 8A
fig, ax = plt.subplots(figsize=(90*mm, 70*mm), dpi=250)
ax.scatter(a, b, color = c, s = 5)
ax.plot( xn, yn, "--", color ="k")
ax = annotate_reg(fig, ax, a, b, lr, s =  f"R = {lr.rvalue:.3f} ") #f"R = {spearmanr(a,b ).correlation:.3f}" )#

ax.set_xlabel("y cell length (mm)")
ax.set_ylabel(f"DC Pola. Length (mm)")
# ax.set_ylabel(f"{freq} Hz Pola. L (mm)")
if savefig:
    plt.savefig(f"results/figures/reg/cellL_ploa/linreg_y_{freq}Hz_Lp.jpg",bbox_inches='tight', dpi=250,)
    plt.savefig(f"results/figures/reg/cellL_ploa/linreg_y_{freq}Hz_Lp.svg",bbox_inches='tight', dpi=250)
plt.tight_layout(pad=0)
plt.show()

from scipy.stats import pearsonr, spearmanr
print(pearsonr(a,b ))
print(spearmanr(a,b ))
# %% Correlation over spectrum betwen y and Lp (FIGURE 8B)
################################################
from scipy.stats import pearsonr, spearmanr

with open("results/cell_bounds.pkl", "rb") as f:
    bounds= pkl.load(f)


with open("results/AllCell_polarisationLength.pkl", "rb") as f:
    Lp = pkl.load(f)

freqs = Lp["freqs"]

cl = allCells
dy = np.abs(np.array([ np.diff(bounds[c][0:2,1])[0]/1e3 for c in cl ]))
Lps = np.array([ Lp[c] for c in cl ])
c = np.array([color_by_cellname(c) for c in cl], dtype=object)

corr = "spearman"
corr = "pearson"

# cf = np.array( [linregress(dy, Lps[:,i]).rvalue for i in range(freqs.size) ])
if corr == "spearman":
    cf = np.array([ spearmanr(dy, Lps[:,i]).correlation for i in range(freqs.size) ])
    pval_f = np.array([ spearmanr(dy, Lps[:,i]).pvalue for i in range(freqs.size) ])
elif corr == "pearson":
    cf = np.array([ pearsonr(dy, Lps[:,i]).correlation for i in range(freqs.size) ])
    pval_f = np.array([ pearsonr(dy, Lps[:,i]).pvalue for i in range(freqs.size) ])

fig, ax = plt.subplots(figsize=(90*mm, 70*mm), dpi=250)
ax.plot(freqs,cf, lw=1.2, ls="--", color="black")
ax.set_xlim((-0.5, freqs.max()) )
ax.set_xlabel("frequency (Hz)")
ax.set_ylabel(r"$R_{spearman}$")

plt.savefig(f"results/figures/reg/cellL_ploa/{corr}Corr_spectrum.jpg",bbox_inches='tight', dpi=250,)
plt.savefig(f"results/figures/reg/cellL_ploa/{corr}Corr_spectrum.svg",bbox_inches='tight')
plt.show()

