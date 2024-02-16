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


dfs = pd.read_pickle(f"results/allCell_w10Hz.pkl")



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


def get_PLVslopes(freq, df = dfs, metric ="PLV"):

    df_freq = df[df["tacs_freq"]==freq]
    cell_names =  np.unique(df["cell_name"])
    d={}
    dpval = {}
    amps = np.array([0,1,2,3,4,5,10])
    for cell_name in cell_names:
        d[cell_name] = df_freq[df_freq["cell_name"]==cell_name].sort_values("tacs_amp")[metric].values
        dpval[cell_name] = df_freq[df_freq["cell_name"]==cell_name].sort_values("tacs_amp")[metric].values
    d = pd.DataFrame(d, columns=cell_names, index=amps)
    dpval = pd.DataFrame(dpval, columns=cell_names, index=amps)


    ####
    from scipy.stats import linregress
    lrs = {}
    for c in d.columns:
        lr = linregress(d.index.values, d[c].values)
        lrs[c] = [lr.slope, lr.intercept, lr.rvalue, lr.pvalue]
    lrs = pd.DataFrame(lrs, index=["slope", "intercept", "rvalue", "pvalue"])

    return lrs

#%%
freq=5
mean_slopes = {"L5":[], "L23":[], "L6_TPC":[], "SST":[], "VIP":[] , "PV":[]}

for i, freq in enumerate([5,10,20,40]):
    
    lrs = get_PLVslopes(freq=freq)
    # lrs.to_csv(f"results/LinRegPLV_w10Hz_{freq}Hz_tACS.csv")
    # lrs.to_pickle(f"results/LinRegPLV_w10Hz_{freq}Hz_tACS.pkl")
    slopes = lrs.loc["slope"]
    for cl in mean_slopes.keys():
        mean_slopes[cl].append( np.mean(slopes[eval(cl)] ) )

df_slopes = pd.DataFrame(mean_slopes, index=['5', '10', '20', '40'])


#%% Linear regression between polarisation length at one freq and cell size in y dir
###############################################################################
savefig=1
freq = 40
lrs = pd.read_pickle(f"results/LinRegPLV_w10Hz_{freq}Hz_tACS.pkl")


with open("results/cell_bounds.pkl", "rb") as f:
    bounds= pkl.load(f)

with open("results/AllCell_polarisationLength.pkl", "rb") as f:
    Lp = pkl.load(f)
i_freq = np.where(Lp["freqs"]==freq)[0][0]

L, S, y,c = [], [], [], []
# for k in bounds.keys():
for k in PC + SST+ PV:
    if not np.isnan( lrs[k].slope):#  and k not in VIP:
        y.append( np.diff(bounds[k][0:2,1])[0]/1e3 )
        S.append(lrs[k].slope)
        L.append(Lp[k][i_freq])
        c.append(color_by_cellname(k))
L, S, y = np.array(L), np.array(S), np.array(y)
c = np.array(c, dtype=object)


## for correlation between cell size // to EF & polarization Length 
# a,b = y,L 
# xlab, ylab = "y cell length (mm)", f"{freq} Hz Pola. L (mm)"
# figname = f"cellL_ploa/linreg_y_{freq}Hz_Lp"

## for correlation between polarization Length & PLVslopes
a,b = L,S 
xlab, ylab = f"{freq}Hz Pola. Length (mm)", f"{freq}Hz PLV slope (mm/mV)"
# figname = f"PLVxLP/linreg_{freq}Hz_Lp_PLVslope_w10Hz"
figname = f"withoutVIP/linreg_{freq}Hz_Lp_PLVslope_w10Hz"

lr = linregress(a, b)

xn = np.linspace(a.min(), a.max(),100)
yn = xn*lr.slope + lr.intercept

## PLOT
fig, ax = plt.subplots(figsize=(90*mm, 70*mm), dpi=250)
ax.scatter(a, b, color = c, s = 5)
ax.plot( xn, yn, "--", color ="k")
ax.set_xlabel(xlab)
ax.set_ylabel(ylab)
ax = annotate_reg(fig, ax, a, b, lr, s = f"R = {lr.rvalue:.3f} ") #, s = f"y = {round(lr.slope,2)}x + {round(lr.intercept,3)}")

if savefig:
    plt.savefig(f"results/figures/reg/{figname}.jpg",bbox_inches='tight', dpi=250,)
    plt.savefig(f"results/figures/reg/{figname}.svg",bbox_inches='tight', dpi=250)
plt.tight_layout(pad=0)
plt.show()


from scipy.stats import pearsonr, spearmanr
print(pearsonr(a,b ))
print(spearmanr(a,b ))



################################################


#####################################################################
#%% Correaltion for natural activity
# lrs = pd.read_pickle(f"results/LinRegPLV_natural_{10}Hz_tACS.pkl")

df_nat = pd.read_pickle(f"results/AllCell_natural.pkl")

lrs = get_PLVslopes(freq=10, df=df_nat)
# lrs.to_csv(f"results/LinRegPLV_natural_{10}Hz_tACS.csv")
# lrs.to_pickle(f"results/LinRegPLV_natural_{10}Hz_tACS.pkl")

# df_nat["ubPLV"] = df_nat["ubPLV"] .fillna(0)
# lrs = get_PLVslopes(freq=10, df=df_nat, metric="ubPLV")

i_freq = np.where(Lp["freqs"]==10)[0][0]

L, S, c = [], [], []
# for k in bounds.keys():
for k in PC + VIP + SST+ PV:
    if not np.isnan( lrs[k].slope):
        S.append(lrs[k].slope)
        L.append(Lp[k][i_freq])
        c.append(color_by_cellname(k))
L, S, y = np.array(L), np.array(S), np.array(y)
c = np.array(c, dtype=object)

## for correlation between polarization Length & PLVslopes
a,b = L,S 
xlab, ylab = f"{freq} Hz Pola. L (mm)", f"{freq} Hz PLV slope (mm/mV)"
figname = f"linreg_{freq}Hz_polaL_PLVslope_natural"

lr = linregress(a, b)

xn = np.linspace(a.min(), a.max(),100)
yn = xn*lr.slope + lr.intercept

## PLOT
fig, ax = plt.subplots(figsize=(90*mm, 70*mm), dpi=250)
ax.scatter(a, b, color = c, s = 5)
ax.plot( xn, yn, "--", color ="k")
ax.set_xlabel(xlab)
ax.set_ylabel(ylab)
ax = annotate_reg(fig, ax, a, b, lr, s = f"R = {lr.rvalue:.3f} ") #, s = f"y = {round(lr.slope,2)}x + {round(lr.intercept,3)}")

if savefig:
    plt.savefig(f"results/figures/reg/{figname}.jpg",bbox_inches='tight', dpi=250,)
    plt.savefig(f"results/figures/reg/{figname}.jpg",bbox_inches='tight', dpi=250)
plt.tight_layout(pad=0)
plt.show()


print(pearsonr(a,b ))
print(spearmanr(a,b ))





#####################################################################
#%% TEST PLOT WITH HIST
freq = 10
# lrs = pd.read_pickle(f"results/LinRegPLV_w10Hz_{freq}Hz_tACS.pkl")
lrs = pd.read_pickle(f"results/LinRegPLV_natural_{10}Hz_tACS.pkl")

with open("results/cell_bounds.pkl", "rb") as f:
    bounds= pkl.load(f)

with open("results/AllCell_polarisationLength.pkl", "rb") as f:
    Lp = pkl.load(f)
i_freq = np.where(Lp["freqs"]==freq)[0][0]

L, S, y,c = [], [], [], []
for k in bounds.keys():
    if not np.isnan( lrs[k].slope):
        y.append( np.diff(bounds[k][0:2,1])[0]/1e3 )
        S.append(lrs[k].slope)
        L.append(Lp[k][i_freq])
        c.append(color_by_cellname(k))
L, S, y = np.array(L), np.array(S), np.array(y)
c = np.array(c, dtype=object)


a,b = L,S 

lr = linregress(a, b)
xn = np.linspace(a.min(), a.max(),100)
yn = xn*lr.slope + lr.intercept


xlab, ylab = f"{freq} Hz Pola. L (mm)", f"{freq} Hz PLV slope (mm/mV)"
figname = f"linreg_{freq}Hz_polaL_PLVslope_natural"

# Start with a square Figure.
fig = plt.figure(figsize=(6, 6))
# Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
# the size of the marginal axes and the main axes in both directions.
# Also adjust the subplot parameters for a square plot.
gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)
# Create the Axes.
ax = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
# Draw the scatter plot and marginals.
ax.scatter(a, b, color = c, s = 5)
ax.plot( xn, yn, "--", color ="k")
ax.set_xlabel(xlab)
ax.set_ylabel(ylab)
ax = annotate_reg(fig, ax, a, b, lr) #, s = f"y = {round(lr.slope,2)}x + {round(lr.intercept,3)}")

ax_histy.hist(b, orientation='horizontal', bins = 30)
ax_histx.hist(a, bins = 30)

ax_histx.tick_params(axis="x", labelbottom=False)
ax_histy.tick_params(axis="y", labelleft=False)

[ax_histx.spines[side].set_visible(False) for side in ["top", "bottom", "right", "left"]]
[ax_histy.spines[side].set_visible(False) for side in ["top", "bottom", "right", "left"]]

plt.tight_layout(pad=0)
plt.show()
