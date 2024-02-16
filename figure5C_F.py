#%%
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, time,efel
import pickle as pkl
from astropy.stats import rayleightest
import pandas as pd
from listcells import *

plt.rcParams.update({        
    "xtick.direction": 'in',
    "ytick.direction": 'in',
})
dt = 0.025
# tstart = time.time()
amps = np.array([0,1,2,3,4,5,10])

def polar_twin(ax):
    ax2 = ax.figure.add_axes(ax.get_position(), projection='polar', 
                             label='twin', frameon=False,
                             theta_direction=ax.get_theta_direction(),
                             theta_offset=ax.get_theta_offset())
    ax2.xaxis.set_visible(False)
    # There should be a method for this, but there isn't... Pull request?
    ax2._r_label_position._t = (22.5 + 180, 0.0)
    ax2._r_label_position.invalidate()
    # Ensure that original axes tick labels are on top of plots in twinned axes
    # for label in ax.get_yticklabels():
    #     ax.figure.texts.append(label)
    from matplotlib.ticker import MaxNLocator
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=len(ax.get_yticks())))

    #Set the last tick as the plot limit
    ax2.set_ylim(0, ax2.get_yticks()[-1])

    #Remove the tick label at zero
    ax2.yaxis.get_major_ticks()[0].label1.set_visible(False)
    return ax2


## Phases Histogram
def PLV_polar_plot(ax, phases, PLVs, cmap = mpl.colormaps["Blues"], 
                   rticks_loc= 330,rticks =None, hide_ticks=False, marker="o",amps=amps,  legend=True):

    # values = values/np.sum(values)
    # fig, ax = plt.subplots(figsize = (10,12),subplot_kw={'projection': 'polar'})
    scat = ax.scatter(phases,PLVs, marker= marker, #markers[i%5],  
                c=amps,cmap = cmap, edgecolors="k", linewidths=0.2)
    if legend:
        ax.legend(scat.legend_elements()[0], [f"{amp} V/m" for amp in amps],
                    title=" ", bbox_to_anchor=(1.35, 0.55),loc='center', frameon=False,borderaxespad=0.)
        for ha in ax.legend_.legendHandles:
            ha.set_markeredgecolor("k") 
            ha.set_markeredgewidth(0.8) 
            ha.set_markersize(10)
    ## settings
    ax.set_rlabel_position(rticks_loc)  # Move radial labels away from plotted line
    # ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2], ["0", r"$\pi/2$",r"$\pi$",r"$3\pi/2$"])
    ax.spines["polar"].set_color('lightgray') ## contour in gray
    ax.tick_params(axis='both', which='major')
    ax.grid(linestyle=":",alpha=0.8)
    # ax.set_axisbelow(True)  ## grid behing bars
    if rticks is not None:
        ax.set_rticks(rticks)
    if hide_ticks:
        ax.set_xticklabels([])
        ax.set_yticklabels([])    
    return ax, scat

def hist_plot(ax, phi, dtheta=10, meanVec = False, p_val = False, 
              ticks = None, polar=False, hide_ticks = False, rticks_loc= 330,
              kwargs = {"color":"skyblue", "linewidth":0.5,  "edgecolor":"w", }, pvalsize=20):

    ## creating histograms from data (angles)
    theta = np.arange(0 ,360 + dtheta,dtheta)*np.pi/180
    values, th = np.histogram(phi, bins =theta, density=True)
    th = (th[1:] + th[:-1])/2
    values /=  values.max()
    # fig, ax = plt.subplots(figsize = (10,12),subplot_kw={'projection': 'polar'})
    ax.bar(th, values,   width = dtheta/180*np.pi,rasterized=True,  **kwargs)  #,alpha = 0.8)

    # mini, maxi = np.min(values), np.max(values)
    # avg, std = np.mean(values), np.std(values)
    # ax.plot(theta,(mini+maxi)/2 + (maxi-mini)/2*np.sin(theta), color ="orange")
    # ax.plot(np.arange(0,2*np.pi+0.1, 0.1),avg + 1.96*std*np.sin(np.arange(0,2*np.pi+0.1, 0.1)), color ="orange")
    if meanVec:
        # if polar:
        #     ax2 = polar_twin(ax)
        #     ax2.set_rlabel_position(200.5)
        #     plt.setp(ax2.get_yticklabels(), color='red')
        # else:
        #     ax2 = ax.twinx()
        #     ax2.set_ylabel("PLV")
        # # ax.bar(pref_phase, 30,width=0.03, color="red",align='center')
        ax2 = ax
        mVec = np.mean(np.exp(1j*phi))
        (x1,x2) = ([np.angle(mVec), np.angle(mVec)], [0,np.abs(mVec)]) if not polar else ([np.angle(mVec), np.angle(mVec)], [0.005, np.abs(mVec)]) 
        # ax2.arrow(x1[0], x2[0], x1[1]-x1[0], x2[1], color="r", width=0.05)
        import matplotlib.patheffects as pe
        ax2.plot(x1, x2, color="k", linewidth=5, solid_capstyle='round')
        ax2.plot(x1, x2, color="red", linewidth=4, solid_capstyle='round')

        # ax2.plot(0, 0, marker="o",color = "red", markersize=4)
        # ax2.set_ylim((0,1)) #if not polar else ax2.rlim((0,1))

    if p_val:
        pval = rayleightest(phi)
        txt_stat = "***" if pval<0.001 else "**" if pval<0.01 else "*" if pval<0.05 else ""
        pos = (50/180*np.pi, 1.3*values.max()) if polar else (270/180*np.pi, 0.95*values.max())
        ax.text(pos[0], pos[1], txt_stat, 
            fontsize=pvalsize,  ha='center' )


    ## settings
    ax.set_rlabel_position(rticks_loc)  # Move radial labels away from plotted line
    # ax.spines["polar"].set_color('lightgray') ## contour in gray

    ax.set_xlim((0,2*np.pi))
    if not polar:
        ax.set_xlabel(r"$\phi$ (rad)")    
        ax.set_ylabel(r"$Pr. spike$ ")    
        ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi], ["0", r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$"])
    else:
        ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2], ["0", r"$\pi/2$",r"$\pi$",r"$3\pi/2$"])

    ax.set_axisbelow(True)  ## grid behing bars
    ax.grid(linestyle=":",alpha=0.8)

    if hide_ticks:
        ax.set_xticklabels([])
        ax.set_yticklabels([]) 

    return ax, (th, values)
## END hist_plot()



df = pd.read_pickle(f"results/allCell_w10Hz.pkl")

#%% model simplifiés
PLV_L5simp= np.array([ 0.02498163, 0.03809009, 0.05759349, 0.0829909, 
              0.10911719, 0.14227304, 0.16693263, 0.18051044, 0.2042792,  0.22796528, 0.24679369])
Angle_L5simp = np.array([-40.9749756 ,   5.25746317,  26.61683929,  32.47714836,
        38.99471026,  42.53389491,  40.04830171,  41.01991015,
        41.89432673,  42.91074682,  42.22147784])*np.pi/180

PLV_PVsimp = np.array([0.044425,   0.03841095, 0.04885117, 0.06458866, 0.09071947, 0.11788052,
               0.15013993, 0.17664711, 0.21321898, 0.229326,   0.25997988])
Angle_PVsimp = np.array([-234.99449811, -268.04249624,   43.77545736,   37.20041833,
         34.68327354,   32.56862158,   29.20976699,   26.01770186,
         21.02716242,   16.78992187,   15.51133091])*np.pi/180


#%% PLot of one cell mean VECTOR in polar plot for diff amp (scatter as Figure 5C)
from plot_utils import colors, get_continuous_cmap, mm,color_by_cellname
savefig = 1
suf = "w10Hz"
amps = np.array([0,1,2,3,4,5,10])

freq =10
cell_name = L5[0]

col = mpl.colors.to_rgb(color_by_cellname(cell_name))
cmap = get_continuous_cmap([[1,1,1],col] )
cmap = get_continuous_cmap([cmap(0.1),col] )

PLVs = df[(df["cell_name"]==cell_name)&(df["tacs_freq"]==freq)]["PLV"].values 
phases = df[(df["cell_name"]==cell_name)&(df["tacs_freq"]==freq)]["PPh"].values 

fig, ax = plt.subplots(figsize = (85*mm,60*mm) ,subplot_kw={'projection': 'polar'})
ax,sc = PLV_polar_plot(ax, phases, PLVs, cmap=cmap,  rticks_loc=20, legend=False)

ax,_ = PLV_polar_plot(ax, Angle_L5simp, PLV_L5simp, cmap=cmap,  rticks_loc=20,
                       marker="^",amps=np.arange(0,11,1), legend=False)

fig.colorbar(sc, fraction= 0.06, pad=0.12, aspect= 15,ticks=[0,10],label = "EF intensity (V/m)")
ax.set_rticks([0.1,0.2,0.3])
# ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2], ["0", r"$\pi/2$",r"$\pi$",r"$3\pi/2$"])

## saving
if savefig:
    figdir = f"figures/figure5"
    plt.savefig(f"{figdir}/5B.svg",bbox_inches='tight', dpi=250)
    plt.savefig(f"{figdir}/5B.jpg", bbox_inches='tight')
    # plt.close()

plt.show()


#%% FIGURE 5F

freq =10
cell_name = L5[0]
# cell_name = PV[4]

col = mpl.colors.to_rgb(color_by_cellname(cell_name))
cmap = get_continuous_cmap([[1,1,1],col] )
cmap = get_continuous_cmap([cmap(0.1),col] )


PLVs = df[(df["cell_name"]==cell_name)&(df["tacs_freq"]==freq)]["PLV"].values 
phases = df[(df["cell_name"]==cell_name)&(df["tacs_freq"]==freq)]["PPh"].values 

fig, ax = plt.subplots(figsize = (85*mm,60*mm) ,subplot_kw={'projection': 'polar'})
ax,sc = PLV_polar_plot(ax, phases, PLVs, cmap=cmap,  rticks_loc=330, legend=False)

ax,_ = PLV_polar_plot(ax, Angle_PVsimp, PLV_PVsimp, cmap=cmap,  rticks_loc=330,
                     marker="^",amps=np.arange(0,11,1), legend=False)

fig.colorbar(sc, fraction= 0.06, pad=0.12, aspect= 15,ticks=[0,10],label = "EF intensity (V/m)")
ax.set_rticks([0.1,0.2,0.3])
# ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2], ["0", r"$\pi/2$",r"$\pi$",r"$3\pi/2$"])

## saving
if savefig:
    figdir = f"figures/figure5"
    plt.savefig(f"{figdir}/5F.svg",bbox_inches='tight', dpi=250)
    plt.savefig(f"{figdir}/5F.jpg", bbox_inches='tight')
    # plt.close()

plt.show()