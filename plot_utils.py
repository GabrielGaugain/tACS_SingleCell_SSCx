#%%
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import patheffects
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np
from listcells import *
from astropy.stats import rayleightest

mm = 1/25.4

def get_continuous_cmap(rgb_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''
    import matplotlib.colors as mcolors
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

rmap = mpl.colormaps["Reds"]
colors = {"L23":"xkcd:carolina blue", "L5":"xkcd:blue", "L6_TPC":"purple", "PC": ["blue","xkcd:carolina blue","purple"],
          "VIP":["orange"]*5 + ["gold"]*5, 
          "SST":"xkcd:dark pastel green",
          "PV":[rmap(0.8)] + 2*[rmap(0.99)] + 3*[rmap(0.6)] + 3*[rmap(0.3)],
          "L4_LBC":"red",
          "L1_NGC":'gray'
          }

def interp(x,y, xi, type_interpolation = "spline"):

    if type_interpolation == "spline":
        from scipy.interpolate import make_interp_spline
        interpolator = make_interp_spline(x, y)
    else:
        from scipy.interpolate import interp1d
        try:
            interpolator  = interp1d(x, y, kind = type_interpolation)
        except:
            print(f"Error: type_interpolation = {type_interpolation} does not exist or is not coded, exiting without interpolation...")
            return

    yi = interpolator(xi)
    return yi

# colors_PV = ["red"] + 2*["xkcd:deep red"] + 3*["xkcd:salmon"] +3*[ "xkcd:tomato"]
# top = mpl.colormaps['Oranges'].resampled(128)
# bottom = mpl.colormaps['Blues_r'].resampled(128)
# newcolors = np.vstack( ( np.vstack((1-bottom(np.linspace(0, 0.5, 64)), gray)),
#                          1-top(np.linspace(0.5, 1, 64)) ) )
# newcmp = ListedColormap(newcolors, name='OrangeGrayBlue')
newcmp = get_continuous_cmap([[1, 1, 1],
                              mpl.colors.to_rgb(colors["L5"])
                              ])

def color_by_cellname(cell_name):
    assert type(cell_name)==str, "cell_name must be a string"

    for cl in colors.keys():
        l = eval(cl)
        if cell_name in l:
            if type(colors[cl])==str:
                return colors[cl]
            else:
                return colors[cl][l.index(cell_name) ]
    return "k"




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