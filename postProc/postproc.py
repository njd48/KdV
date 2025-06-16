# -*- coding: utf-8 -*-
"""

"""


import matplotlib.pyplot as plt
from matplotlib import cm
#from matplotlib import colors as mcolors
from matplotlib.collections import LineCollection
import numpy as np

def waterfall_plot(fig,ax,X,Y,Z,**kwargs):
    '''
    Make a waterfall plot
    Input:
        fig,ax : matplotlib figure and axes to populate
        Z : n,m numpy array. Must be a 2d array even if only one line should be plotted
        X,Y : n,m array
        kwargs : kwargs are directly passed to the LineCollection object
    '''
    # Set normalization to the same values for all plots
    norm = plt.Normalize(Z.min().min(), Z.max().max())
    # Check sizes to loop always over the smallest dimension
    n,m = Z.shape
    if n>m:
        X=X.T; Y=Y.T; Z=Z.T
        m,n = n,m

    for j in range(n):
        # reshape the X,Z into pairs 
        points = np.array([X[j,:], Z[j,:]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)  
        # The values used by the colormap are the input to the array parameter
        lc = LineCollection(segments, cmap='plasma', norm=norm, array=(Z[j,1:]+Z[j,:-1])/2, **kwargs)
        line = ax.add_collection3d(lc,zs=(Y[j,1:]+Y[j,:-1])/2, zdir='y') # add line to axes

    fig.colorbar(lc) # add colorbar, as the normalization is the same for all
    # it doesent matter which of the lc objects we use
    ax.auto_scale_xyz(X,Y,Z) # set axis limits

fname = "../simData/u.txt"

data = np.loadtxt( fname, skiprows=8, delimiter="," )

time = data[:,0]

u    = data[:,1:]

S    = u.shape

x    = np.linspace( -np.pi, np.pi, S[1]+1 )
x    = np.delete(x, S[1] )


X,T = np.meshgrid( x, time )

fig  = plt.figure()
ax   = fig.subplots( subplot_kw ={ "projection":"3d"} )
#surf = ax.plot_surface( X,T, u, cmap=cm.Blues,
#                       linewidth=0, antialiased=False)
surf = waterfall_plot( fig, ax, X,T, u )
ax.set_aspect('equalxz')
ax.view_init(elev=20., azim=280)
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('u')
plt.show()