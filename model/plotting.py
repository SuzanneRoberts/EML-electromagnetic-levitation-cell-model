# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 09:10:53 2015

@author: Suzanne
"""
import numpy as np, pylab as pl
# For: Hack found online to remove the perspective view from the plot
from mpl_toolkits.mplot3d import proj3d
# For plotting
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.collections import PolyCollection
import matplotlib as mpl
import itertools
import pickle


# Hack found online to remove the perspective view from the plot
def orthogonal_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    return np.array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,a,b],
                        [0,0,-0.0001,zback]])
                            
                
# Function to plot computed fields without interpolation
def patchSampleField(patchVerts, fieldvec, titleString, colorbarString):
    fig, ax = plt.subplots()

    # Make the collection and add it to the plot.
    coll = PolyCollection(patchVerts*1000.0, array=fieldvec, cmap=mpl.cm.jet, edgecolors='none')
    #coll = PolyCollection(patchVerts*1000.0, array=fieldvec, cmap=mpl.cm.jet, edgecolors='none',norm=mpl.colors.LogNorm())
    ax.add_collection(coll)
    ax.autoscale_view()
    #coll.set_clim(0, 5E6)

    # Add a colorbar for the PolyCollection
    #fig.colorbar(coll, ax=ax)
    #plt.show()

    # Labels and colorbar settings
    cb = fig.colorbar(coll, shrink=0.5, aspect=5)
    cb.set_label(colorbarString)
    ax.set_title(titleString)
    ax.set_xlabel('mm')
    ax.set_ylabel('mm')
    
    pl.axis('equal')
    #plt.show()

                            
# Function to plot a field over the sample
def plotSampleField(xvec, yvec, fieldvec, titleString, colorbarString):
    
    # Hack found online to remove the perspective view from the plot
    proj3d.persp_transformation = orthogonal_proj
    
    # Plotting the field
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf=ax.plot_trisurf(xvec*1000, yvec*1000, fieldvec, cmap=cm.jet, linewidth=0.0)
    
    ax.view_init(90,-90)
    ax.set_zticks([])
    cb = fig.colorbar(surf, shrink=5, aspect=5)
    cb.set_label(colorbarString)
    ax.set_title(titleString)
    ax.set_xlabel('mm')
    ax.set_ylabel('mm')


# Function to plot a field over the sample - mirrored along the symmetry axis
def plotSampleFieldFull(xvec, yvec, fieldvec, titleString, colorbarString):
    
    # Hack found online to remove the perspective view from the plot
    proj3d.persp_transformation = orthogonal_proj
    
    # Plotting the field
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    fullX = np.concatenate((xvec,-xvec))
    fullY = np.concatenate((yvec,yvec))
    fullZ = np.concatenate((fieldvec,fieldvec))
    surf=ax.plot_trisurf(fullX*1000, fullY*1000, fullZ, cmap=cm.jet, linewidth=0.0)
    
    ax.view_init(90,-90)
    ax.set_zticks([])
    cb = fig.colorbar(surf, shrink=0.5, aspect=5)
    cb.set_label(colorbarString)
    ax.set_title(titleString)
    ax.set_xlabel('mm')
    ax.set_ylabel('mm')
    

# Function to plot nodes
def plotNodes(x,y):
    pl.figure()
    pl.plot(x,y,'-bo')
    pl.axis('equal')
    

# Function to format plots for a dissertation
def plotDissertation(plotTuple, titleString, xString='', yString='', legendList=[], locString='upper right'):
    # add later: error catch if len(plotTuple) is not divisible by two
    marker = itertools.cycle(('-ko','-k^','-k*','-kv','-ks','-k+'))
    colourMarker = itertools.cycle(('-bo','-r^','-g*','-mv','-cs','-k+'))

    if pl.get_fignums() == []:
        newFigNum = 1
    else:
        newFigNum = np.max(pl.get_fignums()) + 1
        
    for i in np.arange(len(plotTuple)/2):
        pl.figure(newFigNum)
        pl.plot(np.array(plotTuple[int(2*i)]), np.array(plotTuple[int((2*i)+1)]), next(marker))
        pl.figure(newFigNum+1)
        pl.plot(np.array(plotTuple[int(2*i)]), np.array(plotTuple[int((2*i)+1)]), next(colourMarker))
    
    for n in newFigNum, newFigNum+1:
        pl.figure(n)
        pl.xlabel(xString)
        pl.ylabel(yString)
        pl.legend(legendList, loc=locString)
        pl.grid(True)  
        pl.savefig("outputFigures/'+titleString+str(n)+'.pdf", format="pdf", bbox_inches="tight")
    
    # np.savetxt('outputFigures/'+titleString+'.txt', plotTuple) # only works for 2D arrays... therefore rather pickle
    with open('outputFigures/'+titleString+'.pickle', 'wb') as f:
        pickle.dump(plotTuple, f)
    
# Polar contour plot - copied from:
def plot_polar_contour(values, azimuths, zeniths):
    """Plot a polar contour plot, with 0 degrees at the North.
 
    Arguments:
 
     * `values` -- A list (or other iterable - eg. a NumPy array) of the values to plot on the
     contour plot (the `z` values)
     * `azimuths` -- A list of azimuths (in degrees)
     * `zeniths` -- A list of zeniths (that is, radii)
 
    The shapes of these lists are important, and are designed for a particular
    use case (but should be more generally useful). The values list should be `len(azimuths) * len(zeniths)`
    long with data for the first azimuth for all the zeniths, then the second azimuth for all the zeniths etc.
 
    This is designed to work nicely with data that is produced using a loop as follows:
 
    values = []
    for azimuth in azimuths:
      for zenith in zeniths:
        # Do something and get a result
        values.append(result)
 
    After that code the azimuths, zeniths and values lists will be ready to be passed into this function.
 
    """
    #theta = np.radians(azimuths)
    azimuths = np.array(azimuths)
    zeniths = np.array(zeniths)
 
    values = np.array(values)
    values = values.reshape(len(zeniths), len(azimuths))
 
    #r, theta = np.meshgrid(zeniths, np.radians(azimuths))
    r, theta = np.meshgrid(zeniths, azimuths)
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    #plt.autumn()
    cax = ax.contourf(theta, r, np.transpose(values), 50, norm=mpl.colors.LogNorm())
    #cax = ax.contour(theta, r, np.transpose(values), 20)
    #plt.autumn()
    cb = fig.colorbar(cax)
    cb.set_label("Pixel reflectance")
    
    #fig = plt.figure()
    return fig, ax, cax
    
    
testContourPlot = 0
if testContourPlot == 1:    
    a = np.array([0, np.pi/6, 2*np.pi/6, np.pi/2, 4*np.pi/6, 5*np.pi/6, np.pi])
    b = np.array([1, 2])
    # = np.array([5, 11, 5, 11, 5, 11, 5, 11, 5, 11, 5, 11, 5, 11])
    c = np.array([5, 5, 5, 5, 5, 5, 5, 11, 11, 11, 11, 11, 11, 11])
           
    plot_polar_contour(c, a, b)
    
testPlotDissertation = 0
if testPlotDissertation == 1:
    a = np.linspace(0,9,10)
    b = np.sin(a)
    c = a/2
    d = np.cos(c)
    plotDissertation((a,b,c,d), 'testPlotDissertation', 'xlabel', 'ylabel', ['sin(x)','cos(x/2)'])
    #plotDissertation(np.loadtxt('outputFigures/testPlotDissertation.txt'), 'testPlotDissertation2', 'xlabel', 'ylabel', ['sin(x)','cos(x/2)'])
    with open('outputFigures/testPlotDissertation.pickle', 'rb') as f:
        testPlotDissertationLoadedData = pickle.load(f)
    plotDissertation(testPlotDissertationLoadedData, 'testPlotDissertationLoadedData', 'xlabel', 'ylabel', ['sin(x)','cos(x/2)'])
