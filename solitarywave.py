'''Program to (contour) plot solitary wave (u) formed as a result of exitation of FHN equations'''
import numpy as np
import matplotlib.pyplot as pyplot
from myFHN import timeStep
from mpl_toolkits.axes_grid1 import make_axes_locatable


#epsilon, beta, gamma = 0.3, 0.7, 0.5 

dt,totalT = 0.08, 40.0
numbertimeSteps = totalT / dt
timeArray = np.arange(0,totalT, dt)

xMin,xMax,dx = -50.0, 50.0, 0.5
xAxis = np.arange(xMin , xMax , dx) #set x axis 

uData = np.zeros(( int(numbertimeSteps), np.size(xAxis) ), dtype=float) #create 2d array to record variation in u over time
vData = np.zeros(( int(numbertimeSteps), np.size(xAxis) ), dtype=float) #create 2d array to record variation in u over time
u = np.zeros(np.size(xAxis))
v = np.zeros(np.size(xAxis))



if __name__ == "__main__":
    uData, vData = timeStep(u,v)
    
    fig = pyplot.figure()
    ax1 = fig.add_subplot(111)
    
    plot = ax1.imshow(uData[::-1], cmap = 'viridis', aspect = 'auto', extent = (xMin, xMax, 0.0, totalT))
    divider = make_axes_locatable(ax1)     #enable positioning of colourbar
    cax1 = divider.append_axes("right", size="5%", pad=0.15)
    cbar = pyplot.colorbar(plot, orientation = 'vertical',cax=cax1)
    cbar.set_label('u', size=24)
    ax1.set_xlabel('x (m)',size=24)
    ax1.set_ylabel('t (s)', size=24)

    pyplot.show()
    

