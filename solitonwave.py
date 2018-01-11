'''Program to plot sech-soliton collision from milestone'''

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable

from milestoneSoliton import wavefunctionDeltaT

xMin, xMax, dx = -20, 20, 0.05
xAxis = np.arange(xMin , xMax , dx) #set x axis 
numberIterations = 800.0
nHarmonic = np.zeros(np.size(xAxis))
g = -4.0
zeta , v = 1.0, ((xMax - xMin) / 4)
totalT = 5.0


def waveFunction(xAxis):
    '''Function to return wavefunction'''
    psi = np.sqrt(2.0 * zeta) * (np.exp(1j * v * xAxis) / np.cosh((xAxis+20) * np.sqrt(zeta)) + np.exp(-1j * v * xAxis) / np.cosh((xAxis-20) * np.sqrt(zeta)))
    return psi

def getPotentialE(psi):
    potentialE = g * np.real(np.real(psi) ** 2 + np.imag(psi) ** 2) #remove 0j component
    return potentialE

def timestep(psi):
    for i in range(int(numberIterations)):
        potentialE = getPotentialE(psi)
        psiSquared = (np.real(psi) ** 2 + np.imag(psi) ** 2) #=conj(psi) * psi.Factor of 0.5
        wavefunctionProbability[int(numberIterations - 1) - i , :] = np.real(psiSquared)     #populate array from bottom up. This is is real anyway; just exclude 0*j
        
        #check N conservation
        nHarmonic[i] = np.sum(np.real(psiSquared))                     
        
        psi = wavefunctionDeltaT(psi, potentialE, xAxis)
    
    return wavefunctionProbability
    

if __name__ == "__main__":
    fig = pyplot.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('x (m)', size=24)
    ax1.set_ylabel('t (s)', size=24)
    
    psi = waveFunction(xAxis)        
    wavefunctionProbability = np.zeros(( int(numberIterations), np.size(xAxis) ), dtype=float) #create 2d array to record variation in psi over time

    wavefunctionProbability = timestep(psi)
    print "Psi squared conserved to within " + str(100.0 * ((nHarmonic[0] - nHarmonic[-1]) / nHarmonic[0])) + " %"
    
    
    plot = ax1.imshow(wavefunctionProbability, cmap = 'viridis', aspect = 'auto', extent = (xMin, xMax, 0.0, totalT))
    divider = make_axes_locatable(ax1)     #enable positioning of colourbar
    cax1 = divider.append_axes("right", size="5%", pad=0.15)
    cbar = pyplot.colorbar(plot, orientation = 'vertical',cax=cax1)
    cbar.set_label(r'$|\psi(x,t)|^2$', size=24)
    ax1.set_ylim([0,2])
    pyplot.show()

    
    
