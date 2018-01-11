'''Program to propagate sech soliton confining and non-confining potential'''

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable


#Global Variables
totalT = 5.0
numberIterations = 800.0
deltaT = totalT / numberIterations 
timeArray = np.linspace(0.0, totalT, 800.0)
myLambda = 0.05 
xMin, xMax, dx = -20, 20, 0.05
xAxis = np.arange(xMin , xMax , dx) #set x axis 
zeta, v = 1.0, ((xMax - xMin) / 4)
dk = 2 * np.pi / (xMax - xMin)
kDiff = np.size(xAxis) * dk
k = np.arange(-kDiff / 2, kDiff / 2, dk)
#k = 2 * np.pi * np.fft.fftfreq((np.size(xAxis)), myLambda) #equaivalent to above
nHarmonic = np.zeros(np.size(xAxis))
g=-1.0

    
    
def waveFunction(xAxis):
    '''Function to return wavefunction'''
    psi = np.sqrt(2.0 * zeta) * np.exp(1j * v * xAxis) / np.cosh(xAxis * np.sqrt(zeta))
    return psi
    
def getPotentialE(psi):
    potentialE = g * np.real(np.real(psi) ** 2 + np.imag(psi) ** 2) #remove 0j component
    return potentialE
    
def wavefunctionDeltaT(psi, V, xAxis):
    '''Function to return wavefunction at time deltaT later using split step fourier'''     
    step1 = np.exp(-1j * V * deltaT / 2) * psi
    step2 = np.fft.fftshift(step1)
    step3 = np.fft.fft(step2)  
    step4 = np.fft.fftshift(step3) * dx
    
    step5 = np.exp(-1j * k ** 2 * deltaT) * step4
    step6 = np.fft.fftshift(step5)
    step7 = np.fft.ifft(step6)
    step8 = np.fft.fftshift(step7) * 1.0 / dx
    
    step9 = np.exp(-1j * V * deltaT / 2) * step8
    
    
    return step9
    

def animate(i, wavefunctionProbabilityConserved, wavefunctionProbabilitySpread):
    x = timeArray
    y1 = wavefunctionProbabilityConserved[int(numberIterations - 1) - i,:]
    y2 = wavefunctionProbabilitySpread[int(numberIterations - 1) - i,:]
    line1.set_data(x,y1)
    line2.set_data(x,y2)
    return line1,line2

def init():
        line1.set_data([],[])
        line2.set_data([],[])
        return line1, line2
        
if __name__ == "__main__":   
    
    #set up plots
    fig = pyplot.figure()
    ax1, ax2 , ax3 = fig.add_subplot(221), fig.add_subplot(222) , fig.add_subplot(212)
    ax1.set_xlabel('x (m)', size=24)
    ax1.set_ylabel('t (s)', size=24)
    ax1.set_title('Propagation sech soliton with confining potential',size=26)
    ax2.set_xlabel('x (m)', size=24)
    ax2.set_ylabel('t (s)', size=24)
    ax3.xaxis.set_ticklabels([])
    #ax2.set_ylabel(r'N = $\int_a^b |\psi(x,t)|^2 dx$')
    ax2.set_title('Propagation sech soliton with non-negative potential',size=26)
    ax3.set_xlabel('x',size=24)
    ax3.set_ylabel('$|\psi(x,t)|^2$', size=24)
    myAxes=(ax1,ax2,ax3)
   
    
    for ax in myAxes:
        psi = waveFunction(xAxis)        
        wavefunctionProbability = np.zeros(( int(numberIterations), np.size(xAxis) ), dtype=float) #create 2d array to record variation in psi over time
    
        #set form of potential
        if ax==ax1:
            g = -1.0
        elif ax==ax2:
            g = 1.0
        
        #timestep
        for i in range(int(numberIterations)):
            potentialE = getPotentialE(psi)
            psiSquared = (np.real(psi) ** 2 + np.imag(psi) ** 2) #=conj(psi) * psi.Factor of 0.5
            wavefunctionProbability[int(numberIterations - 1) - i , :] = np.real(psiSquared)     #populate array from bottom up. This is is real anyway; just exclude 0*j
            
            #check N conservation
            nHarmonic[i] = np.sum(np.real(psiSquared))                     
            
            psi = wavefunctionDeltaT(psi, potentialE, xAxis)      #find new psi dt later                
            
        
        if (ax==ax1 or ax==ax2):
            #plot wavefunction        
            #seems to be a normalisation factor of 1/2 missing from sech soliton plot
            plot = ax.imshow(wavefunctionProbability, cmap = 'viridis', aspect = 'auto', extent = (xMin, xMax, 0.0, totalT))
            divider = make_axes_locatable(ax)     #enable positioning of colourbar
            cax1 = divider.append_axes("right", size="5%", pad=0.15)
            cbar = pyplot.colorbar(plot, orientation = 'vertical',cax=cax1)
            cbar.set_label(r'$|\psi(x,t)|^2$', size=24)
                        
            if (ax==ax1):
                wavefunctionProbabilityConserved = wavefunctionProbability
                print "Psi squared conserved to within " + str(100.0 * ((nHarmonic[0] - nHarmonic[-1]) / nHarmonic[0])) + " %"
            else:
                wavefunctionProbabilitySpread = wavefunctionProbability
                
                        
        #plot N conservation
        #ax2.plot(timeArray, nHarmonic) 
        #print "Psi squared conserved to within " + str(100.0 * ((nHarmonic[0] - nHarmonic[-1]) / nHarmonic[0])) + " %"
        
        
        if (ax==ax3):
            #live plot
            ax3.set_xlim([0.0, totalT])
            ax3.set_ylim([0.0, 2.2])
            line1, = ax3.plot([], [], lw=1)
            line2, = ax3.plot([], [], lw=1)
        
            #live plot
            anim = animation.FuncAnimation(fig, animate, init_func=init, frames=int(numberIterations), interval=0,fargs=(wavefunctionProbabilityConserved, wavefunctionProbabilitySpread,), repeat_delay=2000)
        

    
    
    
    pyplot.show()
    
    
   
    
    
    
    
     




