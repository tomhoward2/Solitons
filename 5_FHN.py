'''Program to model excitation of excitable media according to FHN equations'''
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.animation as animation
import scipy.optimize


epsilon, beta, gamma = 0.3, 0.7, 0.5

dt,totalT = 0.08, 40.0
numbertimeSteps = totalT / dt
timeArray = np.arange(0,totalT, dt)

xMin,xMax,dx = -50.0, 50.0, 0.5
xAxis = np.arange(xMin , xMax , dx) #set x axis 

uData = np.zeros(( int(numbertimeSteps), np.size(xAxis) ), dtype=float) #create 2d array to record variation in u over time
vData = np.zeros(( int(numbertimeSteps), np.size(xAxis) ), dtype=float) #create 2d array to record variation in u over time
u = np.zeros(np.size(xAxis))
v = np.zeros(np.size(xAxis))

zeta = 1.0


def setinitialData():    
    u = 2 *( np.exp(-(xAxis+50)**10/5**10) + np.exp(-(xAxis-50)**10/5**10) ) -1.2
    v[:] = -0.6
             
    return u,v

def f(u,v):
    f = u - (1.0 / 3.0) * u ** 3 - v
    return f

def g(u,v):
    g = u + beta - gamma * v
    return g

def dxxU(u):
    uxx = np.zeros(np.size(u))
    
    for i in range(1,np.size(uxx)-1):       
        uxx[i] = (u[i+1] - 2*u[i] + u[i-1]) / dx ** 2  

    return uxx   

def timeStep(u_i,v_i):
    u_n, v_n = setinitialData()
    
    #no flux boundary conditions
    u_n[0] = u_n[1]
    u_n[-1] = u_n[-2] 
    
    for i in range(int(numbertimeSteps)):
        u_i, v_i =  u_n.copy(), v_n.copy()
        uData[i, :], vData[i, :] = u_n , v_n
        u_n = u_i + (1 / epsilon) * f(u_i,v_i) * dt + dxxU(u_i) * dt
        v_n = v_i + epsilon * g(u_i,v_i) * dt
    
        #hardcode no flux boundary conditions each timetstep
        u_n[0] = u_n[1]
        u_n[-1] = u_n[-2] 

    return uData, vData
        
        
def animate(i, uData, vData):
    x = xAxis
    y1 = uData[i,:]
    y2 = vData[i,:]
    line1.set_data(x,y1)
    line2.set_data(x,y2)
    line3.set_data(uData[i,:], vData[i,:])
    return line1, line2, line3

def init():
    line1.set_data([],[])
    line2.set_data([],[])
    line3.set_data([],[])

    return line1, line2, line3
                
if __name__ == "__main__":       
  uData, vData = timeStep(u,v) 
  
   
  fig = pyplot.figure()
  ax1,ax2, ax3 = fig.add_subplot(131), fig.add_subplot(132), fig.add_subplot(133)     #1-uData, 2-vData, 3-vData.uData
  
  #live plot
  line1, = ax2.plot([], [], lw=1)
  line2, = ax3.plot([], [], lw=1)
  line3, = ax1.plot([], [], lw=2)
  anim = animation.FuncAnimation(fig, animate, init_func=init, frames=int(numbertimeSteps), interval=0,fargs=(uData, vData,), repeat_delay=1000)

  #params
  ax1.set_ylim([-2.0,2.0])
  ax2.set_ylim([-5.0, 5.0])
  ax3.set_ylim([-5.0, 5.0])
  ax2.set_xlim([-50.0,50.0])
  ax3.set_xlim([-50.0,50.0])
  ax1.set_ylabel('v')
  ax1.set_xlabel('u')
  ax2.set_ylabel('u')
  ax2.set_xlabel('x')
  ax3.set_ylabel('v')
  ax3.set_xlabel('x')
  u_nullcline, v_nullcline = np.arange(-2.2,2.2,0.01), np.arange(-2.2,2.2,0.01)
  ax1.plot(u_nullcline, u_nullcline-(1./3.)*u_nullcline**3) 
  ax1.plot(v_nullcline,(1./gamma)*(u_nullcline + beta))

  pyplot.show()
    
