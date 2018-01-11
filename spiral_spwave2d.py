'''Program to model spiral waves using FHN equations - spiral perturbations contour plot'''
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable



epsilon, beta, gamma = 0.3, 0.7, 0.5 

dt,totalT = 0.05, 10.0
numbertimeSteps = totalT / dt
timeArray = np.arange(0,totalT, dt)

xMin,xMax,dx = -30.0, 30.0, 0.5
xAxis = np.arange(xMin , xMax , dx) #set x axis 
yMin,yMax,dy = -30.0, 30.0, 0.5
yAxis = np.arange(yMin , yMax , dy) #set y axis

u = np.zeros((np.size(yAxis),np.size(xAxis)), dtype=float)
v = np.zeros((np.size(yAxis),np.size(xAxis)), dtype=float)
uMin, uMax = np.zeros(int(numbertimeSteps)), np.zeros(int(numbertimeSteps))


#uData = np.zeros(( int(numbertimeSteps),  ))  #create 2d array to record variation in u over time
#vData = np.zeros(( int(numbertimeSteps), )) #create 2d array to record variation in u over time

uData, vData = [], []

zeta = 1.0
c = 0.2

def setinitialData():       
    x, y = np.meshgrid(xAxis,yAxis)
    u = c * x
    v = c * y 
    
    return u,v

def f(u,v):
    f = u - (1.0 / 3.0) * u ** 3 - v
    return f

def g(u,v):
    g = u + beta - gamma * v
    return g

def grad2U(u):    
    grad2u = np.zeros( (np.size(u[:,0]), np.size(u[0,:])) )

    for j in range(1,np.size(grad2u[:,0])-1): 
        for i in range(1,np.size(grad2u[0,:])-1):       
            grad2u[j,i] = (u[j,i+1] + u[j,i-1] + u[j+1,i] + u[j-1,i] - 4*u[j,i]) / dx ** 2   #dx and dy are same size

    return grad2u

def timeStep(u_i,v_i):    
    u_n, v_n = setinitialData()

    #no flux boundary conditions
    u_n[0,:], u_n[-1,:] = u_n[1,:], u_n[-2,:]
    u_n[:,0], u_n[:,-1] = u_n[:,1], u_n[:,-2]
    v_n[0,:], v_n[-1,:] = v_n[1,:], v_n[-2,:]
    v_n[:,0], v_n[:,-1] = v_n[:,1], v_n[:,-2]
    
    
    for i in range(int(numbertimeSteps)):
        u_i, v_i =  u_n.copy(), v_n.copy()
        #uData[i, :], vData[i, :] = u_n , v_n
        uData.append(u_n)
        vData.append(v_n)
        u_n = u_i + (1 / epsilon) * f(u_i,v_i) * dt + grad2U(u_i) * dt
        v_n = v_i + epsilon * g(u_i,v_i) * dt
    
        #hardcode no flux boundary conditions each timetstep
        u_n[0,:], u_n[-1,:] = u_n[1,:], u_n[-2,:]
        u_n[:,0], u_n[:,-1] = u_n[:,1], u_n[:,-2]
        v_n[0,:], v_n[-1,:] = v_n[1,:], v_n[-2,:]
        v_n[:,0], v_n[:,-1] = v_n[:,1], v_n[:,-2] 
        
        uMax[i] = np.max(u_n)
        

    return uData, vData
        
        

def animate(i, uData):
    ax1.clear()
    x = xAxis
    y = yAxis
    X, Y = np.meshgrid(x,y)
    #cont = pyplot.pcolormesh(X,Y,uData[i], cmap=cm.jet)
    cont = ax1.imshow(uData[i+1], cmap = 'viridis')
    ax1.set_xlabel('x', size=24)
    ax1.set_ylabel('y', size=24)
    ax1.legend()
    
    
    return cont


                
if __name__ == "__main__":       
  uData, vData = timeStep(u,v) 
     
  fig = pyplot.figure()
  ax1 = fig.add_subplot(111)
  
  cont = ax1.imshow(uData[0], cmap='viridis')
  
  anim = animation.FuncAnimation(fig, animate, frames=int(numbertimeSteps-1), interval=0,fargs=(uData,), repeat_delay=1000)

  #Z = uData[np.argmax(uMax)] 
  #m = cm.ScalarMappable(cmap=cm.jet)
  #m.set_array(Z)
  #fig.colorbar(m)
  
  divider = make_axes_locatable(ax1)     #enable positioning of colourbar
  cax1 = divider.append_axes("right", size="5%", pad=0.15)
  cbar = pyplot.colorbar(cont, orientation = 'vertical',cax=cax1)
  cbar.set_label('u', size=26)
  
  


  pyplot.show()
    
