#d2x/dt2 + (w^2)*x=0 given initial condition x(t=0)=0 and dx/dt(t=0)=0
#dz/dt=-(w^2)*x, dx/dt=z
#z(t0+h)=z(t0)-h*(w^2)*x(t0), x(t0+h)=x(t0)+h*z(t0) and so on
#so the itterative formula for euler's method is
#                                               z(t(n+1))=z(t(n))-h*(w^2)*(x(tn)) and x(t(n+1))=x(tn)+h*z(tn)

import numpy as np
import matplotlib.pyplot as plt

t0=0#<-initial time
tn=4*(np.pi)#<-final time
n=1000#<-total no. of points-1
h=(tn-t0)/n
time=[]
for i in range(n+1):
    t=t0+i*h
    time.append(t)

#print(time)
#initial conditions
x0=1#<-initial position
z0=0#<-initial velocity
w=1#<-angular vel

def sinw(w,xx):
    sinw=[]
    for i in range(len(xx)):
        s=np.sin(w*xx[i])
        sinw.append(s)
    return(sinw)

def eulersoln(xx,w,x0,z0):#x0=initial position z0=initial velocity,w=angular vel
    eulersoln=[]
    eulersolnvel=[]
    xold=x0
    velold=z0
    for i in range(len(xx)):
        xnew=xold+h*velold
        velnew=velold-h*(w**2)*xold
        eulersoln.append(xnew)
        xold=xnew
        velold=velnew

    return(eulersoln)

def eulervelsoln(xx,w,x0,z0):#x0=initial position z0=initial velocity,w=anglular vel
    eulersoln=[]
    eulersolnvel=[]
    xold=x0
    velold=z0
    for i in range(len(xx)):
        xnew=xold+h*velold
        velnew=velold-h*(w**2)*xold
        eulersoln.append(xnew)
        eulersolnvel.append(velnew)
        xold=xnew
        velold=velnew

    return(eulersolnvel)

def rk2soln(eulersoln,eulervelsoln,w,x0,z0):        # eulers formula z(t(n+1))=z(t(n))-h*(w^2)*(x(tn)) and x(t(n+1))=x(tn)+h*z(tn)
    
    xrksoln=[]
    xrkold=x0
    velrkold=z0

    for i in range(len(eulersoln)):
        velrknew=velrkold - (h/2)*(((w**2)*xrkold)+ (w**2)*eulersoln[i])
        xrknew=xrkold +(h/2)*(velrkold+eulervelsoln[i])
        xrksoln.append(xrknew)
        xrkold=xrknew
        velrkold=velrknew

    return(xrksoln)

#print(rk2soln(eulersoln(time,1,0,1),eulervelsoln(time,1,0,1),1,0,1))
#print(eulersoln(time,1,0,1))

plt.xlabel('time')
plt.ylabel('displacement')
plt.plot(time,np.cos(time),color='g',label='actual soln')
plt.plot(time,eulersoln(time,w,x0,z0),color='r',label='eulers soln')
plt.plot(time,rk2soln(eulersoln(time,w,x0,z0),eulervelsoln(time,w,x0,z0),w,x0,z0),color='b',label='rk2 soln')
plt.title('d2x/dt2 + x = 0,initial position=1,initial vel=0 h=0.01256')
#plt.plot(time,(eulersoln(time,w,x0,z0)-np.cos(time)),color='r',label='euler')
#plt.plot(time,(rk2soln(eulersoln(time,w,x0,z0),eulervelsoln(time,w,x0,z0),w,x0,z0)-np.cos(time)),color='b',label='rk2')
#plt.title('accuracy of two methods h=0.01256')
plt.legend()
plt.grid()
plt.show()
        
        

    
    
    
