#d2x/dt2 + 2*a*(dx/dt) + (w^2)*x=0
#z=dx/dt(velocity) and dz/dt + 2*a*z + (w^2)*x=0
#x(t0+h)-x(t0)=h*z(t0) and z(t0+h)-z(t0)=-h*(2*a*z(t0) + (w^2)*(x(t0)))
#x(t1)=x(t0)+h*z(t0) and z(t1)=z(t0)-h*(2*a*z(t0) + (w^2)*(x(t0))) and so on
#the iteretive formula is
#                   x(n+1)=x(n)+h*z(n) and z(n+1)=z(n)-h*(2*a*z(n) + (w^2)*(x(n)))

import numpy as np
import matplotlib.pyplot as plt

t0=0#<-initial time
tn=(2.5)*np.pi#<-final time
n=10000#<-total no. of points-1
h=(tn-t0)/n
time=[]
for i in range(n+1):
    t=t0+i*h
    time.append(t)
#print(time)

def xdiffsoln(a,w,x0,z0,tt):
    xdiffsoln=[]
    velsoln=[]
    xsolnold=x0
    velsolnold=z0

    for i in range(len(tt)):
        xsolnnew=xsolnold+h*velsolnold
        velsolnnew=velsolnold-(h*((2*a*velsolnold)+(w**2)*(xsolnold)))
        xdiffsoln.append(xsolnnew)
        velsoln.append(velsolnnew)
        xsolnold=xsolnnew
        velsolnold=velsolnnew

    return(xdiffsoln)

def velsoln(a,w,x0,z0,tt):
    xdiffsoln=[]
    velsoln=[]
    xsolnold=x0
    velsolnold=z0

    for i in range(len(tt)):
        xsolnnew=xsolnold+h*velsolnold
        velsolnnew=velsolnold-(h*((2*a*velsolnold)+(w**2)*(xsolnold)))
        xdiffsoln.append(xsolnnew)
        velsoln.append(velsolnnew)
        xsolnold=xsolnnew
        velsolnold=velsolnnew

    return(velsoln)

#print(xdiffsoln(0,2,1,1,time))
plt.xlabel('time')
plt.ylabel('displacement')
#plt.plot(time,xdiffsoln(0,2,1,0,time),color='r',label='w=2')
plt.plot(time,xdiffsoln(1,6,1,0,time),color='r',label='w=6 a=1 underdamped a^2<w^2')
plt.plot(time,xdiffsoln(6,6,1,0,time),color='b',label='w=6 a=6 crticallydamped a^2=w^2')
plt.plot(time,xdiffsoln(10,6,1,0,time),color='g',label='w=6 a=10 overdamped a^2>w^2')
plt.title('harmonic osscilator d2x/dt2 +2*a*dx/dt+(w^2)x=0 initial vel=0 initial position=1')
plt.legend()
plt.grid()
plt.show()


    
