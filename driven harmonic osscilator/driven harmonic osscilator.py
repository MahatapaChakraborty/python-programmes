#driven harmonic osscilator
#d2x/dt2 + 2*a*dx/dt + (w0^2)*x =f0*cos(wd*t)
#original sol. x=f0/(((wd^2-w0^2)^2 + 4*a^2*wd^2)^0.5) cos(wd*t-phi)<-particular soln.
#               tan(phi)=2*a*wd/(wd^2-w0^2)
#resonat freq(max amp)=wd^2=wR^2=w0^2-2*(a^2)
#dz/dt=(f0*cos(wd*t))-(w0^2)*x-2*a*z and dx/dt=z
#z(n+1)=z(n)+h*(f0*cos(wd*tn)-(w0^2)*xn-2*a*zn) and x(n+1)=x(n)+h*zn 

import numpy as np
import matplotlib.pyplot as plt

time=[]
t0=0#initial time
tn=10#final time
n=10000#total no. of points-1
h=(tn-t0)/n
for i in range(n+1):
    t=t0+i*h
    time.append(t)

#print(time)
#initial conditions
x0=2#<-initial position
z0=0#<-initial velocity
w0=10#<-natural frequency
#wd=20#<-driven frequency
a=4#<-damping coeff/m
f0=50

def driven(a,z,w0,x,f0,wd,t):
    p1=f0*(np.cos(wd*t))
    p2=(w0**2)*x
    p3=2*a*z
    return(p1-p2-p3)

def xsoleuler(time,x0,z0,a,w0,f0,wd):
    xsoleuler=[]
    veuler=[]
    xinitial=x0
    vinitial=z0
    xold=xinitial
    vold=vinitial
    for i in range(len(time)):
        xnew=xold+h*vold
        vnew=vold+h*driven(a,vold,w0,xold,f0,wd,time[i])
        veuler.append(vnew)
        xsoleuler.append(xnew)
        xold=xnew
        vold=vnew

    return(xsoleuler)

#print(xsoleuler(time,x0,z0,a,w0,f0,wd))

def veleuler(time,x0,z0,a,w0,f0,wd):
    xsoleuler=[]
    veuler=[]
    xinitial=x0
    vinitial=z0
    xold=xinitial
    vold=vinitial
    for i in range(len(time)):
        xnew=xold+h*vold
        vnew=vold+h*driven(a,vold,w0,xold,f0,wd,time[i])
        veuler.append(vnew)
        xsoleuler.append(xnew)
        xold=xnew
        vold=vnew

    return(veuler)

#plt.xlabel('time(t)')
#plt.ylabel('position(x)')
#plt.plot(time,xsoleuler(time,x0,z0,a,10,f0,4),color='b',label='wd=4')
#plt.plot(time,xsoleuler(time,x0,z0,a,10,f0,10),color='g',label='wd=10')
#plt.plot(time,xsoleuler(time,x0,z0,a,10,f0,8.24),color='r',label='wd=8.24 resonant freq')
#plt.plot(time,xsoleuler(time,x0,z0,a,10,f0,13),color='k',label='wd=13')
#plt.title('d2x/dt2+2*a*dx/dt+(w0^2)*x=f0*cos(wd*t),x0=2,z0=0,a=4,f0=50,w0=10,h=0.001')
#plt.legend()
#plt.grid()
#plt.show()

#resonence curve
frequencies=np.linspace(1,15,250)

def rmsamplitude(frequencies,x0,z0,a,w0,f0):
    rmsamplitudes=[]
    for i in range(len(frequencies)):
        amp=xsoleuler(time,x0,z0,a,w0,f0,frequencies[i])
        ampavg=0
        for j in range(len(amp)):
            ampavg+=(amp[j])**2

        rmsamp=ampavg**0.5
        
        
        rmsamplitudes.append(rmsamp)

    return(rmsamplitudes)

plt.xlabel('driving frequencies')
plt.ylabel('rms amp at that frequency')
plt.plot(frequencies,rmsamplitude(frequencies,x0,z0,4,10,f0),color='k',label='w0=10 a=4')
plt.plot(frequencies,rmsamplitude(frequencies,x0,z0,3.07,10,f0),color='r',label='w0=10 a=3.07')
plt.plot(frequencies,rmsamplitude(frequencies,x0,z0,5,10,f0),color='b',label='w0=10 a=5')
plt.title('resonence curve for two different resonant frequencies hw=0.056 ht=0.001')
plt.grid()
plt.legend()
plt.show()
        
        

    
