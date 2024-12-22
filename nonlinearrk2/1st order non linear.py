#dy/dt=((sin(t))^2)*y,y(t=0)=2,original sol=2*exp(0.5*(t-sint*cost))

import numpy as np
import matplotlib.pyplot as plt

y0=2#<-initial y
t0=0#<-initial t
tn=20#<-final t
n=200#<-total no. of points-1
h=(tn-t0)/n#<-step length
time=[]

for i in range(n+1):
    t=t0+i*h
    time.append(t)

#print(time)

def sin2(x):
    sin21=np.sin((x)**2)
    return(sin21)

#print(sin2(45*(np.pi/180)))
def originalsol(xx):
    originalsol=[]
    for i in range(len(xx)):
        originalsol1=2*np.exp(0.5*(xx[i]-(np.sin(xx[i])*np.cos(xx[i]))))
        originalsol.append(originalsol1)
    return(originalsol)


def eulersol(time,y0):
    eulersol=[]
    yold=y0

    for i in range(len(time)):
        ynew=yold+h*((np.sin(time[i]))**2)*yold
        eulersol.append(ynew)
        yold=ynew

    return(eulersol)

#plt.plot(time,eulersol(time,y0))
#plt.show()

def rk2sol(time,eulersol,y0):
    rk2sol=[]
    rk2old=y0

    for i in range(len(eulersol)-1):
        rk2new=rk2old+(h/2)*(((((np.sin(time[i]))**2)*rk2old)+(((np.sin(time[i+1]))**2)*eulersol[i])))
        rk2sol.append(rk2new)
        rk2old=rk2new

    return(rk2sol)

#print(rk2sol(time,eulersol(time,y0),y0))

timenew=[]
for i in range(len(time)-1):
    tnew=time[i]
    timenew.append(tnew)

eulersolnew=[]
for i in range(len(eulersol(time,y0))-1):
    eulersolnew1=eulersol(time,y0)[i]
    eulersolnew.append(eulersolnew1)

plt.xlabel('t')
plt.ylabel('difference from original soln')
#plt.plot(timenew,originalsol(timenew),color='k',label='original soln')
#plt.plot(timenew,eulersolnew,color='r',label='euler soln')
#plt.plot(timenew,rk2sol(time,eulersol(time,y0),y0),color='b',label='rk2 soln')
#plt.title('dy/dt=((sint)^2)*y, y(0)=2,timestep=0.1')
plt.plot(timenew,np.array(eulersolnew)-np.array(originalsol(timenew)),color='r',label='euler')
plt.plot(timenew,np.array(rk2sol(time,eulersol(time,y0),y0))-np.array(originalsol(timenew)),color='b',label='rk2')
plt.text(0,-2000,'rk2 better than euler as difference from original at large time for rk2 is lower')
plt.title('difference between original and method soln. dy/dt=((sint)^2)*y,h=0.1')
plt.legend()
plt.grid()
plt.show()
                             

    



    
    
