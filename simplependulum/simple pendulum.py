#solving d2y/dt2 =-(A^2).siny with initial conditions y(t=0)=y0 and dy/dt(t=0)=0(=z0)
#we break the differential eqn. in two simultaneous first order equations
#z=dy/dt<-angular velocity, and dz/dt=-(A^2)*sin(y)
#(y(t0+h)-y(t0))/h =z(t0) and (z(t0+h)-z(t0))/h = -(A^2)*sin(y(t0))
#so y(t1)=y(t0)+h*z(t0) and z(t1)=z(t0)-(A^2)*h*sin(y(t0)) now
#dy/dt(t=t1)=z(t=t1) and dz/dt(t=t1)=-(A^2)*sin(y(t1))
#(y(t1+h)-y(t1))=h*z(t1) and z(t1+h)-z(t1)=-(A^2)*h*sin(y(t1)) so y(t2)=y(t1)+h*z(t1) and z(t2)=z(t1)-(A^2)*h*sin(y(t1)) and so on
#the iterative formula is
#                       y(t(n+1))=y(t(n))+h*z(tn) and z(t(n+1))=z(tn)-(A^2)*h*sin(y(tn))

import matplotlib.pyplot as plt
import numpy as np

#initial conditions
#z0=0#initial velocity
#y0=5#deg)#initial angle
A=4.9

def radian(x):
    return(x*(np.pi)/180)

t0=0#initial time
tn=4*(np.pi)#final time
n=10000#total no. of points-1
h=(tn-t0)/n
time=[]
for i in range(n+1):
    t=t0+i*h
    time.append(t)

#print(time)

def sin(x):
    return(np.sin(x))

def resultangles(xx,z0,y0):
    initialangle=radian(y0)
    initialvel=z0
    resultangles=[]
    angleold=initialangle
    resultangleold=initialangle
    resultvelocities=[]
    velold=initialvel

    for i in range(len(xx)):                        #python enters the loop with old variables as intialvalues and after each itaration
                                                    #it stores the new values at resvel[] and resangl[] and then the newvalues of current
                                                    #itteration becomes the old values for the next itteration
        anglenew=angleold + h*velold
        velnew=velold+(-(A**2))*h*sin(angleold)
        resultangles.append(anglenew)
        resultvelocities.append(velnew)
        velold=velnew
        angleold=anglenew

    return(resultangles)

#print(resultangles(time,0,5))

def resultvelocities(xx,z0,y0):
    initialangle=radian(y0)
    initialvel=z0
    resultangles=[]
    angleold=initialangle
    resultangleold=initialangle
    resultvelocities=[]
    velold=initialvel

    for i in range(len(xx)):                          
        anglenew=angleold + h*velold
        velnew=velold+(-(A**2))*h*sin(angleold)
        resultangles.append(anglenew)
        resultvelocities.append(velnew)
        velold=velnew
        angleold=anglenew

    return(resultvelocities)

plt.xlabel('time(s)')
plt.ylabel('angular displacement(radians)')
plt.plot(time,resultangles(time,0,5),color='r',label='angular displacement(rdians)')
plt.plot(time,resultvelocities(time,0,5),color='b',label='angular velocity(rad/s)')
#plt.plot(time,resultangles(time,0,30),color='g',label='initial angle 30deg')
#plt.plot(time,resultangles(time,0,90),color='b',label='initial angle 90deg')
#plt.plot(time,resultangles(time,0,120),color='m',label='initial angle 120deg')
#plt.plot(time,resultangles(time,0,179),color='k',label='initial angle 179deg')
plt.title('angular displacement vs time and angular velocity vs time with g=9.8m/s^2 and L=2m')
plt.legend()
plt.grid()
plt.show()

        
        
    






