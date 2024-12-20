#d2x/dt2 + x=0 x0=0 z0=1
#dz/dt = -x and dx/dt=z
#z(t0+h) - z(t0) = -h*x(t0) and x(t0+h)-x(t0)=h*z(t0)

import numpy as np
import matplotlib.pyplot as plt

time=[]
tn=2*np.pi
t0=0
n=500
h=(tn-t0)/n

for i in range(n+1):
	t=t0+i*h
	time.append(t)

#print(time)

def eulervel(time,x0,z0):
	euler=[]
	eulervel=[]
	vold=z0
	xold=x0
	
	for i in range(len(time)):
		vnew=vold-h*xold
		xnew=xold+h*vold
		euler.append(xnew)
		eulervel.append(vnew)
		vold=vnew
		xold=xnew
		
	return(eulervel)
	
def euler(time,x0,z0):
	euler=[]
	eulervel=[]
	vold=z0
	xold=x0
	
	for i in range(len(time)):
		vnew=vold-h*xold
		xnew=xold+h*vold
		euler.append(xnew)
		eulervel.append(vnew)
		vold=vnew
		xold=xnew
		
	return(euler)
	
#print(euler(time,0,1))
#plt.plot(time,euler(time,0,1))
#plt.show()

def rk2(time,euler,eulervel,x0,z0):
	rk2=[]
	rk2vel=[]
	xrk2old=x0
	velrk2old=z0
	
	for i in range(len(time)):
		velrk2new=velrk2old-(h/2)*(xrk2old+euler[i])
		xrk2new=xrk2old+(h/2)*(velrk2old+eulervel[i])
		rk2.append(xrk2new)
		rk2vel.append(velrk2new)
		xrk2old=xrk2new
		velrk2old=velrk2new
		
	return(rk2)
	
def rk2vel(time,euler,eulervel,x0,z0):
	rk2=[]
	rk2vel=[]
	xrk2old=x0
	velrk2old=z0
	
	for i in range(len(time)):
		velrk2new=velrk2old-(h/2)*(xrk2old+euler[i])
		xrk2new=xrk2old+(h/2)*(velrk2old+eulervel[i])
		rk2.append(xrk2new)
		rk2vel.append(velrk2new)
		xrk2old=xrk2new
		velrk2old=velrk2new
		
	return(rk2vel)
	
plt.xlabel("time")
plt.ylabel("position")	
plt.plot(time,rk2(time,euler(time,0,1),eulervel(time,0,1),0,1),color="b",label="rk2 solnn")
#plt.plot(time,rk2vel(time,euler(time,0,1),eulervel(time,0,1),0,1),color='k',label='rk2 velocity')
plt.plot(time,euler(time,0,1),color="r",label="euler soln")
#plt.plot(time,eulervel(time,0,1),color='g',label='euler velocity')
plt.plot(time,np.sin(time),color="g",label="original soln")
plt.title('simple harmonic osscilator result using rk2 intial pos 0 intial vel 1 h=0.0125')
plt.grid()
plt.legend()

plt.show()


	
	

