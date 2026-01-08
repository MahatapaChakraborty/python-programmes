#Numerical calculation of Berry phase for a particle in magnetic field
#|up ncap>=(cos theta/2 sin theta/2*e^iphi)

import numpy as np

def upeigenvector(theta,phi):
    a=np.cos(theta/2)
    b=np.sin(theta/2)*np.exp(1j*phi)
    up=np.array([a,b])
    return(up)

#print(upeigenvector(0,np.pi/2))
N=700
theta1=np.linspace(0,np.pi/2,N)
#print(theta1)
phi1=np.linspace(0,np.pi/2,N)
#print(phi1[0],phi1[N])
prod=1
for i in range(N-1):
    a=np.vdot(upeigenvector(theta1[i],0),upeigenvector(theta1[i+1],0))
    prod=prod*a#from north pole to x=1,y=0,z=0

#print(prod)
m=prod
for i in range(1,N-1):
    b=np.vdot(upeigenvector(np.pi/2,phi1[i]),upeigenvector(np.pi/2,phi1[i+1]))
    m=m*b#from x=1,y=0,z=0 to x=0,y=1,z=0

#print(m)
p=m                 
phi1=np.linspace(np.pi/2,0,N)

for i in range(N-1):
    c=np.vdot(upeigenvector(phi1[i],np.pi/2),upeigenvector(phi1[i+1],np.pi/2))
    p=p*c#from x=0,y=1,z=0 to x=0,y=0,z=1 aka the north pole

#print(p)
totalprod=p
lntotalprod=np.log(p)
berryphase=-lntotalprod.imag
print(berryphase/np.pi)
#result -0.24964234620886974 theoritically expected -0.25

    
