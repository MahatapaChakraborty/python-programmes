#dispersion relation for kitaev model
#E(kx,ky)=((Jz+Jxcos(kvec.M2)+ Jycos(kvec.M1))^2 + (Jxsin(kvec.M2)+ Jysin(kvec.M1))^2)^0.5
#M1 and M2 are real lattice basis vectors in direction of y and x bonds.
#M2=(sqrt3/2, 3/2) M1=(-sqrt3/2,3/2)
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def Ek(x,y,jx,jy,jz):
    k1= -(np.sqrt(3))*(x/2) + 3*(y/2)
    k2=(np.sqrt(3))*(x/2) + 3*(y/2)
    sink1=np.sin(k1)
    sink2=np.sin(k2)
    cosk2=np.cos(k2)
    cosk1=np.cos(k1)

    A=jz+jx*cosk2 + jy*cosk1
    B=jx*sink2 + jy*sink1
    c=A**2 + B**2
    return(c**0.5)

jx=1
jy=1
jz=2.5

x=np.linspace(-np.pi, np.pi ,2000)
y=np.linspace(-np.pi/1,np.pi/1,2000)
X,Y=np.meshgrid(x,y)
ek=Ek(X,Y,jx,jy,jz)
eneg=-ek

fig=plt.figure(figsize=(12,8))
ax=fig.add_subplot(111,projection='3d')

ax.plot_surface(X,Y, ek ,color='red')
ax.plot_surface(X,Y, eneg ,color='blue')

ax.set_xlabel("$k_x$")
ax.set_ylabel("$k_y$")
ax.set_zlabel("$Ek$")
plt.show()


    
