import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path


#lattice basis vectors
a1=np.array([0.5,np.sqrt(3)/2])
a2=np.array([-0.5,np.sqrt(3)/2])

#reciprocal lattice vectors
b1=np.array([2*np.pi ,2*np.pi/3])
b2=np.array([-2*np.pi ,2*np.pi/np.sqrt(3)])

#BZ vertices #all three work
#vertices=np.array([(b1+b2)/3,-(b1+b2)/3,(2*b1-b2)/3,-(2*b1-b2)/3,(b1-2*b2)/3,-(b1-2*b2)/3])
vertices=np.array([[-2.53699287, -2.86099692],[ 2.53699287,  2.86099692],[-2.53699287,  2.86099692],[-3.74619244,  0.76660181],[ 2.53699287, -2.86099692],[ 3.74619244, -0.76660181]])
#vertices=np.array([[2*np.pi/3,2*np.pi/np.sqrt(3)],
                   #[-2*np.pi/3,2*np.pi/np.sqrt(3)],
                   #[2*np.pi/3,-2*np.pi/np.sqrt(3)],
                   #[-2*np.pi/3,-2*np.pi/np.sqrt(3)],
                   #[4*np.pi/3,0],[-4*np.pi/3,0]])
                   
#plt.scatter(vertices[:,0],vertices[:,1])
#plt.show()
#print(vertices)
#print(vertices.shape)
BZboundary=Path(vertices)#creates a shape from that path
#latticegrid
#kx=np.linspace(-4*np.pi/3 -1,4*np.pi/3 + 1,100)
#ky=np.linspace(-2*np.pi/np.sqrt(3)-1,2*np.pi/np.sqrt(3)+1,100)
kxmax,kymax=vertices.max(axis=0)
kxmin,kymin=vertices.min(axis=0)
kx=np.linspace(kxmin -2,kxmax + 2,200)
ky=np.linspace(kymin-2,kymax+2,200)
Kx,Ky=np.meshgrid(kx,ky)

#reclatticepoints=np.column_stack([Kx.ravel(),Ky.ravel()])
#from the mesh grid makes array of all possible 2D points
#print(reclatticepoints)

#insideBZ=BZboundary.contains_points(reclatticepoints)
#gives boolean array (true/false) if the reclatticepoints are inside BZ or not
#print(insideBZ)

#kxin=reclatticepoints[insideBZ,0]#picks out the xvalues from
#reclatticepoints when insideBZ true
#kyin=reclatticepoints[insideBZ,1]
#print(kxin)
#kxinordered=np.array(kxin.sort())
#kyinordered=np.array(kyin.sort())

Kxin,Kyin=np.meshgrid(kx,ky)

k1=Kx*(a1[0])+Ky*(a1[1])
k2=Kx*(a2[0])+Ky*(a2[1])

Jx=1
Jy=1
Jz=1

def refk(k1,k2,Jx,Jy,Jz):
    refk=Jz+Jx*np.cos(k1)+Jy*np.cos(k2)
    return(refk)

def Imfk(k1,k2,Jx,Jy,Jz):
    imfk=Jx*np.sin(k1)+Jy*np.sin(k2)
    return(imfk)

def xy(k1,k2,Jx,Jy,Jz):
    x=refk(k1,k2,Jx,Jy,Jz)
    y=Imfk(k1,k2,Jx,Jy,Jz)
    xy=np.array([x,y])
    return xy


def real(i,j):
    real=xy(k1[i,j],k2[i,j],Jx,Jy,Jz)[0]
    return(real)

def im(i,j):
    im=xy(k1[i,j],k2[i,j],Jx,Jy,Jz)[1]
    return(im)

def Hk(g,h):
    hk=np.array([[0,real(g,h)+(im(g,h))*(1j)],[real(g,h)-(im(g,h))*(1j),0]])
    return(hk)

def eigenvaluesgreater(g,h):
    eigenvalue=np.linalg.eigvals(Hk(g,h))[0]
    return(eigenvalue)

def eigenvalueslesser(g,h):
    eigenvalue1=np.linalg.eigvals(Hk(g,h))[1]
    return(eigenvalue1)

#print(eigenvalueslesser(9,7))


greatereigval=np.zeros((len(kx),len(ky)), dtype=complex)
for i in range(len(kx)):
    for j in range(len(ky)):
        greatereigval[i,j]=eigenvaluesgreater(i,j)

lessereigval=np.zeros((len(kx),len(ky)), dtype=complex)
for i in range(len(kx)):
    for j in range(len(ky)):
        lessereigval[i,j]=eigenvalueslesser(i,j)


#print(greatereigval)

#print(seteiggreatervalues)

#fig=plt.figure(figsize=(10,7))
#ax=fig.add_subplot(111,projection='3d')

#ax.plot_surface(Kxin,Kyin,greatereigval-lessereigval,color='red')
#ax.plot_surface(Kxin,Kyin,lessereigval,color='blue')

#ax.set_xlabel('Kx')
#ax.set_ylabel('Ky')
#ax.set_zlabel('Ek')

#plt.show()
#delE=np.column_stack((Kxin.ravel(),Kyin.ravel(),(np.real(greatereigval-lessereigval)).ravel()))
#print(delE)
#np.savetxt("delEdata",delE)
#import os
#print(os.getcwd())

plt.figure(figsize=(10,7))
plt.pcolormesh(Kxin,Kyin,np.real(greatereigval-lessereigval),shading='auto',cmap='coolwarm')
plt.colorbar(label="Difference between two bands")
plt.xlabel("k_x values")
plt.ylabel("k_y values")
plt.show()





