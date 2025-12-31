#2D BHZ model Berry phase plotted on each plaquette

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

N=50#creates an (N-1)^2 grid
m=0
Kx=np.linspace(0,2*np.pi,N,endpoint=False)
Ky=np.linspace(0,2*np.pi,N,endpoint=False)

#plaquette boundaries
plaquettes=[]
for i in range(N-1):
    for j in range(N-1):
        oneplaquette=[(Kx[i],Ky[j]),(Kx[i+1],Ky[j]),(Kx[i+1],Ky[j+1]),(Kx[i],Ky[j+1])]
        plaquettes.append(oneplaquette)

p=plaquettes
#print(p)

def eigenvectorlowerH(kx,ky,m):
    H12 = np.sin(kx) - 1j*np.sin(ky)
    H21 = np.sin(kx) + 1j*np.sin(ky)
    H11 = m + np.cos(kx) + np.cos(ky)
    H22 = -H11
    H = np.array([[H11, H12], [H21, H22]])
    eigval, eigvect = np.linalg.eig(H)
    idx = np.argmin(eigval)          
    v = eigvect[:, idx]
    return v / np.linalg.norm(v)

def berryphaseoneplaquette(oneplaquette,m):
    A = eigenvectorlowerH(oneplaquette[0][0],oneplaquette[0][1], m)
    B = eigenvectorlowerH(oneplaquette[1][0],oneplaquette[1][1], m)
    C = eigenvectorlowerH(oneplaquette[2][0],oneplaquette[2][1], m)
    D = eigenvectorlowerH(oneplaquette[3][0],oneplaquette[3][1], m)
    W = np.vdot(A, B) * np.vdot(B, C) * np.vdot(C, D) * np.vdot(D, A)
    phase = np.log(W / np.abs(W))
    phase=-phase.imag
    return(phase)

#print(berryphaseoneplaquette([(0,0),(1,0),(1,1),(0,1)],2))

berryphasevalues=[]
for i in range((N-1)**2):
    values=berryphaseoneplaquette(p[i],m)
    berryphasevalues.append(values)

#print(berryphasevalues)
b=berryphasevalues

collection=PolyCollection(p,array=b,cmap='viridis',edgecolors='k')
fig,ax=plt.subplots()
ax.add_collection(collection)
ax.autoscale()
ax.set_xlabel("k_x values")
ax.set_ylabel("k_y values")
ax.set_aspect('equal')

cbar=plt.colorbar(collection,ax=ax)
cbar.set_label("Berry phase value")
ax.set_title("Berry phase over BZ for 2D BHZ model m=0")
plt.show()




    
