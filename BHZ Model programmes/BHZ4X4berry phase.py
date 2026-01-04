#2D BHZ model Berry phase plotted on each plaquette

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
#brakes the plot into many polygons each is colored for a scaler no.
N=71#creates an (N-1)^2 grid
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

def eigenvectors(kx,ky,m):
    H12=np.sin(kx)- (1j)* (np.sin(ky))
    H21=np.sin(kx)+ (1j)* (np.sin(ky))
    H11=m + np.cos(kx) + np.cos(ky)
    H22=-H11
    H13=0
    H14=0
    H23=0
    H24=0
    H31=0
    H32=0
    H41=0
    H42=0
    H34=np.sin(kx)- (1j)* (np.sin(ky))
    H43=np.sin(kx)+ (1j)* (np.sin(ky))
    H33=m + np.cos(kx) + np.cos(ky)
    H44=-H33

    H=np.array([[H11,H12,H13,H14],
                [H21,H22,H23,H24],
                [H31,H32,H33,H34],
                [H41,H42,H43,H44]])

    E,v=np.linalg.eig(H)
    idx=np.argsort(E)#sorts E and gives their positions from smallest to largest
    V=v[:,idx]
    V0=V[:,2]
    V2=V[:,3]#two upper bands taken
    return(V0/np.linalg.norm(V0) , V2/np.linalg.norm(V2))
    
def berryphaseoneplaquette(oneplaquette,m):
    A0 = eigenvectors(oneplaquette[0][0], oneplaquette[0][1], m)[0]#A,B,C,D diff. points around the plaquette
    A1 = eigenvectors(oneplaquette[0][0], oneplaquette[0][1], m)[1]#0,1 different band state vectors at the same (kx,ky) point
    B0 = eigenvectors(oneplaquette[1][0], oneplaquette[1][1], m)[0]
    B1 = eigenvectors(oneplaquette[1][0], oneplaquette[1][1], m)[1]
    C0 = eigenvectors(oneplaquette[2][0], oneplaquette[2][1], m)[0]
    C1 = eigenvectors(oneplaquette[2][0], oneplaquette[2][1], m)[1]
    D0 = eigenvectors(oneplaquette[3][0], oneplaquette[3][1], m)[0]
    D1 = eigenvectors(oneplaquette[3][0], oneplaquette[3][1], m)[1]
    M12=np.array([[np.vdot(A0,B0),np.vdot(A0,B1)],[np.vdot(A1,B0),np.vdot(A1,B1)]])
    M23=np.array([[np.vdot(B0,C0),np.vdot(B0,C1)],[np.vdot(B1,C0),np.vdot(B1,C1)]])
    M34=np.array([[np.vdot(C0,D0),np.vdot(C0,D1)],[np.vdot(C1,D0),np.vdot(C1,D1)]])
    M41=np.array([[np.vdot(D0,A0),np.vdot(D0,A1)],[np.vdot(D1,A0),np.vdot(D1,A1)]])
    M12M23=np.matmul(M12,M23)
    M34M41=np.matmul(M34,M41)
    M12M23M34M41=np.matmul(M12M23,M34M41)
    detoneplaquette=np.linalg.det(M12M23M34M41)
    lndet=np.log(detoneplaquette/np.abs(detoneplaquette))
    phase=lndet.imag
    return(phase)
#print(berryphaseoneplaquette([(0,0),(1,0),(1,1),(0,1)],2))

berryphasevalues=[]
for i in range((N-1)**2):
    values=berryphaseoneplaquette(p[i],m)
    berryphasevalues.append(values)

#print(berryphasevalues)
b=berryphasevalues

collection=PolyCollection(p,array=b,cmap='coolwarm',edgecolors='k')
fig,ax=plt.subplots()
ax.add_collection(collection)
ax.autoscale()
ax.set_xlabel("k_x values")
ax.set_ylabel("k_y values")
ax.set_aspect('equal')

cbar=plt.colorbar(collection,ax=ax)
cbar.set_label("Berry phase value")
ax.set_title("Berry phase over BZ for 4X4 BHZ model (70X70 grid) m=0")
plt.show()

