import numpy as np
import matplotlib.pyplot as plt
N=200
m=-1
kxpoints=np.linspace(0,2*np.pi,N)
kypoints=np.linspace(0,2*np.pi,N)

Kx,Ky=np.meshgrid(kxpoints,kypoints)

#eigenvalue calculation

def eigenvalueH(kx,ky,m):
    H12=np.sin(kx)- (1j)* (np.sin(ky))
    H21=np.sin(kx)+ (1j)* (np.sin(ky))
    H11=m + np.cos(kx) + np.cos(ky)
    H22=-H11

    H=np.array([[H11, H12],[H21, H22]])
    E=np.linalg.eigvals(H)
    Ereal=[]
    for i in range(len(E)):
        Ereal.append(E[i].real)
    Ereal.sort()
    return(Ereal)

#print(eigenvalueH(3,3,2))
greatereigval=np.zeros((N,N), dtype=complex)
for i in range(N):
    for j in range(N):
        greatereigval[i,j]=eigenvalueH(kxpoints[i],kypoints[j],m)[1]

lessereigval=np.zeros((N,N), dtype=complex)
for i in range(N):
    for j in range(N):
        lessereigval[i,j]=eigenvalueH(kxpoints[i],kypoints[j],m)[0]
        
fig=plt.figure(figsize=(10,7))
ax=fig.add_subplot(111,projection='3d')

ax.plot_surface(Kx,Ky,greatereigval,color='red')#greater eigenvalue
ax.plot_surface(Kx,Ky,lessereigval,color='blue')#lesser eigenvalue

ax.set_xlabel('$K_x$')
ax.set_ylabel('$K_y$')
ax.set_zlabel('$E_k$')
ax.set_title('Band structure 2X2 BHZ model m=-1')

plt.show()
#plt.figure(figsize=(10,7))
#plt.pcolormesh(Kx,Ky,np.real(greatereigval-lessereigval),shading='auto',cmap='coolwarm')
#plt.colorbar(label="Difference between two bands")
#plt.xlabel("k_x values")
#plt.ylabel("k_y values")
#plt.title("m=0 difference between bands BHZ 2X2 model")
#plt.show()
