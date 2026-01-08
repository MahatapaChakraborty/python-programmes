#BHZ model 4X4 with offdiagonal unit matrices
import numpy as np
import matplotlib.pyplot as plt
N=100
kx=np.linspace(0,2*np.pi,N)
ky=np.linspace(0,2*np.pi,N)
Kx,Ky=np.meshgrid(kx,ky)
m=0
def eigenvalues(kx,ky,m):
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
    #E=np.array((E.real).sort())
    Ereal=[]
    for i in range(len(E)):
        Ereal.append((E[i]).real)

    Ereal.sort()
    return(Ereal)
#print(eigenvalues(5,5,0))
eigval0=np.zeros((len(kx),len(ky)), dtype=complex)
for i in range(len(kx)):
    for j in range(len(ky)):
        eigval0[i,j]=(eigenvalues(kx[i],ky[j],m)[0])

eigval1=np.zeros((len(kx),len(ky)), dtype=complex)
for i in range(len(kx)):
    for j in range(len(ky)):
        eigval1[i,j]=(eigenvalues(kx[i],ky[j],m)[1])

eigval2=np.zeros((len(kx),len(ky)), dtype=complex)
for i in range(len(kx)):
    for j in range(len(ky)):
        eigval2[i,j]=(eigenvalues(kx[i],ky[j],m)[2])


eigval3=np.zeros((len(kx),len(ky)), dtype=complex)
for i in range(len(kx)):
    for j in range(len(ky)):
        eigval3[i,j]=(eigenvalues(kx[i],ky[j],m)[3])


fig=plt.figure(figsize=(10,7))
ax=fig.add_subplot(111,projection='3d')

ax.plot_surface(Kx,Ky,eigval0,color='red')
ax.plot_surface(Kx,Ky,eigval1,color='blue')
ax.plot_surface(Kx,Ky,eigval2,color='green')
ax.plot_surface(Kx,Ky,eigval3,color='yellow')
ax.set_title('4X4 BHZ model with offdiagonal unit matrix band structure m=0')

ax.set_xlabel('Kx')
ax.set_ylabel('Ky')
ax.set_zlabel('Ek')

plt.show()

