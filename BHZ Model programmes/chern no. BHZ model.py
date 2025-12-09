#Chern no. BHZ model
#import math
import numpy as np
import matplotlib.pyplot as plt
def chernnoBHZ(N,m):
    kx=np.linspace(-np.pi,np.pi,N)
    ky=np.linspace(-np.pi,np.pi,N)
    dkx=kx[1]-ky[0]
    dky=ky[1]-ky[0]

    def eigenvectlowerH(kx,ky,m):
        H12=np.sin(kx)- (1j)* (np.sin(ky))
        H21=np.sin(kx)+ (1j)* (np.sin(ky))
        H11=m + np.cos(kx) + np.cos(ky)
        H22=-H11
        H=np.array([[H11, H12],[H21, H22]])
        eigval,eigvect=np.linalg.eig(H)
        return(eigvect[1]/np.sqrt(np.linalg.norm(eigvect[1])))

    sumofallsqrs=0
    for i in range(N-1):
        A=eigenvectlowerH(kx[i],ky[i],m)
        B=eigenvectlowerH(kx[i],ky[i+1],m)
        C=eigenvectlowerH(kx[i+1],ky[i+1],m)
        D=eigenvectlowerH(kx[i+1],ky[i],m)
        prodardsqr=np.vdot(A,B)*np.vdot(B,C)*np.vdot(C,D)*np.vdot(D,A)
        angleprodarsqr=np.angle(prodardsqr)
        #modulo2piangleprodarsqr=angleprodarsqr % (2*(np.pi))
        minusangleprodarsqr=-angleprodarsqr
        sumofallsqrs+=minusangleprodarsqr

    #modulo2pisum=sumofallsqrs % (2*(np.pi))

    #C=modulo2pisum/(dkx*dky)
    return(np.round(sumofallsqrs/(2*(np.pi))))

    
M=np.linspace(-5,5,100)
#plt.plot(M,chernnoBHZ(30,M))
C=[]
for i in range(len(M)):
    C.append(chernnoBHZ(100,M[i]))

plt.xlabel('m')
plt.ylabel('Chern Number')
plt.title('BHZ Model Variation of Chern No. (100X100 grid taken)')
plt.plot(M,C)

plt.show()

#print(chernnoBHZ(300,2))

    
    
    
