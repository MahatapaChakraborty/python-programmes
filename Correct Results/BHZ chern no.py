import numpy as np
import matplotlib.pyplot as plt

def chernnoBHZ(N, m):
    kmin, kmax = -np.pi, np.pi
    kxlist = np.linspace(kmin, kmax, N, endpoint=False)#Creates array of equidistant
    #points including kmax and excluding kmin
    kylist = np.linspace(kmin, kmax, N, endpoint=False)

    #Returns lower energy eigen vector
    #at a given BZ point
    def eigenvectlowerH(kx, ky, m):
        H12 = np.sin(kx) - 1j*np.sin(ky)
        H21 = np.sin(kx) + 1j*np.sin(ky)
        H11 = m + np.cos(kx) + np.cos(ky)
        H22 = -H11
        H = np.array([[H11, H12], [H21, H22]])
        eigval, eigvect = np.linalg.eig(H)
        idx = np.argmin(eigval)          
        v = eigvect[:, idx]
        return v / np.linalg.norm(v)     

    sumsq = 0.0

    for i in range(N):
        for j in range(N):

            A = eigenvectlowerH(kxlist[i], kylist[j], m)
            B = eigenvectlowerH(kxlist[i], kylist[(j+1) % N], m)
            C = eigenvectlowerH(kxlist[(i+1) % N], kylist[(j+1)%N], m)
            D = eigenvectlowerH(kxlist[(i+1)%N], kylist[j], m)#j runs from 0 to N-1, total N steps
            #the Nth value and 0th value in discretized kx axis must be one and the same as both kx and ky
            #wraps around themselves

            W = np.vdot(A, B) * np.vdot(B, C) * np.vdot(C, D) * np.vdot(D, A)

            phase = np.log(W / np.abs(W))   # = i*theta
            sumsq += -phase.imag

    ch = sumsq / (2*np.pi)
    return np.round(ch)

#print(chernnoBHZ(101, -1))
M=np.linspace(-5,5,100)
#plt.plot(M,chernnoBHZ(30,M))
C=[]
for i in range(len(M)):
    C.append(chernnoBHZ(25,M[i]))

plt.xlabel('m (100 points between +-5)')
plt.ylabel('Chern Number')
plt.title('BHZ Model Variation of Chern No. (25X25 grid taken)')
plt.plot(M,C)
plt.axhline(y=1,color='b',linestyle='--')
plt.axhline(y=-1,color='b',linestyle='--')

plt.show()
