#Haldane model edge spectrum

import numpy as np
import matplotlib.pyplot as plt

t1=3#nearest neighbor interaction
t2=0.2#next nearest neighbor interaction
M=0.15#mass(inversion symmetry breaking) term
phi=np.pi/2

dely=np.sqrt(3)/2
N=120#no. of unit cells taken in x direction

kyvals=np.linspace(-np.pi,np.pi,500)

def edgehamiltonian(ky):

    H=np.zeros((2*N,2*N),dtype=complex)

    tNNright=-2*t1*np.cos(ky*dely)
    tNNleft=-t1

    tNNNb=-2*t2*np.cos(ky*dely-phi)
    tNNNa=-2*t2*np.cos(ky*dely+phi)
    tNNNaa=-4*t2*np.cos(2*ky*dely+phi)
    tNNNbb=-4*t2*np.cos(2*ky*dely-phi)

    for x in range(N):
        A=2*x #a sublattice points
        B=2*x + 1#b sublattice points

        #diagonal terms
        H[A,A]=M+tNNNaa
        H[B,B]=-M+tNNNbb

        #NN terms
        if x<N-1:
            H[A+2,B]=tNNright
            H[B,A+2]=tNNright

        if x>0:
            H[A-2,B]=tNNleft
            H[B,A-2]=tNNleft

        #NNN terms

        if x<N-2:
            H[A+4,A]=tNNNa
            H[A,A+4]=tNNNa
            H[B+4,B]=tNNNb
            H[B,B+4]=tNNNb

    return(H)

#print(edgehamiltonian(0.5))

energies=[]
for ky in kyvals:
    H=edgehamiltonian(ky)
    eigvals=np.linalg.eigvalsh(H)
    energies.append(eigvals)

energies=np.array(energies)


plt.figure(figsize=(7,5))
for i in range(2*N):
    plt.plot(kyvals,energies[:,i],'k',lw=0.5)

plt.xlabel(r"$k_y$")
plt.ylabel("$Energy$")
plt.axvline(2*np.pi/(3*np.sqrt(3)),color='k',linestyle='--')
plt.axvline(-2*np.pi/(3*np.sqrt(3)),color='k',linestyle='--')
plt.title(f"Haldane Model Edge Spectrum (t1={t1},t2={t2},M={M},phi={phi:.2f})")
plt.grid(True)
plt.show()
#this is wrong for tuple counting taken wrong

            

    
