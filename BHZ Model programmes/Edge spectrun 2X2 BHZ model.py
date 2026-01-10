import numpy as np
import matplotlib.pyplot as plt

Nx=100# no. of sites in x-direction
Nk=200# we will plot energies wrt ky this
     # is ky line reselution
m=4

kypoints=np.linspace(0,2*np.pi,Nk)

sigmaz=np.array([[1,0],[0,-1]])
sigmax=np.array([[0,1],[1,0]])
sigmay=np.array([[0,-1j],[1j,0]])

H=np.zeros((2*Nx,2*Nx),dtype=complex)
energies=np.zeros((Nk,2*Nx),dtype=complex)

for i,ky in enumerate(kypoints):
    H1=(m+2*np.cos(ky))*sigmaz - np.sin(ky)*sigmay
    H2=sigmaz - 0.5j*sigmax

    for x in range(Nx):
        H[2*x:2*x+2,2*x:2*x+2]=H1
        #sets the diagonal 2X2
        # matrices in the larger H matrix

        if x<Nx-1:
            H[2*x:2*x+2,2*(x+1):2*(x+1)+2]=H2#for each diagonal H1
            #the next two coloumns the same two rows has an H2
            H[2*(x+1):2*(x+1)+2,2*x:2*x+2]=H2.conj().T#for each diagonal H1
            #in that same coulmn the lower 2 rows have a H2dagger

    energies[i]=np.linalg.eigvalsh(H) #each specific row has stored energies corresponding
    # to a specific ky point and all x points

plt.figure(figsize=(6,4))
for n in range(2*Nx):
    plt.plot(kypoints,energies[:,n],'k',lw=0.6)

plt.xlabel("$K_y$")
plt.ylabel("Energies")
plt.title("Edge spectrum 2X2 BHZ model m=0")
plt.show()
                                 

                                 

        
        



