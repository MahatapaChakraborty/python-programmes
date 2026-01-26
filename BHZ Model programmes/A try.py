import numpy as np
import matplotlib.pyplot as plt

Nk = 50
m = 3.0

kx = np.linspace(0, 2*np.pi, Nk)#, endpoint=False)
ky = np.linspace(0, 2*np.pi, Nk)# , endpoint=False)
#kx=np.append(kx,0)
#ky=np.append(ky,0)

KX, KY = np.meshgrid(kx, ky,indexing='ij')#this sets kx as rows (axis=0)
# and ky as coulomns(axis=1)


def eigenvectH(kx, ky, m):
    H12 = np.sin(kx) - 1j*np.sin(ky)
    H21 = np.sin(kx) + 1j*np.sin(ky)
    H11 = m + np.cos(kx) + np.cos(ky)
    H22 = -H11

    H = np.array([[H11, H12],
                  [H21, H22]], dtype=complex)

    eigval, eigvec = np.linalg.eig(H)
    idx = np.argsort(eigval.real)      
    v = eigvec[:, idx[0]]
    #v=np.array([1,v[1]/v[0]])

    return v #/ np.linalg.norm(v)


u = np.zeros((Nk, Nk, 2), dtype=complex)

for i in range(Nk):
    for j in range(Nk):
        u[i, j] = eigenvectH(kx[i], ky[j], m)
dukx = np.gradient(u, kx, axis=0)
duky = np.gradient(u, ky, axis=1)

Ax = np.zeros((Nk, Nk))
Ay = np.zeros((Nk, Nk))

for i in range(Nk):
    for j in range(Nk):
        Ax[i, j] = np.imag(np.vdot(u[i, j], dukx[i, j]))
        Ay[i, j] = np.imag(np.vdot(u[i, j], duky[i, j]))
        p=(Ax[i, j]**2 + Ay[i, j]**2)**0.5
        Ax[i, j]=Ax[i, j]/(p)
        Ay[i, j]=Ay[i, j]/(p)

def eigenvectHanalytic(kx, ky, m):
    v1=-(np.sin(kx)-1j*np.sin(ky))
    v2=m+np.cos(kx)+np.cos(ky)
    v3=v2**2 +((np.sin(kx))**2 + (np.sin(ky))**2)
    v=np.array([v1/(v2 + v3**0.5) ,1])
    v1anothergauge=-v1
    v2anothergauge=v2 - v3**0.5
    vanothergauge=np.array([1,v1anothergauge/v2anothergauge])

    return v #/ np.linalg.norm(v)


uanalytic = np.zeros((Nk, Nk, 2), dtype=complex)

for i in range(Nk):
    for j in range(Nk):
        uanalytic[i, j] = eigenvectHanalytic(kx[i], ky[j], m)
duanalytickx = np.gradient(uanalytic, kx, axis=0)
duanalyticky = np.gradient(uanalytic, ky, axis=1)

Aanalyticx = np.zeros((Nk, Nk))
Aanalyticy = np.zeros((Nk, Nk))

for i in range(Nk):
    for j in range(Nk):
        Aanalyticx[i, j] = np.imag(np.vdot(uanalytic[i, j], duanalytickx[i, j]))
        Aanalyticy[i, j] = np.imag(np.vdot(uanalytic[i, j], duanalyticky[i, j]))
        panalytic=(Aanalyticx[i, j]**2 + Aanalyticy[i, j]**2)**0.5
        Aanalyticx[i, j]=Aanalyticx[i, j]/(panalytic)
        Aanalyticy[i, j]=Aanalyticy[i, j]/(panalytic)
        
        
        
plt.figure(figsize=(6, 6))
plt.quiver(KX, KY, Ax, Ay,label='numerical result',color='b', scale=30, width=0.0009)
plt.quiver(KX, KY, Aanalyticx, Aanalyticy,label='analytical result',color='r', scale=30, width=0.0009)
plt.legend()
plt.xlabel(r"$k_x$")
plt.ylabel(r"$k_y$")
plt.title("Berry Potential Field m=2")
plt.tight_layout()
plt.show()
