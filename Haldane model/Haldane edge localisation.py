import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


# PARAMETERS 

t1 = 1.0
t2 = 0.24
M  = 0.7
phi = np.pi/2

Ly = 20
kx_vals = np.linspace(-np.pi, np.pi, 400)


# FUNCTIONS

def f_plus(kx):
    return 2*t2*np.cos(np.sqrt(3)*kx + phi)

def f_minus(kx):
    return 2*t2*np.cos(np.sqrt(3)*kx - phi)

def g_plus(kx):
    return 2*t2*np.cos(np.sqrt(3)/2 * kx + phi)

def g_minus(kx):
    return 2*t2*np.cos(np.sqrt(3)/2 * kx - phi)

def g0(kx):
    return 2*t1*np.cos(np.sqrt(3)/2 * kx)


# BUILD HAMILTONIAN

def build_H(kx):

    H = np.zeros((2*Ly, 2*Ly), dtype=complex)

    fp = f_plus(kx)
    fm = f_minus(kx)
    gp = g_plus(kx)
    gm = g_minus(kx)
    g0v = g0(kx)

    for n in range(Ly):

        A = 2*n
        B = 2*n + 1

        
        H[A, A] = M + fp
        H[B, B] = -M + fm

        
        H[A, B] = g0v
        H[B, A] = g0v

        
        if n > 0:
            H[A, 2*(n-1)] = gm
            H[2*(n-1), A] = gm

        if n < Ly-1:
            H[A, 2*(n+1)] = gm
            H[2*(n+1), A] = gm

        
        if n > 0:
            H[B, 2*(n-1)+1] = gp
            H[2*(n-1)+1, B] = gp

        if n < Ly-1:
            H[B, 2*(n+1)+1] = gp
            H[2*(n+1)+1, B] = gp

        
        if n > 0:
            H[A, 2*(n-1)+1] = t1
            H[2*(n-1)+1, A] = t1

        
        if n < Ly-1:
            H[B, 2*(n+1)] = t1
            H[2*(n+1), B] = t1

    return H


energies = []
upperedge_weights = []
loweredge_weights = []

for kx in kx_vals:

    H = build_H(kx)

    eigvals, eigvecs = np.linalg.eigh(H)

    energies.append(eigvals)

    # Edge projections
    wU = np.abs(eigvecs[0, :])**2
    wL = np.abs(eigvecs[-1, :])**2

    upperedge_weights.append(wU)
    loweredge_weights.append(wL)

energies = np.array(energies)
upperedge_weights = np.array(upperedge_weights)
loweredge_weights = np.array(loweredge_weights)


#EDGE SPECTRUM

'''plt.figure(figsize=(6,6))

for i in range(2*Ly):
    plt.plot(kx_vals, energies[:, i], 'k', linewidth=0.5)

plt.xlabel(r"$k_x$")
plt.ylabel("Energy")
plt.title("Haldane Model Edge Spectrum (Zigzag)")
plt.ylim(-3, 3)
plt.grid(alpha=0.3)
plt.show()'''


#UPPER EDGE WEIGHT PLOT

'''plt.figure(figsize=(6,6))

for i in range(2*Ly):
    plt.plot(kx_vals, upperedge_weights[:, i],'k', linewidth=0.7)

plt.xlabel(r"$k_x$")
plt.ylabel(r"$|\langle \psi(n=0)|eigvecs \rangle|^2$")
plt.title(f"Upper Edge Weight(Ly={Ly},t1={t1},t2={t2},M={M},phi={phi:.2f})")
plt.grid(alpha=0.3)
plt.show()'''


#LOWER EDGE WEIGHT

plt.figure(figsize=(6,6))

for i in range(2*Ly):
    plt.plot(kx_vals, loweredge_weights[:, i],'k', linewidth=0.7)

plt.xlabel(r"$k_x$")
plt.ylabel(r"$|\langle \psi(n=2Ly-1)|eigvecs \rangle|^2$")
plt.title(f"Lower Edge Weight(Ly={Ly},t1={t1},t2={t2},M={M},phi={phi:.2f})")
plt.grid(alpha=0.3)
plt.show()


# PLOT 4: SPECTRUM COLORED BY UPPER EDGE LOCALIZATION

'''plt.figure(figsize=(6,6))

for i in range(2*Ly):
    plt.scatter(kx_vals, energies[:, i],
                c=upperedge_weights[:, i],
                cmap='viridis',norm=LogNorm(vmin=1e-6,vmax=1), s=3)

plt.colorbar(label=r"$|\langle \psi(n=0)|eigvecs \rangle|^2$")
plt.xlabel(r"$k_x$")
plt.ylabel("Energy")
plt.title(f"Finite Sample Band Spectrum with Upper Edge Localization(t1={t1},t2={t2},M={M},phi={phi:.2f})")
#plt.ylim(-3, 3)
plt.show()'''

