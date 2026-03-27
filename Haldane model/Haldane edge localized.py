import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize


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


# ============================================================
# 🔥 ONLY CHANGE IS HERE (GLOBAL NORMALIZATION)
# ============================================================

plt.figure(figsize=(6,6))

# Global min/max across ALL bands
global_min = upperedge_weights.min()
global_max = upperedge_weights.max()

# Avoid log(0)
vmin = max(global_min, 1e-10)

#norm = Normalize(vmin=vmin, vmax=global_max)
norm = LogNorm(vmin=vmin, vmax=global_max)
for i in range(2*Ly):
    plt.scatter(kx_vals, energies[:, i],
                c=upperedge_weights[:, i],
                cmap='viridis',
                norm=norm,   # <-- global normalization applied
                s=3)

plt.colorbar(label=r"$|\langle \psi(n=0)|eigvecs \rangle|^2$")

plt.xlabel(r"$k_x$")
plt.ylabel("Energy")
plt.title(f"Finite Sample Band Spectrum with Upper Edge Localization(t1={t1},t2={t2},M={M},phi={phi:.2f})")

plt.show()