import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize


# PARAMETERS 

t1 = 1.0
t2 = 0.3333
M  = 1
phi = np.pi/3

Ly = 25
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
# COLOR ONLY EDGE BANDS (CORRECT APPROACH)
# ============================================================

plt.figure(figsize=(6,6))

# Combine edge weights
edge_weight = upperedge_weights 

# Global normalization
global_min = edge_weight.min()
global_max = edge_weight.max()
vmin = max(global_min, 1e-10)

norm = LogNorm(vmin=vmin, vmax=global_max)

threshold = 0.5   # the larger the threshold the better the
#for loop will work in choosing the edge bands

for i in range(2*Ly):

    #Deciding if this BAND is an edge band
    max_weight = np.max(edge_weight[:, i])

    if max_weight > threshold:#this will just
        #color code those bands whose maximum possible
        #is greater than the threshold
        # EDGE BAND → color entire band(cause only
        #edge bands can have large edge weights)
        plt.scatter(kx_vals,
                    energies[:, i],
                    c=edge_weight[:, i],
                    cmap='viridis',
                    norm=norm,
                    s=5)
    else:
        # BULK BAND → plain black
        plt.plot(kx_vals,
                 energies[:, i],
                 color='black',
                 linewidth=0.7)

plt.colorbar(label="Edge weight(lognormalized")

plt.xlabel(r"$k_x$")
plt.ylabel("Energy")
plt.title(f"Edge vs Bulk Bands (t1={t1}, t2={t2}, M={M}, phi={phi:.2f})")
plt.show()
