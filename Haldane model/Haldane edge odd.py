import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm



t1 = 1.0#parameters
t2 = 0.24
M  = 0.7
phi = np.pi/2

Ly = 26
N = 2*Ly - 1   

kx_vals = np.linspace(-np.pi, np.pi, 400)




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




def build_H_odd(kx):

    H = np.zeros((N, N), dtype=complex)

    fp = f_plus(kx)
    fm = f_minus(kx)
    gp = g_plus(kx)
    gm = g_minus(kx)
    g0v = g0(kx)

    for n in range(Ly):

        A = 2*n
        B = 2*n + 1

        # -------- A SITE --------
        if A < N:

            H[A, A] = M + fp

            # A-A hopping (g_-)
            if n > 0:
                A_prev = 2*(n-1)
                if A_prev < N:
                    H[A, A_prev] = gm
                    H[A_prev, A] = gm

            if n < Ly-1:
                A_next = 2*(n+1)
                if A_next < N:
                    H[A, A_next] = gm
                    H[A_next, A] = gm

            # A-B (same cell)
            if B < N:
                H[A, B] = g0v
                H[B, A] = g0v

            # A-B (n-1)
            if n > 0:
                B_prev = 2*(n-1) + 1
                if B_prev < N:
                    H[A, B_prev] = t1
                    H[B_prev, A] = t1


        
        if B < N:

            H[B, B] = -M + fm

            # B-B hopping (g_+)
            if n > 0:
                B_prev = 2*(n-1) + 1
                if B_prev < N:
                    H[B, B_prev] = gp
                    H[B_prev, B] = gp

            if n < Ly-1:
                B_next = 2*(n+1) + 1
                if B_next < N:
                    H[B, B_next] = gp
                    H[B_next, B] = gp

            # B-A (n+1)
            if n < Ly-1:
                A_next = 2*(n+1)
                if A_next < N:
                    H[B, A_next] = t1
                    H[A_next, B] = t1

    return H




energies = []
upperedge_weights = []
loweredge_weights = []

for kx in kx_vals:

    H = build_H_odd(kx)
    eigvals, eigvecs = np.linalg.eigh(H)

    energies.append(eigvals)

    # Edge localization (first & last site)
    wU = np.abs(eigvecs[0, :])**2
    wL = np.abs(eigvecs[-1, :])**2

    upperedge_weights.append(wU)
    loweredge_weights.append(wL)

energies = np.array(energies)
upperedge_weights = np.array(upperedge_weights)
loweredge_weights = np.array(loweredge_weights)


# PLOT 1: EDGE SPECTRUM

'''plt.figure(figsize=(6,6))

for i in range(N):
    plt.plot(kx_vals, energies[:, i], 'k', linewidth=0.5)

plt.xlabel(r"$k_x$")
plt.ylabel("Energy")
plt.title(f"Edge Spectrum (Odd N={N},t1={t1},t2={t2},M={M})")
#plt.ylim(-3, 3)
plt.grid(alpha=0.3)
plt.show()'''



# PLOT 2: UPPER EDGE WEIGHT

'''plt.figure(figsize=(6,6))

for i in range(N):
    plt.plot(kx_vals, upperedge_weights[:, i], 'k', linewidth=0.6)

plt.xlabel(r"$k_x$")
plt.ylabel(r"$|\psi(0)|^2$")
plt.title(f"Upper Edge Weight(N={N},t1={t1},t2={t2},M={M},phi={phi:.2f}")
plt.grid(alpha=0.3)
plt.show()'''



# PLOT 3: LOWER EDGE WEIGHT


'''plt.figure(figsize=(6,6))

for i in range(N):
    plt.plot(kx_vals, loweredge_weights[:, i], 'k', linewidth=0.6)

plt.xlabel(r"$k_x$")
plt.ylabel(r"$|\psi(2*Ly-2)|^2$")
plt.title(f"Lower Edge Weight(N={N},t1={t1},t2={t2},M={M},phi={phi:.2f})")
plt.grid(alpha=0.3)
plt.show()'''



# PLOT 4: EDGE LOCALIZATION (ONLY EDGE STATES COLORED)

plt.figure(figsize=(6,6))

# Combine edge weights
edge_weight = upperedge_weights + loweredge_weights

# Global normalization
global_min = edge_weight.min()
global_max = edge_weight.max()
vmin = max(global_min, 1e-10)

norm = LogNorm(vmin=vmin, vmax=global_max)

threshold = 0.2  # the larger the threshold the better the
#for loop will work in choosing the edge bands

for i in range(2*Ly-1):

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

plt.colorbar(label="Edge weight(upper+lower)(lognormalized)")

plt.xlabel(r"$k_x$")
plt.ylabel("Energy")
plt.title(f"Edge Weight(Upper+Lower) (N={N},t1={t1}, t2={t2}, M={M}, phi={phi:.2f})")
plt.show()
