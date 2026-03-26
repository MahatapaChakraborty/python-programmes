import numpy as np
import matplotlib.pyplot as plt


t1 = 1.0
t2 = 1
M  = 4.2
phi = np.pi/2

Ly = 25   # number of unit cells along y
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
        if n < Ly-1:
            H[A, 2*(n+1)] = gm

        
        if n > 0:
            H[B, 2*(n-1)+1] = gp
        if n < Ly-1:
            H[B, 2*(n+1)+1] = gp

        
        if n > 0:
            H[A, 2*(n-1)+1] = t1

        
        if n < Ly-1:
            H[B, 2*(n+1)] = t1

    return H


energies = []

for kx in kx_vals:
    H = build_H(kx)
    eigvals = np.linalg.eigvalsh(H)
    energies.append(eigvals)

energies = np.array(energies)


plt.figure(figsize=(6,6))

for i in range(2*Ly):
    plt.plot(kx_vals, energies[:, i], 'k', linewidth=0.5)

plt.xlabel(r"$k_x$")
plt.ylabel("Energy")
#plt.axvline(4*np.pi/(3*np.sqrt(3)),color='k',linestyle='--',lw=0.8)
#plt.axvline(-4*np.pi/(3*np.sqrt(3)),color='k',linestyle='--',lw=0.8)
plt.title(f"Haldane Model Edge Spectrum (Zigzag),t1={t1},t2={t2},M={M},phi={phi:.2f}")
#plt.ylim(-3, 3)
plt.grid(alpha=0.3)

plt.show()
