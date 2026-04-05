import numpy as np
import matplotlib.pyplot as plt


# PARAMETERS

t = 1.0
lambda_SO = 0.3
lambda_R =0.1
lambda_v = 0.1
a = 1.0


# PAULI MATRICES

sigma_x = np.array([[0,1],[1,0]], dtype=complex)
sigma_y = np.array([[0,-1j],[1j,0]], dtype=complex)
sigma_z = np.array([[1,0],[0,-1]], dtype=complex)
I2 = np.eye(2)

s_x = sigma_x.copy()
s_y = sigma_y.copy()
s_z = sigma_z.copy()

def kron(a, b):
    return np.kron(a, b)


# DIRAC MATRICES

Gamma1 = kron(sigma_x, I2)
Gamma2 = kron(sigma_z, I2)
Gamma3 = kron(sigma_y, s_x)
Gamma4 = kron(sigma_y, s_y)
Gamma5 = kron(sigma_y, s_z)

def comm(A, B):
    return (A @ B - B @ A) / (2j)

Gamma12 = comm(Gamma1, Gamma2)
Gamma15 = comm(Gamma1, Gamma5)
Gamma23 = comm(Gamma2, Gamma3)
Gamma24 = comm(Gamma2, Gamma4)


# HAMILTONIAN

def H_k(kx, ky):
    x = kx * a / 2
    y = np.sqrt(3) * ky * a / 2

    d1 = t * (1 + 2*np.cos(x)*np.cos(y))
    d2 = lambda_v
    d3 = lambda_R * (1 + np.cos(x)*np.cos(y))
    d4 = np.sqrt(3) * lambda_R * np.sin(x)*np.sin(y)

    d12 = -2*t*np.cos(x)*np.sin(y)
    d15 = lambda_SO * (2*np.sin(2*x) - 4*np.sin(x)*np.cos(y))
    d23 = -lambda_R*np.cos(x)*np.sin(y)
    d24 = np.sqrt(3)*lambda_R*np.sin(x)*np.cos(y)

    H = (d1*Gamma1 + d2*Gamma2 + d3*Gamma3 + d4*Gamma4
         + d12*Gamma12 + d15*Gamma15 + d23*Gamma23 + d24*Gamma24)

    return H


# RECIPROCAL LATTICE

b1 = (2*np.pi/a) * np.array([1, 1/np.sqrt(3)])
b2 = (2*np.pi/a) * np.array([-1, 1/np.sqrt(3)])

Gamma = np.array([0, 0])
K = (b1 - b2) / 3
Kp = (2*b1 + b2)/3
M = b1 / 2

k_path = [Gamma, K, M, Kp, Gamma]
labels = [r'$\Gamma$', r'$K$', r'$M$', r"$K'$", r'$\Gamma$']


# INTERPOLATION

def interpolate_kpath(points, n_per_segment=250):
    k_list = []
    k_dist = [0]
    node_positions = [0]

    for i in range(len(points)-1):
        k_start = points[i]
        k_end = points[i+1]

        for j in range(n_per_segment):
            t_interp = j / n_per_segment
            k = (1 - t_interp)*k_start + t_interp*k_end
            k_list.append(k)

            if len(k_list) > 1:
                dk = np.linalg.norm(k_list[-1] - k_list[-2])
                k_dist.append(k_dist[-1] + dk)

        node_positions.append(k_dist[-1])

    return np.array(k_list), np.array(k_dist), node_positions

k_list, k_dist, node_positions = interpolate_kpath(k_path)


# BAND TRACKING + SPIN

Sz = kron(I2, s_z)

num_k = len(k_list)
num_bands = 4

bands = np.zeros((num_k, num_bands))
spin_vals = np.zeros((num_k, num_bands))

prev_vecs = None

for i, (kx, ky) in enumerate(k_list):
    H = H_k(kx, ky)
    eigvals, eigvecs = np.linalg.eigh(H)

    if i == 0:
        idx = np.argsort(eigvals)
    else:
        overlap = np.abs(prev_vecs.conj().T @ eigvecs)

        idx = np.zeros(num_bands, dtype=int)
        used = set()

        for n in range(num_bands):
            best = np.argmax(overlap[n])
            while best in used:
                overlap[n, best] = -1
                best = np.argmax(overlap[n])
            idx[n] = best
            used.add(best)

    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    bands[i] = eigvals

    for n in range(num_bands):
        psi = eigvecs[:, n]
        spin_vals[i, n] = np.real(np.vdot(psi, Sz @ psi))

    prev_vecs = eigvecs.copy()


# PLOT (FIXED COLORBAR)

fig, ax = plt.subplots(figsize=(8,6))

for n in range(num_bands):
    for i in range(len(k_dist)-1):
        c = spin_vals[i, n]
        ax.plot(k_dist[i:i+2], bands[i:i+2, n],
                color=plt.cm.coolwarm((c+1)/2),
                lw=2)

# symmetry lines
for pos in node_positions:
    ax.axvline(pos, color='k', linestyle='--', linewidth=0.6)

ax.set_xticks(node_positions)
ax.set_xticklabels(labels)
ax.set_ylabel("Energy")
ax.set_title(f"Kane–Mele Band Structure (⟨s_z⟩ colored)(M={lambda_v },t2={lambda_SO},trashba={lambda_R})")

# colorbar (FIXED)
sm = plt.cm.ScalarMappable(cmap='coolwarm',
                           norm=plt.Normalize(vmin=-1, vmax=1))
fig.colorbar(sm, ax=ax, label=r'$\langle s_z \rangle$')

ax.grid(alpha=0.3)
plt.tight_layout()
plt.show()
