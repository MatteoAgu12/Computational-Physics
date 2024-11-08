import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import eigh

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

""" Solving Helium atom (with pre-defined ALPHAs) """

Nlevels = 4

# Optima 4 gaussian basis
alpha = np.array([0.298073, 1.242567, 5.782948, 38.474970])  # α (suggested α values)

# Automatic selection
# alpha = 0.01 * ((1e6) ** np.linspace(0, 1, 50)) # α runs from 0.01 to 100 (geometric progression)

# Save number of basis
N = len(alpha)

""" Construction of the matrices involved """

# Create the matrices amul_{ij} = a_i * a_j and aplu_{ij} = a_i + a_j
amul = np.array(np.meshgrid(alpha, alpha)).prod(0)
aplu = np.array(np.meshgrid(alpha, alpha)).sum(0)

# Create the matrices aplumul_{ijml} = (a_i + a_m) * (a_j + a_l) and apluplu_{ijml} = a_i + a_j + a_m + a_l
aplumul = np.swapaxes(np.array(np.meshgrid(aplu, aplu)).prod(0).reshape((N, N, N, N)), 1, 2)
apluplu = np.array(np.meshgrid(aplu, aplu)).sum(0).reshape((N, N, N, N))

# Compute the H, S and Q matrices directly
H = 3*(np.pi**1.5) * amul / np.power(aplu, 2.5) - 4*np.pi / aplu
S = np.power(np.pi / aplu, 1.5)
Q = 2*np.pi**2.5/(aplumul* np.sqrt(apluplu))

""" Self consistent loop """

# Create random wavefunction
C = np.random.randn(1, len(alpha))
C /= np.sqrt(C @ S @ C.T)

# Placeholders for energy
E, prevE, count = 0, 1, 0

# Self consistent loop
print("Iter.  %11s  %11s"%("Energy [eV]", "ΔE [eV]"))
while 27.211386245988*np.abs(E - prevE) > 1e-6:
  # Start of a new cycle
  prevE = E

  # Compute the Fock matrix
  F = H + np.einsum('prqs,ir,is->pq', Q, C, C) # Q[p][r][q][s]*C[r]*C[s]

  # Solve the generalized eigenvalue problem
  eigenVals, eigenVects = eigh(F, S, subset_by_index=[0, Nlevels-1])

  # Take the ground state solution
  C = eigenVects[:, 0].real.reshape((1,N))

  # Evaluate the average energy
  E = 2 * C @ H @ C.T + np.einsum('prqs, ip, iq, ir, is', Q, C, C, C, C) # Q[p][r][q][s]*C[p]*C[q]*C[r]*C[s]

  count += 1
  print("%5d  %11.7f  %11.7f"%(count, 27.211386245988*E, 27.211386245988*(E - prevE)))

print("Final results:")
print("Energy = %11.7f eV\t\tExpected = %11.7f eV  [Took %4d iterations]"%(27.211386245988*E, -79.005151, count))

""" Print results """

# Define r domain [0, 20]
r = np.linspace(0, 20, 1000)

# Construct the basis set
basis = np.exp(- r*r * alpha[:, np.newaxis])

# Compute eigenfunctions
eigen = np.sum(eigenVects[:, :, np.newaxis] * np.repeat(basis[:, np.newaxis, :], Nlevels, axis=1), axis=0)

# Real plotting
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(5,3), dpi=200)
plt.subplots_adjust(wspace=0.5)

ax1.plot(np.repeat(r[np.newaxis, :], Nlevels, axis=0).T, (r *r * eigen * eigen).T)
ax2.plot([eigenVals * 27.211386246] * 2)

ax1.set(xlabel=r"$r$ [A]", ylabel=r"$r^2\vert \psi\vert^2$")
ax2.set(ylabel=r"$\varepsilon$ [eV]", xticks=[])

plt.show()
