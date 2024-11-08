import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import eigh

""" Solving Hydrogen atom (with pre-defined ALPHAs) """

# Number of levels to compute
Nlevels = 5  # 1s 2s 2p 3s 3p 3d

# Basis function definition
# alpha = np.array([13.00773, 1.962079, 0.444529, 0.1219492]) # Standard values form Thjissen for ground state 4 basis function
# alpha = np.array([1])
alpha = 0.0001 * ((1e6) ** np.linspace(0, 1, 40)) # Geometric series

""" Solving The eigenvalue problem """

# Create the matrices amul_{ij} = a_i * a_j and aplu_{ij} = a_i + a_j
amul = np.repeat(alpha[np.newaxis, :], len(alpha), axis=0) * alpha[:, np.newaxis]
aplu = np.repeat(alpha[np.newaxis, :], len(alpha), axis=0) + alpha[:, np.newaxis]

# Compute the H and S matrices directly
H = 3*(np.pi**1.5) * amul / np.power(aplu, 2.5) - 2*np.pi / aplu
S = np.power(np.pi / aplu, 1.5)

# Solve the eigenvalue problem
eigenVals, eigenVects = eigh(H, S, subset_by_index=[0, Nlevels-1])

"""## Plot the result"""

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

ax1.set(xlabel=r"$r$", ylabel=r"$r^2\vert \psi\vert^2$")
ax2.set(ylabel=r"$\varepsilon$ [eV]", xticks=[])

plt.show()
