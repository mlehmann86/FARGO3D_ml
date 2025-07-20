import numpy as np
import matplotlib.pyplot as plt
import glob

# Gather all rank files
filelist = sorted(glob.glob('viscosity_profile_rank*.dat'))

# Load and concatenate
r_total = []
nu_total = []

for filename in filelist:
    r, nu = np.loadtxt(filename, unpack=True)
    r_total.append(r)
    nu_total.append(nu)

r_total = np.concatenate(r_total)
nu_total = np.concatenate(nu_total)

# Sort by radius if needed
idx = np.argsort(r_total)
r_total = r_total[idx]
nu_total = nu_total[idx]

# Plot
plt.plot(r_total, nu_total, marker='o', markersize=2, linestyle='-')
plt.xlabel('Radius r')
plt.ylabel('Viscosity ν')
plt.yscale('log')
plt.grid(True, which='both', ls='--')
plt.tight_layout()
plt.show()
