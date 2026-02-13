import numpy as np
import matplotlib.pyplot as plt

# Load the .dat file
data = np.loadtxt("../data/after.dat", delimiter=',')

# Column 1 = radius, Column 2 = angle (radians)
r = data[:, 0]
theta = data[:, 1]   # assumed 0 → pi/2

# # ---- Extend solution symmetrically ----
# # Mirror across pi/2  → gives 0 → pi
# theta_half = np.concatenate([theta, np.pi - theta[::-1]])
# r_half = np.concatenate([r, r[::-1]])

# # Mirror across pi → gives 0 → 2pi
# theta_full = np.concatenate([theta_half, theta_half + np.pi])
# r_full = np.concatenate([r_half, r_half])

# ---- Minimalistic polar plot ----
fig, ax = plt.subplots(
    subplot_kw={'projection': 'polar'},
    figsize=(6, 6)
)

# ax.plot(theta_full, r_full, color='black', linewidth=2)
# ax.plot(theta_full, r_full, color='black', linewidth=2)
ax.plot(5*theta, r/r, color='darkgrey', linewidth=1)
ax.plot(5*theta, 3*r/r, color='lightgrey', linewidth=0.5)
# ax.plot(5*theta, 10*r/r, color='darkgrey', linewidth=1)
# ax.plot(5*theta, 30*r/r, color='lightgrey', linewidth=0.5)
# ax.plot(5*theta,100*r/r, color='darkgrey', linewidth=1)
ax.plot(theta, r, color='black', linewidth=1)
# ax.plot(theta + np.pi/2., r[::-1], color='black', linewidth=1)
# ax.plot(theta + np.pi, r, color='black', linewidth=1)
# ax.plot(theta + 3.*np.pi/2., r[::-1], color='black', linewidth=1)

# Minimal styling
ax.set_xticks([])
ax.set_yticks([])
ax.grid(False)
ax.spines['polar'].set_visible(False)

plt.tight_layout()
plt.show()
