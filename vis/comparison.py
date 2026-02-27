import numpy as np
import matplotlib.pyplot as plt

# --- Load data manually ---
data1 = np.loadtxt("../data/relaxed.dat", delimiter=",")
data2 = np.loadtxt("../data/shot.dat", delimiter=",")
data3 = np.loadtxt("../data/spectral.dat", delimiter=",")

# --- Extract columns ---
x1, y1 = data1[:,1], data1[:,0]
x2, y2 = data2[:,1], data2[:,0]
x3, y3 = data3[:,1], data3[:,0]

# --- Plot setup ---
plt.figure(figsize=(8,6))
plt.plot(x1, y1, linestyle='-', label="Relaxation")
plt.plot(x2, y2, linestyle='-', label="Shooting")
# plt.plot(x3, y3, linestyle='-', label="Spectral")


# Vertical line at pi/2
plt.axvline(x=np.pi/2, color='k', linestyle='--', linewidth=1.5, label=r'$\pi/2$')

plt.xlabel(r"$\theta$")
plt.ylabel(r"Horizon radius $r_H$")
plt.title(r"Horizon radius $r_H(\theta)$")
plt.legend()
plt.grid(True)
plt.tight_layout()

# --- Show or save ---
plt.show()
# plt.savefig("plot.png", dpi=300)  # optional
