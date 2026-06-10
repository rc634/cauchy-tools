import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../data/chased_surf.dat", comments="#")

theta = data[:, 0]
f     = data[:, 1]
psi   = data[:, 2]
divS  = data[:, 3]
Ricci = data[:, 4]

trims = [0,2,4,6,8,20,80,160,320]
ti = 1
trim = trims[ti]

N = len(theta)
means = []
for k in trims:
    sl     = slice(k, N - k) if k > 0 else slice(None)
    w      = np.sin(theta[sl])
    mean_R = np.sum(Ricci[sl] * w) / np.sum(w)
    means.append(mean_R)
    print(f"trim={k:2d}  <Ricci> = {mean_R:.6e}")


hw = trim
Ricci_smooth = np.empty(N)
for i in range(hw, N - hw):
    Ricci_smooth[i] = np.mean(Ricci[i - hw : i + hw + 1])
slope_l = Ricci_smooth[2*hw] - Ricci_smooth[2*hw - 1]
for i in range(2*hw - 1, -1, -1):
    Ricci_smooth[i] = Ricci_smooth[i + 1] - slope_l / 4
slope_r = Ricci_smooth[N - 2*hw] - Ricci_smooth[N - 2*hw - 1]
for i in range(N - 2*hw, N):
    Ricci_smooth[i] = Ricci_smooth[i - 1] + slope_r / 4

w_all = np.sin(theta)
area  = np.sum(w_all)
K2    = np.sum(w_all * (1. - Ricci_smooth/means[ti])**2) / area
print(f"root-K2 = {np.sqrt(K2):.6e}")

bg = '#0d1117'
fig, ax = plt.subplots(facecolor=bg)
ax.set_facecolor(bg)
ax.plot(theta, Ricci, linewidth=1.0, color='cyan')
ax.plot(theta, Ricci_smooth, color='r')
ax.set_xlabel(r"$\theta$", color='white')
ax.set_ylabel(r"$R$", color='white')
ax.set_title("Horizon shape", color='white')
ax.tick_params(colors='white')
for spine in ax.spines.values():
    spine.set_edgecolor('#444444')
colors = ['#5bc8ef', '#ffaa44', '#57d68a', '#ff7eb3', '#7eb8f7', '#f0e040', '#b57bee', '#ff8080', '#4ecba8', '#e09070']
for m, k, c in zip(means, trims, colors):
    ax.axhline(m, linewidth=0.7, linestyle='--', color=c, label=f"trim={k}")

rs_min, rs_max = Ricci_smooth.min(), Ricci_smooth.max()
margin = 0.1 * (rs_max - rs_min)
ax.set_ylim(rs_min - margin, rs_max + margin)
ax.set_xlim(theta[0], theta[-1])

ax.text(0.98, 0.12, rf"$K_2 = {K2:.4e}$",
        transform=ax.transAxes, ha='right', va='bottom',
        color='#ff4444', fontsize=13, fontweight='bold')
ax.text(0.98, 0.04, rf"$\langle R \rangle = {means[ti]:.4e}$",
        transform=ax.transAxes, ha='right', va='bottom',
        color='#ff4444', fontsize=13, fontweight='bold')
ax.legend(facecolor='#1a1f2e', edgecolor='#444444', labelcolor='white')
plt.tight_layout()
plt.show()
