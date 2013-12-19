import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


As = np.arange(0.0, 1.0, 0.05)
Bs = np.arange(0.0, 1.0, 0.05)

splarr = np.loadtxt("test_dat/edge_splarr.dat")

A, B = np.meshgrid(As, Bs)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(A, B, splarr)

ax.set_xlabel("a")
ax.set_ylabel("b")

fig.savefig("test_figs/edge_splarr.pdf")
