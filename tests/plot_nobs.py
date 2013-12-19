import numpy as np
from astropy.io import fits


nbar_arr = np.loadtxt("../dat/CMASS/NGC/nbar-cmass-dr11v1-N-Anderson.dat")

f = fits.open("../dat/CMASS/NGC/cmass-dr11v1-N-Anderson.dat.fits")

zbinedges = np.append(np.array([0.0]), nbar_arr[:, 2])
shell_vols = nbar_arr[:, -2]

dat_zs = f[1].data["Z"]

H = np.histogram(dat_zs, bins=zbinedges)

nbar_obs = H[0] / shell_vols


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)


fig = plt.figure()
ax = fig.add_subplot(111)

l1, = ax.plot(nbar_arr[:, 0], nbar_arr[:, 3], label=r"complete gals")
l2, = ax.plot(nbar_arr[:, 0], nbar_obs, label=r"observed gals")

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)

ax.set_title(r"Comparing observed number density to full n")
ax.set_xlabel(r"$z$")
ax.set_ylabel(r"$\bar{n}$ $(\frac{h^{3}}{Mpc^{3}})$")

plt.xlim(0.02, 1.1)
plt.ylim(0., 0.0005)

fig.savefig("/home/kilian/public_html/tinker/nbar_obs_compare.png")
