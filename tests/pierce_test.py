import numpy as np


"""
Here we check the error in estimating the compromised volume of a sphere
by a surface that intersects with the 2D projection of the sphere. The volume
is determined by associating a characteristic area with each point in the
compromising surface based on the number density of random points. This
area is multiplied by the length of intersection through the sphere to give an
approximate volume. The error should reduce with increased number density
of points.
"""

Nptsarr = np.logspace(4.0, 7.0, num=7).astype(int)
errs = []

for Npts in Nptsarr:

    # Set up a grid of bad points
    bad_pts = 2 * np.random.rand(Npts, 2) - 1

    # effective area per point
    Apt = 4. / Npts

    # get distance of each point to the origin
    r = np.sqrt(bad_pts[:, 0] ** 2 + bad_pts[:, 1] ** 2)

    # cut out points outside of the circle
    r = r[np.where(r < 1.0)]

    # get array of volumes compromised by points
    fac = 2 * Apt
    vols = fac * np.sqrt(1.0 - r ** 2)

    # print the percentage error in volume of sphere
    volfac = (4. / 3.) * np.pi
    volnorm = 1. / volfac
    errs.append(volnorm * np.sum(vols) - 1.)

errs = np.array(errs)

import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(Nptsarr, errs)

fig.savefig("/home/kilian/public_html/tinker/pierce_test.pdf")
