import numpy as np
from astropy.cosmology import Planck13, WMAP5
from scipy.interpolate import interp2d
from scipy.spatial import cKDTree


simul_cosmo = "WMAP"


def spherical_cap(h):
    return 0.75 * (h ** 2) * (1 - h / 3)


if simul_cosmo == "Planck":
    # First make h free
    Planck13.__init__(100.0, Planck13.Om0)
    cosmo = Planck13
elif simul_cosmo == "WMAP":
    WMAP5.__init__(100.0, WMAP5.Om0)
    cosmo = WMAP5

Nsph = 10000000
rad = np.arange(5.0, 66.0, 5.0)

As = np.arange(0.0, 1.0, 0.05)
Bs = np.arange(0.0, 1.0, 0.05)

splarr = np.loadtxt("test_dat/edge_splarr.dat")
A, B = np.meshgrid(As, Bs)

inty = interp2d(A[0, :], B[:, 0], splarr)

for r_i, r in enumerate(rad):

    print "Starting radius {0} Mpc...".format(r)

    spheres = 1000 * np.random.rand(Nsph, 2)
    badsphs = spheres[(spheres[:, 0] < r) + (spheres[:, 1] < r) + \
            ((1000 - spheres[:, 0]) < r) + ((1000 - spheres[:, 1]) < r)]

    pickle_bool = ((badsphs[:, 0] ** 2 + badsphs[:, 1] ** 2) < r) + \
           (((1000 - badsphs[:, 0]) ** 2 + (1000 - badsphs[:, 1]) ** 2) < r)
    pickles = badsphs[pickle_bool]

    print "there are {0} corner spheres".format(len(pickles))

    caps = badsphs[~pickle_bool]

    for cap in caps:

        if cap[0] < r:
            a = cap[0] / r
        elif 1000 - cap[0] < r:
            a = (1000 - cap[0]) / r
        else:
            a = None

        if cap[1] < r:
            b = cap[1] / r
        elif 1000 - cap[1] < r:
            b = (1000 - cap[1]) / r
        else:
            b = None

        if a and b:
            badvol = spherical_cap(1 - a) + spherical_cap(1 - b)
        elif a:
                badvol = spherical_cap(1 - a)
        elif b:
                badvol = spherical_cap(1 - b)

        f = open("test_dat/mock_completeness/simul_edge_badvol_{0}.dat".format(r), 'a')
        f.write("{0}\n".format(badvol))
        f.close()

    for pickle in pickles:
        pass

        badvol = inty(pickle[0] / r, pickle[1] / r)[0]

        f = open("test_dat/mock_completeness/simul_edge_badvol_{0}.dat".format(r), 'a')
        f.write("{0}\n".format(badvol))
        f.close()
