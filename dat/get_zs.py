import sys
import numpy as np
from astropy.cosmology import Planck13, WMAP5


if len(sys.argv) != 3:
    print "usage: python get_zs.py <survey (=CMASS/LOWZ/theoretical)>\
           <Planck | WMAP>"

# width around maximum
Q = 0.85

hemis = ["NGC", "SGC"]


# Set the cosmology with h free
# Here the cosmology is based on WMAP (for first MultiDark simulation)
if sys.argv[-2] == "Planck":
    # First make h free
    Planck13.__init__(100.0, Planck13.Om0)
    cosmo = Planck13
elif sys.argv[-2] == "WMAP":
    WMAP5.__init__(100.0, WMAP5.Om0)
    cosmo = WMAP5
comv = cosmo.comoving_distance

if sys.argv[1] == "CMASS":
    north = np.loadtxt("./dat/NGC/nbar-cmass-dr11v1-N-Anderson.dat")
    south = np.loadtxt("./dat/SGC/nbar-cmass-dr11v1-S-Anderson.dat")

    # Cut out the first bit of crap
    ind03_N = np.abs(north[:, 0] - 0.3).argmin()
    ind03_S = np.abs(south[:, 0] - 0.3).argmin()

    north = north[ind03_N:, :]
    south = south[ind03_S:, :]

    zcen_N = north[:, 0]
    zcen_S = south[:, 0]
    nbar_N = north[:, 3]
    nbar_S = south[:, 3]

    # add data from both files using above weights
    zcen = zcen_N  # since both files have the same bins

    for i, nbar in enumerate([nbar_N, nbar_S]):

        # Find nbar peak and index
        max_nbar = np.max(nbar)
        max_i = int(np.where(nbar == max_nbar)[0])

        # get the interval edge indices for 10%
        L = np.abs(nbar[:max_i] - max_nbar * Q).argmin()
        R = max_i + np.abs(nbar[max_i:] - max_nbar * Q).argmin()

        zL = zcen[L]
        zR = zcen[R]

        avg_nbar = np.average(nbar[L:R + 1])

        nf = open("./dat/{}/cmass_nbar.dat".format(hemis[i]), 'w')
        nf.write("# average nbar around peak, nearest redshift of interval, \
                  furthest redshift\n")
        nf.write("{0}\t{1}\t{2}".format(avg_nbar, zL, zR))
        nf.close()
