import os
import numpy as np
from glob import glob


def mk_zfile(dset):

    # width around maximum
    Q = 0.85

    nbar_arr = np.loadtxt(dset)

    # Cut out the first bit of crap (works for CMASS, dunno about LOWZ)
    ind03 = np.abs(nbar_arr[:, 0] - 0.3).argmin()

    nbar_arr = nbar_arr[ind03:, :]

    zcen = nbar_arr[:, 0]
    nbar = nbar_arr[:, 3]

    # Find nbar peak and index
    max_nbar = np.max(nbar)
    max_i = int(np.where(nbar == max_nbar)[0])

    # get the interval edge indices for 10%
    L = np.abs(nbar[:max_i] - max_nbar * Q).argmin()
    R = max_i + np.abs(nbar[max_i:] - max_nbar * Q).argmin()

    zL = zcen[L]
    zR = zcen[R]

    avg_nbar = np.average(nbar[L:R + 1])

    nf = open("{0}/nbar_zrange.dat".format(os.path.dirname(dset)), 'w')
    nf.write("# average nbar around peak, nearest redshift of interval, \
              furthest redshift\n")
    nf.write("{0}\t{1}\t{2}".format(avg_nbar, zL, zR))
    nf.close()


if __name__ == "__main__":

    import sys

    if len(sys.argv) != 2:
        print "usage: python nbar_zrange.py\
               <nbarfile: '[NGC, SGC]/[CMASS, LOWZ]'>"
        assert(False)

    nbar_fn = glob("{0}/nbar*".format(sys.argv[-1]))[0]

    mk_zfile(nbar_fn)
