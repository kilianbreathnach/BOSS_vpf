import os
import json
import numpy as np
from glob import glob
from astropy.io import fits
from h5_funcs import arr2_h5


"""
This file contains functions to take the BOSS data file for corrected number
density as a function of redshift and process it for our purposes.

mk_zfile decides on a width around the redshift of maximum number density
of each survey that we focus on for our analysis (decided by Q below) and
computes the overall average in that redshift range. It stores
information to The Another function uses the range and
"""


def cut_zs(nz_dict, nbar, fitsfile):

    # open files
    f = fits.open(fitsfile)
    nbar_arr = np.loadtxt(nbar)

    radecz = np.zeros((f[1].data.shape[0], 3))

    # get the right columns
    radecz[:, 0] = f[1].data["RA"]
    radecz[:, 1] = f[1].data["DEC"]
    radecz[:, 2] = f[1].data["Z"]

    # Make the redshift cut in the nbar array
    nbar_arr = nbar_arr[(nz_dict["zlo"] < nbar_arr[:, 0]) & \
                        (nbar_arr[:, 0] < nz_dict["zhi"])]

    # Get binning those observed galaxies
    zbinedges = np.append(nbar_arr[0, 1], nbar_arr[:, 2])

    # Find the counts per bin
    H = np.histogram(radecz[:, 2], bins=zbinedges)

    # The number to downsample to in each bin
    # (multiply bin number by the relative fraction determined from
    #  corrected distribution of nbar)
    num_down = (nz_dict["nbar_corr_tophat"] / nbar_arr[:, 3]) * H[0]

    # make a mask for the final array for analysis within the redshift limits
    finmask = np.array(radecz.shape[0] * [False])

    for i, nd in enumerate(num_down):
        """Turn on the right amount of galaxies in each bin."""
        zbin_ids = np.where( ( (zbinedges[i] < radecz[:, 2]) & \
                               (radecz[:, 2] <= zbinedges[i + 1]) ) == True )

        keep = np.random.choice(zbin_ids[0], size=nd, replace=False)

        finmask[keep] = True

    radecz = radecz[(finmask,)]

    arr2_h5(radecz,
            "./{0}/radecz.hdf5".format(os.path.dirname(nbar)),
            "radecz")


def mk_zfile(nbarfile, fitsfile, mk_cut=True):

    # width around maximum
    Q = 0.85

    nbar_arr = np.loadtxt(nbarfile)
    nz_dict = {"Q": Q}

    # Cut out the first bit of crap (works for CMASS, dunno about LOWZ)
    ind03 = np.abs(nbar_arr[:, 0] - 0.3).argmin()

    nbar_arr = nbar_arr[ind03:, :]

    zcen = nbar_arr[:, 0]
    nbar = nbar_arr[:, 3]

    # Find nbar peak and index
    max_nbar = np.max(nbar)
    max_i = int(np.where(nbar == max_nbar)[0])

    nz_dict["max_nbar"] = max_nbar
    nz_dict["nbar_corr_tophat"] = Q * max_nbar
    nz_dict["z_nbar_max"] = zcen[max_i]

    # get the interval edge indices for 10%
    L = np.abs(nbar[:max_i] - max_nbar * Q).argmin()
    R = max_i + np.abs(nbar[max_i:] - max_nbar * Q).argmin()

    nz_dict["zlo"] = zcen[L]
    nz_dict["zhi"] = zcen[R]

    nz_dict["avg_nbar_corr"] = np.average(nbar[L:R + 1])
    nz_dict["shell_vol"] = np.sum(nbar_arr[L:R + 1, -2])

    if mk_cut:
        cut_zs(nz_dict, nbarfile, fitsfile)

    nf = open("{0}/nbar_zrange.json".format(os.path.dirname(nbarfile)), 'w')

    json.dump(nz_dict, nf, sort_keys=True, indent=4, separators=(',', ':\t'))

    nf.close()


if __name__ == "__main__":

    import sys

    if len(sys.argv) != 3:
        print "usage: python nbar_zrange.py\
               <nbarfile + associated fits: '[NGC, SGC]/[CMASS, LOWZ]'>"
        assert(False)

    nbar_fn = glob("{0}/nbar*dat".format(sys.argv[-1]))[0]
    fits_fn = glob("{0}/*fits".format(sys.argv[-1]))[0]

    mk_zfile(nbar_fn, fits_fn)
