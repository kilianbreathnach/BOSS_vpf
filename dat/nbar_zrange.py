import os
import json
import numpy as np
from glob import glob
from astropy.io import fits
from astropy.coordinates import Distance
from astropy.cosmology import Planck13, WMAP5
from h5_funcs import arr2h5, h5_arr


"""
This file contains functions to take the BOSS data file for corrected number
density as a function of redshift and process it for our purposes.

mk_zfile decides on a width around the redshift of maximum number density
of each survey that we focus on for our analysis (decided by Q below) and
computes the overall average in that redshift range. It stores
information to The Another function uses the range and
"""


def calc_shell_vol(dis_func, z_near, z_far, zcen):
    """
    Computes the volume of a spherical shell of at redshift zcen, with edges
    at redshifts z_near and z_far, given a cosmological distance conversion
    function. Results are in (h^-1 Mpc)^3.
    """

    r = dis_func(zcen).value
    r_near = dis_func(z_near).value
    r_far = dis_func(z_far).value
    dr = r_far - r_near

    return 4 * np.pi * (r ** 2) * dr


def cut_zs(nz_dict, nbar, radeczfile, cosmo):


    # open files
    radecz = h5_arr(radeczfile, "radecz")
    nbar_arr = np.loadtxt(nbar)

    # Make the redshift cut in the nbar array
    nbar_arr = nbar_arr[(nz_dict["zlo"] < nbar_arr[:, 0]) * \
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
        zbin_ids = np.where(((zbinedges[i] < radecz[:, 2]) * \
                             (radecz[:, 2] <= zbinedges[i + 1])) == True)

        keep = np.random.choice(zbin_ids[0], size=nd, replace=False)

        finmask[keep] = True


    radecz = radecz[finmask]

    arr2h5(radecz, "{0}/{1}/down_radecz.hdf5".format(os.path.dirname(radeczfile), cosmo), "radecz")


def mk_zfile(nbarfile, radeczfile, cosmology, mk_cut=True):

    # magic number for width around maximum
    Q = 0.65
    # magic number for
    Nfrac = (6769.0358 * np.pi) / 129600

    if cosmology == "Planck":
        Planck13.__init__(100.0, Planck13.Om0)
        cosmo = Planck13
    elif cosmology == "WMAP":
        WMAP5.__init__(100.0, WMAP5.Om0)
        cosmo = WMAP5
    comv = cosmo.comoving_distance

    nbar_arr = np.loadtxt(nbarfile)
    nz_dict = {"tophat height for zrange": Q}

    # Cut out the first bit of crap (works for CMASS, dunno about LOWZ)
    ind03 = np.abs(nbar_arr[:, 0] - 0.3).argmin()

    nbar_arr = nbar_arr[ind03:, :]

    zcen = nbar_arr[:, 0]
    z_near = nbar_arr[:, 1]
    z_far = nbar_arr[:, 2]
    gal_counts = nbar_arr[:, 6]

    nbar = []
    shell_vols = []

    for i in range(len(zcen)):

        shell_vols.append(Nfrac * calc_shell_vol(comv, z_near[i], z_far[i], zcen[i]))
        nbar.append(gal_counts[i] / shell_vols[i])

    nbar = np.array(nbar)

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
    nz_dict["shell_vol"] = np.sum(shell_vols[L:R + 1])

    if mk_cut:
        cut_zs(nz_dict, nbarfile, radeczfile, cosmology)

    nf = open("{0}/{1}/nbar_zrange.json".
                  format(os.path.dirname(radeczfile), cosmology), 'w')

    json.dump(nz_dict, nf, sort_keys=True, indent=4, separators=(',', ':\t'))

    nf.close()


if __name__ == "__main__":

    import sys

    if len(sys.argv) != 3:
        print "usage: python nbar_zrange.py\
               <nbarfile + associated fits: '[CMASS, LOWZ]/[NGC, SGC]'>\
               <WMAP/Planck>"
        assert(False)

    nbar_fn = glob("{0}/nbar*dat".format(sys.argv[-1]))[0]
    fits_fn = glob("{0}/*fits".format(sys.argv[-1]))[0]

    mk_zfile(nbar_fn, fits_fn, cosmology)
