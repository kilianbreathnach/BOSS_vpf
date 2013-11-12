import subprocess
import h5py as h5
from glob import glob
from astropy.io import fits
from dat.cartesian_cosmo import mk_coords


"""
This script runs the full experiment as outlined in the paper:

    <Name of future paper here and arxiv extension>

The data is downloaded and processed in the dat directory. Then the VPF is
calculated for each of the north and south regions of the survey and all the
relevant plots are made and saved in the plots directory.
"""

####################
# Get all the data #
####################

print "Downloading data..."

tinker_site = "http://cosmo.nyu.edu/~tinker/"
sdss_site = "http://data.sdss3.org/sas/bosswork/boss/lss/"
uname_flag = "--user=sdss3"
pw_flag = "--password=4-surveys"

# NGC data

nbar_N = "nbar-cmass-dr11v1-N-Anderson.dat"
cmass_rands_N = "randoms_1E6_cmass_N_dr11v1.dat"
fits_N = "cmass-dr11v1-N-Anderson.dat.fits"

# get the number densities
subprocess.call(['wget', sdss_site + nbar_N, uname_flag, pw_flag])
# move them to the NGC directory
subprocess.call(["mv", nbar_N, "./dat/NGC/"])

# get the random search points and move them
subprocess.call(['wget', tinker_site + cmass_rands_N])
subprocess.call(["mv", cmass_rands_N, "./dat/NGC/"])

# get the ra, dec and redshift information from sdss fits file
subprocess.call(['wget', sdss_site + fits_N, uname_flag, pw_flag])

fN = fits.open(fits_N)
fN_h5 = h5.File("./dat/NGC/cmass-N-radecz.hdf5", "a")

radecz_N = fN_h5.create_dataset("radecz", shape=(fN[1].data.shape[0], 3),
                            dtype=fN[1].data["RA"][0].dtype)
radecz_N[:, 0] = fN[1].data["RA"]
radecz_N[:, 1] = fN[1].data["DEC"]
radecz_N[:, 2] = fN[1].data["Z"]

fN_h5.close()
subprocess.call(["mv", fits_N, "./dat/NGC/"])

# # SGC data (essentially repeat above)
#
# nbar_S = "nbar-cmass-dr11v1-S-Anderson.dat"
# cmass_rands_S = "randoms_1E6_cmass_S_dr11v1.dat"
# fits_S = "cmass-dr11v1-S-Anderson.dat.fits"
#
# subprocess.call(['wget', sdss_site + nbar_S, uname_flag, pw_flag])
# subprocess.call(["mv", nbar_S, "./dat/SGC/"])
#
# subprocess.call(['wget', tinker_site + cmass_rands_S])
# subprocess.call(["mv", cmass_rands_S, "./dat/SGC/"])
#
# subprocess.call(['wget', sdss_site + fits_S, uname_flag, pw_flag])
#
# fS = fits.open(fits_S)
# fS_h5 = h5.File("./dat/SGC/cmass-S-radecz.hdf5", "a")
#
# radecz_S = fS_h5.create_dataset("radecz", shape=(fS[1].data.shape[0], 3),
#                             dtype=fS[1].data["RA"][0].dtype)
# radecz_S[:, 0] = fS[1].data["RA"]
# radecz_S[:, 1] = fS[1].data["DEC"]
# radecz_S[:, 2] = fS[1].data["Z"]
#
# fS_h5.close()
# subprocess.call(["mv", fits_S, "./dat/SGC/"])


# Convert (RA, Dec, z) to cartesian points for the VPF calculation

radec_fns = glob("")

