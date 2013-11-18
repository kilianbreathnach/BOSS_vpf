import os
import subprocess
import h5py as h5
from glob import glob
from astropy.io import fits
from dat.cartesian_cosmo import mk_coords
from dat.h5_funcs import mk_h5
from dat.nbar_zrange import mk_zfile


"""
This script runs the full experiment as outlined in the paper:

    <Name of future paper here and arxiv URL>

The data is downloaded and processed in the dat directory. Then the VPF is
calculated for each of the north and south regions of the survey and all the
relevant plots are made and saved in the plots directory.
"""

####################
# Get all the data #
####################

print "Downloading data..."

broiler_dat = "http://broiler.astrometry.net/~kilian/tinker/proj_dat/"
tinker_site = "http://cosmo.nyu.edu/~tinker/"
sdss_site = "http://data.sdss3.org/sas/bosswork/boss/lss/"

uname_flag = "--user=" + raw_input("please enter the username for the BOSS\
                                    data page: ")
pw_flag = "--password=" + raw_input("and the corresponding password: ")


# CMASS data
# - NGC data

survey_cap = "CMASS/NGC"

nbar = "nbar-cmass-dr11v1-N-Anderson.dat"
fitsfile = "cmass-dr11v1-N-Anderson.dat.fits"
rands = "randoms_1E6_cmass_N_dr11v1.dat"

# get the number densities
subprocess.call(['wget', sdss_site + nbar, uname_flag, pw_flag])
# move them to the NGC directory
subprocess.call(["mv", nbar, "./dat/{0}/".format(survey_cap)])

# get the ra, dec and redshift information from sdss fits file
subprocess.call(['wget', sdss_site + fitsfile, uname_flag, pw_flag])

f = fits.open(fitsfile)
fh5 = h5.File("./dat/{0}/radecz.hdf5".format(survey_cap), "a")

radecz = fh5.create_dataset("radecz", shape=(f[1].data.shape[0], 3),
                            dtype=f[1].data["RA"][0].dtype)
radecz[:, 0] = f[1].data["RA"]
radecz[:, 1] = f[1].data["DEC"]
radecz[:, 2] = f[1].data["Z"]

fh5.close()
subprocess.call(["mv", fitsfile, "./dat/NGC/CMASS/"])


# get the random search points and move them
subprocess.call(['wget', tinker_site + rands])

mk_h5(rands, "./dat/NGC/CMASS/srch_pts.hdf5", "good_pts", usecols=(0, 1))

subprocess.call(["mv", rands, "./dat/NGC/CMASS/"])


# Pull the bad points from broiler
subprocess.call(['wget', "-nd", "-r", "-l", "1", "-P", "./dat/NGC/CMASS/",
                 "-A", ".dat", broiler_dat])
catfiles = glob("./dat/NGC/CMASS/north_block_*")

with open("./dat/NGC/CMASS/full_veto.dat", 'w') as outfile:
    for fname in catfiles:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

# make a unified hdf5 veto file out of them
mk_h5("./dat/NGC/CMASS/full_veto.dat", "./dat/NGC/CMASS/veto.hdf5",
      "bad_pts", usecols=(0, 1))


# get the zrange around the mean nbar

mk_zfile("./dat/NGC/CMASS/nbar-cmass-dr11v1-N-Anderson.dat", fitsfile)


# Convert (RA, Dec, z) to cartesian points for the VPF calculation

radec_fns = glob("./dat/NGC/CMASS/radecz.hdf5")

mk_coords(radec_fns[0], "./dat/NGC/CMASS/WMAP_cart_coords.hdf5", "WMAP")

# Eventually need to include SGC and LOWZ components too...
