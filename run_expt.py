import os
import subprocess
import requests
from glob import glob
import numpy as np
from dat.cartesian_cosmo import mk_coords
from dat.h5_funcs import mk_h5, fits2h5
from dat.nbar_zrange import process_nbar


"""
This script runs the full experiment as outlined in the paper:

    <Name of future paper here and arxiv URL>

The data is downloaded and processed in the dat directory. Then the VPF is
calculated for each of the north and south regions of the survey and all the
relevant plots are made and saved in the plots directory.
"""

while True:
    get_dat = raw_input("Do you need to download the input data?[y/n]")
    if get_dat != 'y' and get_dat != 'n':
        print "Not a valid response"
    else:
        break

if get_dat == 'y':

    ####################
    # Get all the data #
    ####################

    print "\nDownloading data...\n\n"

    # Location of the data on the interwebz
    broiler_dat = "http://broiler.astrometry.net/~kilian/tinker/proj_dat/"
    sdss_site = "http://data.sdss3.org/sas/bosswork/boss/lss/"

    # Prompts for the sdss3 data username and password,
    # please contact the authors for information on these
    username = raw_input("please enter the username for the BOSS\
                                        data page: ")
    password = raw_input("and the corresponding password: ")


    # CMASS data
    # - NGC data
    ############

    # File names of all the input data
    nbar = "nbar-cmass-dr11v1-N-Anderson.dat"
    fitsfile = "cmass-dr11v1-N-Anderson.dat.fits"
    rands = "randoms_1E6_cmass_N_dr11v1.dat"

    # Create necessary directory structure if not in place
    sdss_dir = "./dat/in/sdss3"
    if not os.path.exists(sdss_dir):
        os.makedirs(sdss_dir)

    cmass_dir = "./dat/in/CMASS_DATA"
    if not os.path.exists(cmass_dir):
        os.makedirs(cmass_dir)

    # First, get the corrected number densities
    r_nbar = requests.get(sdss_site + nbar, auth=(username, password))

    # move them to the sdss3 data directory
    f = open("{0}/{1}".format(sdss_dir, nbar), 'w')
    f.write(r_nbar.text)
    f.close()

    # get the ra, dec and redshift information from sdss fits file
    r_radecz = requests.get(sdss_site + fitsfile, auth=(username, password))

    f = open("{0}/{1}".format(sdss_dir, fitsfile), 'w')
    f.write(r_fitsfile.text)
    f.close()

    # get the random search points and move them to CMASS_DATA
    r_rands = requests.get(broiler_dat + rands)

    f = open("{0}/{1}".format(cmass_dir, rands), 'w')
    f.write(r_rands.text)
    f.close()

    # get the bad points
    r_bad1 = requests.get(broiler_dat + "north_block_badfield.dat")
    r_bad2 = requests.get(broiler_dat + "north_block_brightstar.dat")
    r_bad3 = requests.get(broiler_dat + "north_block_outside.dat")

    f = open("{0}/{1}".format(cmass_dir, "north_block_badfield.dat"), 'w')
    f.write(r_bad1.text)
    f.close()

    f = open("{0}/{1}".format(cmass_dir, "north_block_brightstar.dat"), 'w')
    f.write(r_bad2.text)
    f.close()

    f = open("{0}/{1}".format(cmass_dir, "north_block_outside.dat"), 'w')
    f.write(r_bad3.text)
    f.close()


# Now to get nice HDF5 files made of everything of interest in the right place

print "\nProcessing Data...\n\n"

# CMASS
# - NGC

sdss_dir = "./dat/in/sdss3"
CMASS_dir = "./dat/in/CMASS_DATA"

out_dir = "./dat/out/CMASS/NGC"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# Get RA, Dec, and redshifts from the fits file into a HDF5 file
fits2h5("{0}/cmass-dr11v1-N-Anderson.dat.fits".format(sdss_dir),
        1, ["RA", "DEC", "Z"], "{0}/full_radecz.hdf5".format(out_dir), "radecz")

# create HDF5 of the random search RA, Dec values
mk_h5("{0}/randoms_1E6_cmass_N_dr11v1.dat".format(CMASS_dir),
      "{0}/srch_radec.hdf5".format(out_dir), "good_pts", skiprows=1, usecols=(0, 1))

# Make a unified hdf5 veto file out of the bad points files
catfiles = glob("{0}/north_block_*".format(CMASS_dir))

with open("{0}/full_veto.dat".format(out_dir), 'w') as outfile:
    for fname in catfiles:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

mk_h5("{0}/full_veto.dat".format(out_dir), "{0}/veto.hdf5".format(out_dir),
      "bad_pts", usecols=(0, 1))


# Generate the data for specific cosmologies
# - WMAP

cosmo = "WMAP"

cosmo_dir = "{0}/{1}".format(out_dir, cosmo)
if not os.path.exists(cosmo_dir):
    os.makedirs(cosmo_dir)

# get the zrange around the mean nbar
process_nbar("{0}/nbar-cmass-dr11v1-N-Anderson.dat".format(sdss_dir),
             "{0}/{1}/nbar_zrange.json".format(out_dir, cosmo), cosmo,
             radeczfile="{0}/full_radecz.hdf5".format(out_dir))

# Convert (RA, Dec, z) to cartesian points for the VPF calculation

radec_fns = glob("{0}/{1}/radecz_down.hdf5".format(out_dir, cosmo))

mk_coords(radec_fns[0],
          "{0}/{1}/gals_cart_coords.hdf5".format(out_dir, cosmo),
          cosmo)


# Eventually need to include SGC and LOWZ components too...


# Now to run the VPF analysis on our data
