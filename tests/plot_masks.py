import sys
import numpy as np
import h5py as h5
from gen_plot import plotz


"""
This script plots the points in the randoms mask and in the veto mask
together in RA and Dec for the survey and cap given in the call
arguments.
"""


if len(sys.argv) != 2:

    print "usage: python plot_masks.py <[CMASS, LOWZ]/[NGC, SGC]>"

dirname = "../dat/" + sys.argv[1]

fgood = h5.File(dirname + "/srch_pts.hdf5")
fbad = h5.File(dirname + "/veto.hdf5")

good_pts = fgood["good_pts"]
bad_pts = fbad["bad_pts"]

plotz([good_pts[0, :], bad_pts[0, :]], [bad_pts[1, :], bad_pts[1, :]],
      style='scatter', overplot=True, colours=["g", "r"], alpha=0.1,
      name="./test_figs/{0}/masks.png".format(dirname))


fgood.close()
fbad.close()
