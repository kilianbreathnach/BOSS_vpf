import sys
from glob import glob
from mocks_ascii import mk_mock_coords, mock_vpf
from dat.h5_funcs import mk_h5

survey_cap = "CMASS/NGC"
cosmo = "WMAP"
spheresfile = "dat/out/{0}/{1}/mocks/mock_srch_pts.hdf5".format(survey_cap, cosmo)

if len(sys.argv) != 2:
    print "usage: python run_mocks.py <mock no.: int>"

digits = len(sys.argv[1])

i = 0
name = ""
while i < (4 - digits):
    name += "0"

name += sys.argv[1]

# glob the raw mock files
mock_rdzs = glob("dat/mocks/ngc/*{0}*.rdz".format(name))

for i in range(len(mock_rdzs)):

    mock_rdz = mock_rdzs[i]

    # get the mock number
    name = mock_rdz.split(".")[1].split("_")[-1]

    # get the location to put the hdf5 of the mock
    mock_hier_fn = "dat/out/{0}/mocks_hierarchical/{1}.hdf5".format(survey_cap, name)

    mk_h5(mock_rdz, mock_hier_fn, "radecz", mode='w')

    # the location for the transformed coordinates of the mock galaxies
    cart_coords_fn = "dat/out/{0}/{1}/mocks/cart_coords/ascii/{2}.dat".format(survey_cap, cosmo, name)

    mk_mock_coords(mock_hier_fn, cart_coords_fn, cosmo)

        #    mock_vpf(cart_coords_fn, spheresfile, cosmo)
