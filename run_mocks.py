from glob import glob
from mocks import mk_mock_coords, mock_vpf
from dat.h5_funcs import mk_h5

survey_cap = "CMASS/NGC"
cosmo = "WMAP"
spheresfile = "dat/out/{0}/{1}/mocks/mock_srch_pts.hdf5".format(survey_cap, cosmo)

# glob the raw mock files
mock_rdzs = glob("dat/mocks/ngc/*.rdz")

# just hte first 5 for now
for i in range(5):

    mock_rdz = mock_rdzs[i]

    # get the mock number
    name = mock_rdz.split(".")[1].split("_")[-1]

    # get the location to put the hdf5 of the mock
    mock_hier_fn = "dat/out/{0}/mocks_hierarchical/{1}.hdf5".format(survey_cap, name)

    mk_h5(mock_rdz, mock_hier_fn, "radecz", mode='w')

    # the location for the transformed coordinates of the mock galaxies
    cart_coords_fn = "dat/out/{0}/{1}/mocks/cart_coords/{2}.hdf5".format(survey_cap, cosmo, name)

    mk_mock_coords(mock_hier_fn, cart_coords_fn, cosmo)

    mock_vpf(cart_coords_fn, spheresfile, cosmo)
