import numpy as np
from astropy.cosmology import Planck13, WMAP5
from scipy.interpolate import interp2d
from scipy.spatial import cKDTree


def spherical_cap(h):

    return np.pi * (h ** 2) * (1 - h / 3)


def box_completeness(Nsph, simul_cosmo):

    if simul_cosmo == "Planck":
        # First make h free
        Planck13.__init__(100.0, Planck13.Om0)
        cosmo = Planck13
    elif simul_cosmo == "WMAP":
        WMAP5.__init__(100.0, WMAP5.Om0)
        cosmo = WMAP5

    rad = np.arange(1.0, 67.0, 5.0)

    As = np.arange(0.0, 1.0, 0.05)
    Bs = np.arange(0.0, 1.0, 0.05)

    splarr = np.loadtxt("test_dat/edge_splarr.dat")
    A, B = np.meshgrid(As, Bs)

    inty = interp2d(A[0, :], B[:, 0], splarr)

    bad_pts = 1000 * np.random.rand(138621, 2)
    bad_r = 0.21430156766000885  # determined from quick calculation for now
    bad_A = np.pi * bad_r ** 2

    badbaum = cKDTree(bad_pts)

    for r_i, r in enumerate(rad):

        spheres = 1000 * np.random.rand(Nsph, 2)

        bound_bool = (spheres[:, 0] < r) + (spheres[:, 1] < r) + \
                     ((1000 - spheres[:, 0]) < r) + \
                     ((1000 - spheres[:, 1]) < r)
        bad_inds = np.where(bound_bool == True)
        badsphs = spheres[bound_bool]

        pickle_bool = ((badsphs[:, 0] ** 2 + badsphs[:, 1] ** 2) < r) + \
               (((1000 - badsphs[:, 0]) ** 2 + (1000 - badsphs[:, 0]) ** 2)
                   < r)
        pickle_inds = bad_inds[pickle_bool]

        for i, sph in enumerate(spheres):

            badvol = 0.

            pierce_pts = badbaum.query_ball_point(sph, r)

            for pt in pierce_pts:
                # retrieve coordinates of points within sphere
                pt_coord = bad_pts[pt]
                # calculate fractional projected distance from centre
                dis = np.sqrt((sph[0] - pt_coord[0]) ** 2 + \
                              (sph[1] - pt_coord[1]) ** 2) / r
                # calculate length pierced through sphere
                l = 2 * np.sqrt(1 - dis ** 2)

                bad_vol += l * bad_A

            # check if sphere at boundary
            if i in bad_inds:

                if i in pickle_inds:

                    badvol += inty(sph[0] / r, sph[1] / r)

                else:
                    if sph[0] < r:
                        badvol += spherical_cap(1 - sph[0] / r)
                    elif 1000 - sph[0] < r:
                        badvol += spherical_cap(1 - (1000 - sph[0]) / r)

                    if sph[1] < r:
                        badvol += spherical_cap(1 - sph[1] / r)
                    elif 1000 - sph[1] < r:
                        badvol += spherical_cap(1 - (1000 - sph[1]) / r)

            f = open("test_dat/simul_badvol.dat", 'a')
            f.write("{0}\n".format(badvol))
            f.close()


if __name__ == "__main__":

    import sys

    # Get the runtime arguments

    if len(sys.argv) != 3:
        print "usage: python vpf_analysis.py \
               <No. of spheres at each radius> <Planck | WMAP>"

    box_completeness(int(sys.argv[1]), argv[2])
