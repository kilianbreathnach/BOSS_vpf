import numpy as np
import h5py as h5
from astropy import units as u
from astropy.coordinates import Angle, Distance, ICRSCoordinates
from astropy.cosmology import Planck13, WMAP5
from scipy.spatial import cKDTree


Nsph = 1000000
norm = 1. / Nsph


f_good = h5.File("../dat/NGC/CMASS/srch_pts.hdf5")
good_pts = f_good["good_pts"]

f_bad = h5.File("../dat/NGC/CMASS/veto.hdf5")
bad_pts = f_bad["bad_pts"]

nbar_vals = np.loadtxt("../dat/NGC/CMASS/nbar_zrange.dat")
zlo = nbar_vals[1]
zhi = nbar_vals[2]

bad_r = np.arccos(1.0 - (np.pi * 9.8544099e-05) / (2 * 180 ** 2))
bad_r_deg = np.rad2deg(bad_r)

WMAP5.__init__(100.0, WMAP5.Om0)
cosmo = WMAP5
comv = cosmo.comoving_distance


def radec2xyz(radecarr):

    radecarr = np.atleast_2d(radecarr)
    xyzarr = np.zeros((radecarr.shape[0], 3))
    xyzarr[:, 0] = np.cos(np.radians(radecarr[:, 1])) * \
                   np.cos(np.radians(radecarr[:, 0]))
    xyzarr[:, 1] = np.cos(np.radians(radecarr[:, 1])) * \
                   np.sin(np.radians(radecarr[:, 0]))
    xyzarr[:, 2] = np.sin(np.radians(radecarr[:, 1]))

    return xyzarr


def central_angle(coord1, coord2):

    f1 = np.radians(coord1[1])
    f2 = np.radians(coord2[1])

    dl = abs(coord1[0] - coord2[0])

    if dl > 180:
        dl = np.radians(dl - 180)
    else:
        dl = np.radians(dl)
    # angle is (from wikipedia formula...)
    dsig = np.arccos(np.sin(f1) * np.sin(f2) + \
                     np.cos(f1) * np.cos(f2) * np.cos(dl))

    return dsig


bad_xyz = radec2xyz(bad_pts)
veto_baum = cKDTree(bad_xyz)

rad = np.arange(1.0, 62.0, 5.0)

rand_i = 0

for r_i, r in enumerate(rad):

    # start the count of successful voids
    count = 0

    for i in range(Nsph):

        rand_i = rand_i % 999999  # compensate for finite length of mask file

        radec = good_pts[rand_i, :]

        rang = Angle(radec[0], u.deg)
        decang = Angle(radec[1], u.deg)

        randz = (zlo ** 3 + \
                 (zhi ** 3 - zlo ** 3) * np.random.rand(1)[0]) ** (1. / 3.)
        dis = Distance(comv(randz), u.Mpc)

        coord = ICRSCoordinates(rang, decang, distance=dis)

        sph_cen = np.array([coord.x.value, coord.y.value, coord.z.value])

        print "rad: ", r, ", sphere: ", i

        # Get radius of circular projection of sphere
        R = np.arcsin(r / np.sqrt(np.sum(sph_cen[:] ** 2)))

        # Get coordinates of circle centre on unit sphere
        crc_cen = radec2xyz(radec)[0]

        # Compute tree search radius from Cosine rule
        # (include points extending beyond sphere edge to account for
        # finite area around bad points)
        l_srch = np.sqrt(2. - 2. * np.cos(R))

        # Run search
        pierce_l = veto_baum.query_ball_point(crc_cen, l_srch)

        bad_vol = 0.

        R = np.degrees(R)  # need in degrees for bad_vol computation

        for pt in pierce_l:

            pt_ang = bad_pts[pt]
            dis = np.degrees(central_angle(pt_ang, radec))
            l = dis / R

            bad_vol += 1.5 * (bad_r_deg / R) ** 2 * np.sqrt(1.0 - l ** 2)

        f_r = open("./test_dat/NGC/CMASS/completeness_1mil/volfrac_rad{0}.dat".
                       format(r), 'a')
        f_r.write("{0}\n".format(bad_vol))
        f_r.close()

        rand_i += 1
