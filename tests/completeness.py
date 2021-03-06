import time
import json
import numpy as np
import h5py as h5
from astropy import units as u
from astropy.coordinates import Angle, Distance, ICRSCoordinates
from astropy.cosmology import Planck13, WMAP5
from scipy.spatial import cKDTree


Nsph = 10000000
norm = 1. / Nsph


f_good = h5.File("../dat/out/CMASS/NGC/srch_radec.hdf5")
good_pts = f_good["good_pts"]

f_bad = h5.File("../dat/out/CMASS/NGC/veto.hdf5")
bad_pts = f_bad["bad_pts"]

nbar_dict = json.load(open("../dat/out/CMASS/NGC/WMAP/nbar_zrange.json"))
zlo = nbar_dict["zlo"]
zhi = nbar_dict["zhi"]

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

rad = np.arange(5.0, 66.0, 5.0)

rand_i = 0

t_0 = time.time()

for r_i, r in enumerate(rad):

    print "time: {0} seconds".format(time.time() - t_0)
    print " - starting search at radius {0} Mpc".format(r)

    # Custom zrange for sphere size
    dis_near = Distance(comv(zlo).value + r, u.Mpc)
    dis_far = Distance(comv(zhi).value - r, u.Mpc)

    z_a = dis_near.compute_z(cosmology=cosmo)

    z_b = dis_far.compute_z(cosmology=cosmo)

    randz = (z_a ** 3 + \
             (z_b ** 3 - zlo ** 3) * np.random.rand(Nsph)) ** (1. / 3.)

    for i in range(Nsph):

        rand_i = rand_i % 999999  # compensate for finite length of mask file

        radec = good_pts[rand_i, :]

        rang = Angle(radec[0], u.deg)
        decang = Angle(radec[1], u.deg)

        dis = Distance(comv(randz[i]), u.Mpc)

        coord = ICRSCoordinates(rang, decang, distance=dis)

        sph_cen = np.array([coord.x.value, coord.y.value, coord.z.value])

#        print "rad: ", r, ", sphere: ", i

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

        f_r = open("./test_dat/CMASS/NGC/full_completeness_10mil/volfrac_rad{0}.dat".
                       format(r), 'a')
        f_r.write("{0}\n".format(bad_vol))
        f_r.close()

        rand_i += 1
