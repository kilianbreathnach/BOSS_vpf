import numpy as np


As = np.arange(0.0, 1.0, 0.05)
Bs = np.arange(0.0, 1.0, 0.05)

splarr = np.zeros((As.shape[0], Bs.shape[0]))

Nrand = 100000
norm = 1. / Nrand

for i, a in enumerate(As):

    for j, b in enumerate(Bs):

        cube = 2 * np.random.rand(Nrand, 3) - 1

        sph = cube[(cube[:, 0] ** 2 + cube[:, 1] ** 2 + cube[:, 2] ** 2) < 1]

        badarr = np.zeros(sph.shape[0])
        badarr[:] = (sph[:, 0] > a) * (sph[:, 1] > b)
        Nbad = np.sum(badarr)

        print Nbad

        splarr[i, j] = Nbad * norm

np.savetxt("test_dat/edge_splarr.dat", splarr)
