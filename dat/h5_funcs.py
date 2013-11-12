import numpy as np
import h5py as h5


def mk_h5(infile, outfile, dset_name, **kwargs):

    cols = kwargs.clearpop("usecols", None)

    dat = np.loadtxt(infile, usecols=cols)

    f = h5.File(outfile, 'a')

    dset = f.create_dataset(dset_name, shape=dat.shape, dtype=dat[0, 0].dtype)

    dset[:, :] = dat[:, :]

    f.close()
