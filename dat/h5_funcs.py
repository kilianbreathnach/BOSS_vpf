import numpy as np
import h5py as h5


def mk_h5(infile, outfile, dset_name, **kwargs):

    rows = kwargs.pop("skiprows", 0)
    cols = kwargs.pop("usecols", None)

    dat = np.loadtxt(infile, skiprows=rows, usecols=cols)

    f = h5.File(outfile, 'a')

    dset = f.create_dataset(dset_name, shape=dat.shape, dtype=dat[0, 0].dtype)

    dset[:, :] = dat[:, :]

    f.close()
