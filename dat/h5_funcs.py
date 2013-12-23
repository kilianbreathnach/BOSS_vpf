import numpy as np
import h5py as h5


def mk_h5(infile, outfile, dset_name, **kwargs):
    """
    This function takes an ASCII file with space-separated columns of data
    and makes a HDF5 file of the desired columns of data, skipping an
    optional number of rows and assigning a dataset name.
    """
    rows = kwargs.pop("skiprows", 0)
    cols = kwargs.pop("usecols", None)

    dat = np.loadtxt(infile, skiprows=rows, usecols=cols)

    f = h5.File(outfile, 'a')

    dset = f.create_dataset(dset_name,
                            shape=dat.shape,
                            dtype=dat[0, 0].dtype)

    dset[:, :] = dat[:, :]

    f.close()


def arr2h5(arr, outfile, dset_name, mode='a'):
    """
    This function takes a numpy array and saves it to a HDF5 file.
    """
    f = h5.File(outfile, mode)

    dset = f.create_dataset(dset_name,
                            shape=arr.shape,
                            dtype=arr[0, 0].dtype)

    dset[:, :] = arr[:, :]

    f.close()


def h5_arr(h5file, dset):
    """
    This function copies the desired array from a hdf5 file and returns
    it.
    """
    f = h5.File(h5file, 'r')

    arr = f[dset].value

    f.close()

    return arr
