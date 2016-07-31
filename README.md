# Shift
Yet another wrapper around SHIFTX+ for MD trajectories analysis but using MPI and HDF5

## Prerequisites

You need, at a minimum:

* Python 2.7
* NumPy
* Pandas
* MDAnalysis
* H5py
* Matplotlib
* MPI4py

And *SHIFTX2* chemical shift prediction tool (http://www.shiftx2.ca/).

## Installation

I highly recommand you to install the Anaconda distribution (https://www.continuum.io/downloads) if you want a clean python environnment with nearly all the prerequisites already installed (NumPy, Pandas, Matplotlib, MPI).

For the rest (nearly), you just have to do this.
```bash
conda install mpi4py
pip install MDAnalysis
```

In order to use hdf5 format with MPI, you need a special version compiled in parallel mode (with mpi4py).

1 . First, remove hdf5 and h5py from Anaconda
```bash
conda remove hdf5 h5py
```

2 . After, download the latest version of hdf5 and h5py
```bash
wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.17.tar
git clone https://github.com/h5py/h5py.git
```

3 . Compilation of HDF5
```bash
tar -xvf hdf5-1.8.17.tar
cd hdf5-1.8.17
PREFIX=~/Applications/anaconda    # You need to adapt this to your system
./configure --prefix=$PREFIX --enable-linus-lfs --with-zlib=$PREFIX --enable-parallel --enable-shared
make
make install
```

4 . Installation of h5py
```bash
cd ../h5py
export CC=mpicc
python setup.py configure --mpi --hdf5=~/Applications/anaconda # Again, you have to adapt this
python setup.py build
python setup.py install
```


## Citation
1. Beomsoo Han, Yifeng Liu, Simon Ginzinger, and David Wishart. (2011) SHIFTX2: significantly improved protein chemical shift prediction. Journal of Biomolecular NMR, Volume 50, Number 1, 43-57. doi: 10.1007/s10858-011-9478-4.


## License
MIT