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

## Installation on UNIX

I highly recommand you to install the Anaconda distribution (https://www.continuum.io/downloads) if you want a clean python environnment with nearly all the prerequisites already installed (NumPy, Pandas, Matplotlib, MPI).

In order to use hdf5 format with MPI, you need a special version compiled in parallel mode (with mpi4py).

1 . First, you need to install mpi4py
```bash
conda install mpi4py
```

2 . Now, remove hdf5 and h5py from Anaconda (trust me)
```bash
conda remove hdf5 h5py
```

3 . After, download the latest version of hdf5 and h5py
```bash
wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.17.tar # check if it is the latest version
git clone https://github.com/h5py/h5py.git
```

4 . Compilation of HDF5
```bash
tar -xvf hdf5-1.8.17.tar
cd hdf5-1.8.17
PREFIX=~/Applications/anaconda    # You need to adapt this to your system
./configure --prefix=$PREFIX --enable-linus-lfs --with-zlib=$PREFIX --enable-parallel --enable-shared
make
make install
```

5 . Installation of h5py
```bash
cd ../h5py
export CC=mpicc
python setup.py configure --mpi --hdf5=~/Applications/anaconda # Again, you have to adapt this
python setup.py build
python setup.py install
```

6 . For the rest, you just have to do this
```bash
pip install MDAnalysis
```

## How-To

1 . First you need to extract all the chemical shift from MD trajectory (or PDB files). You need at least to specify a directory with some PDB files or a topology and a dcd file (**don't forget to rename all the special residues like HS[EDP] to HIS in the topology file, because SHIFTX+ won't recognize them**). Prediction from PDB files will use only one processor, so you don't have to use MPI. On the contrary, to predict all the chemical shift from MD trajectories I recommand you to use MPI (knowing that SHIFTX+ is not the fastest software alive), especially if they are very long (> 100000 frames). **And last advice, you will need a lot of RAM because SHIFTX+ takes at minimum 700 Mb (for what, I don't know. In comparison, SHIFTX takes nothing) multiply by the number of core you use !**
```bash
python shift.py -p pdb_directory # With a single core
```
```bash
mpiexec -np 4 python shift.py -t topology.psf -d traj.dcd # With MPI using 4 cores
```
**Command line options**
* -p/--pdb: directory with pdb files
* -t/--top: topology file (pdb, psf)
* -d/--dcd: single trajectory of list of trajectories (dcd, xtc)
* --pH: pH (default: 7.4)
* --temperature: temperature (default: 7.4)
* -i/--interval: used frames at this interval (default: 1)
* -s/--shift: shift the resid number to match with experimental data (current resid + shift) (default: 0)

2 . And finally, compare them to experimental data. It will compute the RMSD between the prediction and the experimental data and plot the secondary chemical shift along the sequence for each element (Ca, Cb, etc ...). Finally, if you want, you can plot the chemical shift distribution for each residue (long operation)(see citation). But now, that you have the HDF5 file with the data, you can  play with it and run whatever analysis you want ...
```bash
python analyze.py -c obs.bmrb -h5 shiftx.hdf5 -d dssp.file
```
**Command line options**
* -c/--obs: BMRB file with experimental chemical shift
* -h5: HDF5 file with the chemical shift
* -d/--dssp: DSSP file to plot secondary structure on the plot (Default: None)
* -p/--pdb: pdb file used to put data in bfactors column (Default: None)
* --distribution: if you want to plot the chemical shift distribution for each residue (default: False)

## BMRB file format

The BMRB file that contains all the experimental chemical shift values must be like that.

```
        1 225 SER CA C  58.40 0.0  1
        2 225 SER C  C 174.90 0.0  1
        3 226 ALA H  H   8.31 0.0  1
        4 226 ALA N  N 126.10 0.0  1
        5 226 ALA CA C  53.50 0.0  1
        6 226 ALA CB C  18.20 0.0  1
        7 226 ALA C  C 178.40 0.0  1
        8 227 ASN H  H   8.26 0.0  1
        9 227 ASN N  N 117.00 0.0  1
       10 227 ASN CA C  54.00 0.0  1
       11 227 ASN CB C  38.30 0.0  1
       12 227 ASN C  C 175.90 0.0  1
       13 228 GLU H  H   7.95 0.0  1
       14 228 GLU N  N 119.20 0.0  1
       15 228 GLU CA C  57.10 0.0  1
       16 228 GLU C  C 177.00 0.0  1
```

## Citation
1. Beomsoo Han, Yifeng Liu, Simon Ginzinger, and David Wishart. (2011) SHIFTX2: significantly improved protein chemical shift prediction. Journal of Biomolecular NMR, Volume 50, Number 1, 43-57. doi: 10.1007/s10858-011-9478-4.
2. Robustelli, P., Stafford, K. A., & Palmer III, A. G. (2012). Interpreting protein structural dynamics from NMR chemical shifts. Journal of the American Chemical Society, 134(14), 6365-6374.


## License
MIT
