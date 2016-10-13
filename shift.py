#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Extract all the chemical shift from MD trajectory using SHIFTX+ from SHIFTX2 """

import os
import re
import sys
import glob
import shlex
import shutil
import fnmatch
import argparse
import subprocess

import h5py
import numpy as np
import pandas as pd

from mpi4py import MPI
from MDAnalysis import Universe

__author__ = "Jérôme Eberhardt, Roland H Stote, and Annick Dejaegere"
__copyright__ = "Copyright 2016, Jérôme Eberhardt"
__credits__ = ["Jérôme Eberhardt", "Roland H Stote", "Annick Dejaegere"]

__lience__ = "MIT"
__maintainer__ = "Jérôme Eberhardt"
__email__ = "qksoneo@gmail.com"

#pd.set_option('display.max_rows', None)

class Shiftx:
    def __init__(self, pdb_dir=None, top_file=None, dcd_file=None):

        # Define MPI environnment
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

        if self.rank == 0:
            if pdb_dir and (top_file or dcd_file):
                print 'Warning: cs will be extracted from pdb and not from dcd'
            elif not pdb_dir and not (top_file and dcd_file):
                print 'Error: pdb file or (top and dcd files) are required'
                self.comm.Abort()

        self.pdb_dir = pdb_dir
        self.top_file = top_file
        self.dcd_file = dcd_file

        # And check if SHIFTX2_DIR environnment variable is set
        if self.rank == 0:
            try:
                self.shiftx2_dir = os.environ['SHIFTX2_DIR']
            except KeyError:
                print 'Error: The environ variable SHIFTX2_DIR is not set!'
                self.comm.Abort()
        else:
            self.shiftx2_dir = None

        self.shiftx2_dir = self.comm.bcast(self.shiftx2_dir, root=0)

    def get_fname(self, dir, s):
        return os.path.join(dir, s)

    def get_files(self, dir, pattern):
        
        files = fnmatch.filter(os.listdir(dir), pattern)
        #files = [os.path.abspath(self.getFname(dir, sn)) for sn in files]
        files = [self.get_fname(dir, sn) for sn in files]
        RE_DIGIT = re.compile(r'(\d+)')
        ALPHANUM_KEY = lambda s: [int(g) if g.isdigit() else g for g in RE_DIGIT.split(s)]
        files.sort(key = ALPHANUM_KEY)
        
        return files

    def get_range(self, length, interval):

        # Si process 0
        if self.rank == 0:
            # On genere une liste de la taille de la simulation
            N = np.arange(start=0, stop=length, step=interval, dtype=np.int32)
            # On la coupe en part egale
            N_part = np.array_split(N, self.size)

            # On supprime la liste N, car elle ne sert plus a rien
            del(N)

            # Si on a plus de 1 processus
            if self.size > 1:
                # On envoie un message a chaque process
                for part, dest in zip(N_part[1:], np.arange(1, self.size)):
                    self.comm.send([part[0], part[-1]], dest=dest)

            # Puis on retourne le range du process 0
            tmp = [N_part[0][0], N_part[0][-1]]
            return tmp[0], tmp[1]

        # Si process non 0
        elif self.rank > 0:
            # On retourne le message recu du process 0
            tmp = self.comm.recv(source=0)
            return tmp[0], tmp[1]

    def run_shiftx(self, directory, pH=7.4, temperature=298.15):
        cmd = 'java -cp %s/bin:%s/lib/weka.jar ShiftXp -b \'%s/*.pdb\' -atoms ALL -ph %s -temp %s -dir %s' 
        cmd = cmd % (self.shiftx2_dir, self.shiftx2_dir, directory, pH, temperature, self.shiftx2_dir)

        args = shlex.split(cmd)

        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, errors = p.communicate()

        return output, errors

    def read_sxp_files(self, sxp_files):

        df = None

        for sxp_file in sxp_files:
            if df is None:
                df = pd.read_csv(sxp_file)
            else:
                df = pd.merge(df, pd.read_csv(sxp_file), on=['NUM', 'RES', 'ATOM', 'IUPAC'], 
                              how='outer', suffixes=['', 'x'])

        df = df.sort_values(by='NUM').reset_index(drop=True)
        df.rename(columns=lambda x: x.replace('x', ''), inplace=True)

        return df

    def create_dataset(self, sxp, size, shift=0, filename='shiftx.hdf5'):
        with h5py.File(filename, 'w', driver='mpio', comm=self.comm) as f:
            for index, row in sxp.iterrows():
                nset = '%d_%s/%s' % (row['NUM']+shift, row['RES'], row['ATOM'])
                dset = f.create_dataset(nset, (size,), dtype=np.float32)


    def fill_dataset(self, sxp, start=0, shift=0, filename='shiftx.hdf5'):
        with h5py.File(filename, 'r+', driver='mpio', comm=self.comm) as f:
            for index, row in sxp.iterrows():
                nset = '%d_%s/%s' % (row['NUM']+shift, row['RES'], row['ATOM'])

                try:
                    value = row['PRED_CS'].values.astype(np.float32)
                    stop = start + value.shape[0]

                except AttributeError:
                    value = np.float32(row['PRED_CS'])
                    stop = start + 1

                try:
                    dset = f[nset]
                    dset[start:stop,] = value
                except KeyError:
                    print "Warning: unable to access to %s" % nset
                    pass

    def run(self, pH=7.4, temperature=298.15, interval=1, shift=0):

        if self.pdb_dir:

            # Monoprocess for pdb extraction
            if self.rank == 0:

                # Fire off shiftx2!!
                self.run_shiftx(self.pdb_dir, pH=pH, temperature=temperature)
                # Get list of all the cs file
                sxp_files = self.get_files(self.pdb_dir, '*.sxp')
                # Read all the sxp files
                sxp = self.read_sxp_files(sxp_files)

                # Create dataset and fill it directly
                self.create_dataset(sxp, len(sxp_files), filename='shiftx.hdf5')
                self.fill_dataset(sxp, filename='shiftx.hdf5')

                # And we remove every thing!
                [os.remove(f) for f in glob.glob('%s/*.sxp' % self.pdb_dir)]

        else:

            try:
                # Open trajectories!
                u = Universe(self.top_file, self.dcd_file)
            except:
                print 'Error: failed to load trajectory!'
                self.comm.Abort()
                
            total_frames = len(u.trajectory)
            used_frames = (total_frames / interval)

            if self.rank == 0:
                if not total_frames % interval == 0:
                    print 'Error: The interval is not a multiple of the total number of frame !'
                    self.comm.Abort()

                print "Total number of frames: %s (%s frames used)" % (total_frames, used_frames)

            start, stop = self.get_range(total_frames, interval)

            j = 1
            created = False
            tmp_dir = 'tmp_shiftx'
            tmp_rank = '%s/tmp_%s' % (tmp_dir, self.rank)

            # Save the orginal path
            ori_path = os.getcwd()

            # Create directories
            os.makedirs(tmp_rank)
            os.chdir(tmp_rank)
            
            for i in xrange(start, stop+interval, interval):

                # Go to the frame i
                u.trajectory[i]
                # And extract the PDB
                # The method zfill() pads string on the left with zeros to fill width.
                u.atoms.write("snapshot_%s.pdb" % str(i).zfill(len(str(total_frames))))

                # We extract 5000 structures at a time
                if (j % 5000) == 0 or stop == i:

                    # Fire off shiftx2 !!
                    self.run_shiftx('.', pH=pH, temperature=temperature)

                    # Get list of all the cs file
                    sxp_files = self.get_files('.', '*.sxp')
                    # Read all the sxp files
                    sxp = self.read_sxp_files(sxp_files)

                    # Create dataset
                    if not created:
                        self.create_dataset(sxp, used_frames, shift, '%s/shiftx.hdf5' % ori_path)
                        created = True

                    # And fill it with data !
                    numstruct = np.int(sxp_files[0].split('_')[-1].split('.')[0])
                    start_at = np.int(numstruct / interval)
                    self.fill_dataset(sxp, start_at, shift, '%s/shiftx.hdf5' % ori_path)

                    # And we remove every thing !
                    [os.remove(pdb_file) for pdb_file in os.listdir('.')]

                    j = 1

                j += 1

            # And we go back in dir_output and clean
            os.chdir(ori_path)

            if self.rank == 0:
                shutil.rmtree(tmp_dir)

def parse_options():

    parser = argparse.ArgumentParser(description = 'Run SHIFTX+ with MPI')
    parser.add_argument('-p', "--pdb", dest='pdb_dir',
                        action="store", type=str,
                        help="directory with pdb files")
    parser.add_argument('-t', "--top", dest='top_file',
                        action="store", type=str,
                        help="top file used for simulation (HS[EDP] --> HIS)")
    parser.add_argument('-d', "--dcd", dest='dcd_file',
                        action="store", type=str, nargs='+',
                        help="one dcd file or a list of dcd files")
    parser.add_argument('-i', "--interval", dest='interval',
                        action="store", type=int, default=1,
                        help="used frames at this interval")
    parser.add_argument('-s', "--shift", dest='shift',
                        action="store", type=int, default=0,
                        help="number to shift residu number")
    parser.add_argument('--pH', dest='pH',
                        action="store", type=float, default=7.4,
                        help="pH")
    parser.add_argument('--temp', dest='temperature',
                        action="store", type=float, default=298.15,
                        help="temperature in Kelvin")

    return parser.parse_args()

def main():
    # Parse the passing arguments
    options = parse_options()

    pdb_dir = options.pdb_dir
    top_file = options.top_file
    dcd_file = options.dcd_file
    pH = options.pH
    temperature = options.temperature
    shift = options.shift
    interval = options.interval

    S = Shiftx(pdb_dir, top_file, dcd_file)
    S.run(pH, temperature, interval, shift)

if __name__ == '__main__':
    main()
    sys.exit(0)