#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Simple analyze by computing RMSD and ploting predicted vs experimental chemical shift """

import os
import re
import sys
import argparse

import h5py
import numpy as np
import pandas as pd

from Bio.SeqUtils import seq3
from matplotlib import pyplot as plt

__author__ = "Jérôme Eberhardt, Roland H Stote, and Annick Dejaegere"
__copyright__ = "Copyright 2016, Jérôme Eberhardt"
__credits__ = ["Jérôme Eberhardt", "Roland H Stote", "Annick Dejaegere"]

__lience__ = "MIT"
__maintainer__ = "Jérôme Eberhardt"
__email__ = "qksoneo@gmail.com"

def get_obs(obs_file):
    return pd.read_csv(obs_file, delim_whitespace=True, usecols=[1,2,3,5], 
                       names=['resid', 'resn', 'atom_type', 'obs'])

def get_random_coil(rc_file):
    rc = pd.read_csv(rc_file, usecols=[1,2,3], names=['resn', 'atom_type', 'rc'])
    rc['resn'] = rc['resn'].apply(seq3).str.upper()
    return rc

def get_pred(hdf_file, atom_types):
    columns = ['resid', 'resn', 'atom_type', 'pred', 'std']
    df = pd.DataFrame(np.nan, index=[0], columns=columns)

    with h5py.File(hdf_file, 'r') as f:

        i = 0

        for residue in f.iteritems():
            resid = int(residue[0].split('_')[0])
            resn = str(residue[0].split('_')[1])

            for atom_type in atom_types:

                try:
                    nset = '%s_%s/%s' % (resid, resn, atom_type)
                    dset = f[nset]

                    df.loc[i] = [resid, resn, atom_type, np.nanmean(dset), np.nanstd(dset)]

                    i += 1
                except KeyError:
                    pass

    return df

def rmsd(x, y):
    diff = x - y
    diff = diff.dropna()
    return np.sqrt(np.sum(diff**2)/(float(diff.shape[0])))

def extract_dssp_info(dssp_file, sectype='H', shift = 0):
    
    start_read = False
    list_sec = []
    i = 0
    
    with open(dssp_file) as f:
        for line in f:
            if '  #  RESIDUE' in line:
                start_read = True
            
            if start_read and not "#" in line:
                resnum = int(line[6:10].strip()) + shift
                secstr = str(line[16])
                
                if sectype == secstr:
                    i += 1
                else:
                    if i: list_sec.append((resnum-i, i))
                    i = 0
    
    return list_sec

def plot_dssp_info(ax, resid, min_data, max_data, dssp_file, baseline=None):
    
    range_data = max_data - min_data
    
    if baseline != None:
        y_min = baseline - (range_data / 10.)
    else:
        y_min = min_data - (range_data / 10.)
        
    y_bar = y_min + (np.abs(min_data - y_min) / 2.)
    y_box = y_bar + (np.abs(y_bar - y_min) / 2.)
    y_len = 2. * (y_bar - y_box)
    
    ax.broken_barh(extract_dssp_info(dssp_file, 'H'), (y_box, y_len), facecolors='SteelBlue', zorder=3)
    ax.broken_barh(extract_dssp_info(dssp_file, 'E'), (y_box, y_len), facecolors='Gold', zorder=3)
    ax.plot(resid, len(resid)*[y_bar], color='black', zorder=1, lw=0.5)

def plot_shift(cs, dssp_file=None):

    color = 'chocolate'
    atom_types = cs['atom_type'].unique()

    xmin = cs['resid'].min()
    xmax = cs['resid'].max()

    label = {'CA': 'C_{\\alpha}', 'CB': 'C_{\\beta}'}

    for atom_type in atom_types:

        if label.has_key(atom_type):
            ylabel = label[atom_type]
        else:
            ylabel = atom_type

        fig, ax = plt.subplots(figsize=(12, 4))

        tmp = cs.loc[cs['atom_type'] == atom_type]

        ax.axhline(linewidth=1, color='Gainsboro')
        ax.plot(tmp['resid'], tmp['dpred'], linestyle='-', marker='.', color=color)
        ax.plot(tmp['resid'], tmp['dobs'], linestyle='-', marker='.', color='black')

        std_neg = tmp['dpred'] - tmp['std']
        std_pos = tmp['dpred'] + tmp['std']

        ax.fill_between(tmp['resid'], std_neg, std_pos, color=color, alpha=0.4)

        ymax = np.max(np.abs(tmp['dobs'])) * 1.5
        ymin = -ymax

        if dssp_file:
            plot_dssp_info(ax, tmp['resid'], ymin / 1.5, ymax, dssp_file)

        ax.set_xlabel("Residue", fontsize=20)
        ax.set_ylabel(r"$\Delta\delta %s$ (ppm)" % ylabel, fontsize=20)
        ax.set_xlim(xmin - 1, xmax + 1)
        ax.set_ylim(ymin, ymax)

        plt.savefig("shift_%s.png" % atom_type, dpi=300, format='png', bbox_inches='tight')

        plt.close(fig)

def plot_distribution(cs, hdf_file):

    color = 'chocolate'
    atom_types = cs['atom_type'].unique()

    with h5py.File(hdf_file, 'r') as f:

        for atom_type in atom_types:

            out_dir = 'shift_distribution/%s' % atom_type
            os.makedirs(out_dir)

            cols = 10
            num = cs['atom_type'].value_counts()[atom_type]
            rows = int(num / cols) + (num % cols > 0)

            figB, axB = plt.subplots(rows, cols, figsize=(35, 45), sharex=True, sharey=True)

            icol, irow = 0, 0

            xmax = 0

            # We want the global xmax !
            for residue in f.iteritems():
                resid = int(residue[0].split('_')[0])
                resn = str(residue[0].split('_')[1])

                try:
                    nset = '%s_%s/%s' % (resid, resn, atom_type)
                    dset = f[nset]

                    row = cs[(cs['resid'] == resid) & (cs['resn'] == resn) & (cs['atom_type'] == atom_type)]

                    x = dset - row['rc'].values[0]
                    xmax = np.nanmax([xmax, np.nanmax(np.abs(x)), np.abs(row['dobs']).values[0]])

                except KeyError:
                    pass

            xmax *= 1.5

            for residue in f.iteritems():
                resid = int(residue[0].split('_')[0])
                resn = str(residue[0].split('_')[1])

                try:
                    nset = '%s_%s/%s' % (resid, resn, atom_type)
                    dset = f[nset]

                    y, x = np.histogram(dset, bins=50)

                    row = cs[(cs['resid'] == resid) & (cs['resn'] == resn) & (cs['atom_type'] == atom_type)]

                    # Get Delta_delta CS
                    x = x[:-1] - row['rc'].values[0]
                    # Get density (something is wrong with numpy density ...)
                    y = y / np.float(np.sum(y))
                    
                    # For the small picture !
                    fig, ax = plt.subplots(figsize=(6, 2))
                    ax.plot(x, y, linestyle='-', color=color)
                    ax.plot(row['dobs'], 0, marker='^', markersize=30, color='black')
                    ax.plot(row['dpred'], 0, marker='^', markersize=25, color=color)

                    ax.set_ylim(0, 0.20)
                    ax.set_xlim(-xmax, xmax)
                    ax.set_ylabel('Population', fontsize=20)
                    ax.set_xlabel(r'$\Delta\delta$ (ppm)', fontsize=20)

                    fig_name = "%s/%s_%s.png" % (out_dir, resid, resn)
                    fig.savefig(fig_name, dpi=300, format='png', bbox_inches='tight')

                    plt.close(fig)

                    # And now the BIG picture !
                    axB[irow, icol].plot(x, y, linestyle='-', color=color)
                    axB[irow, icol].plot(row['dobs'], 0, marker='^', markersize=30, color='black')
                    axB[irow, icol].plot(row['dpred'], 0, marker='^', markersize=25, color=color)

                    axB[irow, icol].set_ylim(0, 0.20)
                    axB[irow, icol].set_xlim(-xmax, xmax)
                    axB[irow, icol].set_title("%s - %s" % (resid, resn))

                    if icol == 0:
                        axB[irow, icol].set_ylabel('Population', fontsize=15)

                    if irow == rows - 1:
                        axB[irow, icol].set_xlabel(r'$\Delta\delta$ (ppm)', fontsize=15)

                    icol += 1

                    if icol == cols:
                        icol = 0
                        irow += 1

                except KeyError:
                    pass

            fig_name = "%s/distribution_%s_all.png" % (out_dir, atom_type)
            figB.savefig(fig_name, dpi=300, format='png', bbox_inches='tight')

            plt.close(figB)

def plot_shift_diff(cs, dssp_file=None):

    color = 'chocolate'
    atom_types = cs['atom_type'].unique()
    label = {'CA': 'C_{\\alpha}', 'CB': 'C_{\\beta}'}

    for atom_type in atom_types:

        if label.has_key(atom_type):
            ylabel = label[atom_type]
        else:
            ylabel = atom_type

        tmp = cs.loc[cs['atom_type'] == atom_type]

        std_neg = tmp['diff'] - tmp['std']
        std_neg[std_neg < 0.] = 0.
        std_pos = tmp['diff'] + tmp['std']

        fig, ax = plt.subplots(figsize=(12,4))

        ax.plot(tmp['resid'], tmp['diff'], color=color)
        ax.fill_between(tmp['resid'], std_neg, std_pos, color=color, alpha=0.4)

        if dssp_file:
            plot_dssp_info(ax, tmp['resid'], 0, 10., dssp_file)

        ax.set_ylim(-1, 10.)
        ax.set_xlim(np.min(tmp['resid'] - 1), np.max(tmp['resid']) + 1)

        ax.set_xlabel('Residue', fontsize=20)
        ax.set_ylabel(r'$\delta %s$ |$\delta_{\mathrm{pred}}-\delta_{\mathrm{exp}}$| (ppm)' % ylabel, fontsize=20)

        fig_name = "diff_shift_%s.png" % (atom_type)
        fig.savefig(fig_name, dpi=300, format='png', bbox_inches='tight')

        plt.close(fig)

def replace_bfactors(cs, column, pdb_file, filename='new_pdb', na_value=-1):

    atom_types = cs['atom_type'].unique()

    for atom_type in atom_types:

        new_pdb_file = '%s_%s.pdb' % (filename, atom_type)
        tmp = cs.loc[cs['atom_type'] == atom_type]

        with open(pdb_file, 'r') as f, open(new_pdb_file, 'w') as w:
            for line in f:

                if re.search('^ATOM', line):
                    resid = int(line[22:26])
                    
                    row = tmp.loc[(tmp['resid'] == resid) & (tmp['atom_type'] == atom_type)]

                    try:
                        value = row[column].values[0]

                        if np.isnan(value):
                            value = na_value

                    except:
                        value = na_value

                    line = line[0:60] + '%6.2f' % value + line[66:-1]
                    w.write('%s\n' % line) 

                else:
                    w.write(line)

def main():
    parser = argparse.ArgumentParser(description='Analyze NMR chemical shift')
    parser.add_argument('-c', "--obs", dest='obsfile', required = True, \
                        action="store", type=str, \
                        help = "bmrb file with experimental chemical shift")
    parser.add_argument("-h5", dest='h5file', required = True, \
                        action = "store", type=str, \
                        help = "provide a hdf5 file with all nmr data")
    parser.add_argument("-d", "--dssp", dest='dsspfile', \
                        action = "store", type=str, default = None, \
                        help = "provide a dssp file for secondary structure")
    parser.add_argument("--distribution", dest='distribution', \
                        action = "store_true", default = False, \
                        help = "if we want to plot all the distribution")
    parser.add_argument("-p", '--pdb', dest='pdb_file', \
                        action = "store", default = None, \
                        help = "pdb file (bfactors replaced by shift diff)")

    options = parser.parse_args()
    
    obs_file = options.obsfile
    hdf_file = options.h5file
    dssp_file = options.dsspfile
    distribution = options.distribution
    pdb_file = options.pdb_file

    try:
        shiftx2_dir = os.environ['SHIFTX2_DIR']
        rc_file = shiftx2_dir + '/lib/RandomCoil.csv'
    except KeyError:
        print 'Error: The environ variable SHIFTX2_DIR is not set !'
        sys.exit(1)

    # Get obs and random coil
    cs_obs = get_obs(obs_file)
    # Get mean and std pred
    cs_pred = get_pred(hdf_file, cs_obs['atom_type'].unique())
    # Get Random coil values
    rc = get_random_coil(rc_file)

    # Merge all the data together
    cs = pd.merge(cs_obs, cs_pred, how='right', on=['resid', 'resn', 'atom_type'])
    cs = cs.sort_values(by='resid').reset_index(drop=True)
    cs = pd.merge(cs, rc, how='left', on=['resn', 'atom_type'])

    # Get Delta_delta chemical shift
    cs['dobs'] = cs['obs'] - cs['rc']
    cs['dpred'] = cs['pred'] - cs['rc']
    # Get diff between pred and obs
    cs['diff'] = np.abs(cs['pred'] - cs['obs'])

    # Compute RMSD between obs and pred for each atom type
    for atom_type in cs['atom_type'].unique():
        tmp = cs.loc[cs['atom_type'] == atom_type]
        print '%3s %5.3f' % (atom_type,  rmsd(tmp['dobs'], tmp['dpred']))

    # Save and plot data
    cs.to_csv('shift_resume.csv', index=False, na_rep='NaN')
    plot_shift(cs, dssp_file)
    plot_shift_diff(cs, dssp_file)

    if pdb_file:
        replace_bfactors(cs, 'diff', pdb_file, 'diff_shift', na_value=-1)

    if distribution:
        plot_distribution(cs, hdf_file)

if __name__=="__main__":
    main()
    sys.exit(0)