#!/usr/bin/env python
###############################################################################
#
# gm2esom.py - make groopm bins into an esom class file!
#
#
# You need 2 files: 
#
# 'gmbin' file is the output of groopm print -f minimal XXX.gm
# looks like:
#
# 1    hydrocarbon_scaffold_100496    7188
# 1    hydrocarbon_scaffold_101405    12398
# 1    hydrocarbon_scaffold_101829    8589
# 1    hydrocarbon_scaffold_102672    32830
# 1    hydrocarbon_scaffold_102777    22555
# ...
#
# *** make sure you've trimmed all the header stuff out of this file
#
# *PLUS* yr esom names file
#
#
# Produces an ESOM class file which looks like:
#
# 26860%
# 0%      NOCLASS 255     255     255
# 1%      1       255     255     0
# 2%      2       0       0       255
# 3%      3       255     0       255
# 4%      4       0       255     0
# 5%      5       0       255     255
# 6%      6       255     128     0
# 7%      7       128     0       255
# 8%      8       153     114     63
# ...
#
#
# 1       15
# 2       15
# 3       15
# 4       0
# 5       15
# 6       5
# ...
#
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2012"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################

import argparse
import sys
import os
from colorsys import hsv_to_rgb as htr
import numpy as np

###############################################################################
###############################################################################
###############################################################################
###############################################################################

  # classes here

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def doWork( args ):
    """ Main wrapper"""
    
    esom_name_to_ids = []
    num_esom_names = 0
    all_named_contigs = {}
    # parse names file and make the hash we'll be needing
    try:
        with open(args.names, "r") as fh:
            for line in fh:
                num_esom_names += 1
                parts = line.rstrip().split("\t")
                cid_parts = parts[1].split('_')
                if 'leftover' in cid_parts:
                    base_name = "_".join(cid_parts[:len(cid_parts)-2])
                else:
                    base_name = "_".join(cid_parts[:len(cid_parts)-1])
                esom_name_to_ids.append((base_name, parts)) # ('hydrocarbon_scaffold_82282', ['1', 'hydrocarbon_scaffold_82282_0']) 
                all_named_contigs[base_name] = True
    except: 
        print "Error opening file:", args.names, sys.exc_info()[0]
        raise    

    # parse the GM file
    gm_bin_ids = {}
    cid_2_gmbin = {}
    try:
        with open(args.gmbin, "r") as fh:
            for line in fh:
                parts = line.rstrip().split("\t")
                bid = int(parts[0])
                cid_2_gmbin[parts[1]] = bid 
                gm_bin_ids[bid] = True                 
    except: 
        print "Error opening file:", args.names, sys.exc_info()[0]
        raise    

    num_gm_bins = len(gm_bin_ids.keys())

    # produce a whole heap of colors
    color_steps = [float(i)/num_gm_bins for i in range(num_gm_bins)]
    S = 1       # SAT and VAL remain fixed at 1. Reduce to make
    V = 1       # Pastels if that's your preference...
    raw_cols = np.array([np.array(htr(val, S, V))*255 for val in color_steps])
    bin_cols = [[int(i) for i in j] for j in raw_cols]

    # work out all the bin 0 contigs
    for cid in all_named_contigs:
        if cid not in cid_2_gmbin:
            cid_2_gmbin[cid] = 0
    
    # build the class file
    print "%d" % num_esom_names + "%"
    print "0%\tNOCLASS\t255\t255\t255"
    for i in range(len(bin_cols)):
        print "%d" % (i+1) + "%" + "\t%d\t%d\t%d\t%d" % (i+1, bin_cols[i][0],bin_cols[i][1],bin_cols[i][2])

    for i in range(len(esom_name_to_ids)):
        # ('hydrocarbon_scaffold_82282', ['1', 'hydrocarbon_scaffold_82282_0']) 
        gm_bid = cid_2_gmbin[esom_name_to_ids[i][0]]
        print "%d\t%d" % (i+1,gm_bid)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('gmbin', help="output of groopm print")
    parser.add_argument('names', help="esom names file")
    
    # parse the arguments
    args = parser.parse_args()        

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
