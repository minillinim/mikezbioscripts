#!/usr/bin/env python
###############################################################################
#
# cLens.py - Extract contig lengths and put into a JSON-encoded dict
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
__copyright__ = "Copyright 2014"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "1.0.0"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Done"

###############################################################################

import argparse
import sys
import json

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ContigParser:
    """Main class for reading in and parsing contigs"""
    def __init__(self): pass

    def readFasta(self, fp): # this is a generator function
        header = None
        seq = None
        while True:
            for l in fp:
                if l[0] == '>': # fasta header line
                    if header is not None:
                        # we have reached a new sequence
                        yield header, "".join(seq)
                    # save the header we just saw
                    header = l.rstrip()[1:].partition(" ")[0]
                    seq = []
                else:
                    seq.append(l.rstrip())
            # anything left in the barrel?
            if header is not None:
                yield header, "".join(seq)
            break

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def doWork( args ):
    """ Main wrapper"""

    # parse conting file
    contig_lens = {}
    try:
        with open(args.contigFile, "r") as fh:
            CP = ContigParser()
            for header, seq in CP.readFasta(fh):
               contig_lens[header] = len(seq)
    except:
        print "Error opening file:", args.contigFile, sys.exc_info()[0]
        raise

    try:
        with open(args.outFile, "w") as fh:
            fh.write(json.dumps(contig_lens))
    except:
        print "Error opening file:", args.outFile, sys.exc_info()[0]
        raise

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('contigFile', help="File containing contgs")
    parser.add_argument('outFile', help="File to write to")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
