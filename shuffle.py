#!/usr/bin/env python
###############################################################################
#
# shuffle.py - shuffle unpaird fastq files
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
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################

import argparse
import sys
import os

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ContigParser:
    """Main class for reading in and parsing contigs"""
    def __init__(self): pass

    def readfq(self, fp): # this is a generator function
        """https://github.com/lh3"""
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq
                for l in fp: # search for the start of the next record
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            name, seqs, last = last[1:].split()[0], [], None
            for l in fp: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                yield name, ''.join(seqs), None # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs); # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    yield name, seq, None # yield a fasta record instead
                    break


###############################################################################
###############################################################################
###############################################################################
###############################################################################

def doWork( args ):
    import mimetypes
    try:
        GM_open = open
        try:
            # handle gzipped files
            mime = mimetypes.guess_type(args.forward)
            if mime[1] == 'gzip':
                import gzip
                GM_open = gzip.open
        except:
            print "Error when guessing contig file mimetype"
            raise
        f = GM_open(args.forward, "r")
        r = GM_open(args.reverse, "r")
        if args.shuffled[len(args.shuffled)-2:] == "gz":
            s = gzip.open(args.shuffled, "w")
        else:
            s = open(args.shuffled, "w")
        f_CP = ContigParser()
        r_CP = ContigParser()
        Rs = r_CP.readfq(r)
        for f_cid,f_seq,f_qual in f_CP.readfq(f):
            
            r_cid,r_seq,r_qual = Rs.next()
            s.write("@%s\n" % f_cid)
            s.write("%s\n" % f_seq)
            if f_qual is not None:
                s.write("+\n%s\n" % f_qual)
            s.write("@%s\n" % r_cid)
            s.write("%s\n" % r_seq)
            if r_qual is not None:
                s.write("+\n%s\n" % r_qual)
        f.close()
        r.close()
        s.close()

    except:
        print "Could not parse contig file:",args.forward,sys.exc_info()[0]
        raise

    return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('forward', help="R1")
    parser.add_argument('reverse', help="R2")
    parser.add_argument('shuffled', help="name of shuffled file")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

