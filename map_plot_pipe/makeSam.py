#!/usr/bin/python
import subprocess
import os
from optparse import OptionParser
import sys
###############################################################################
#
#    makesam.py
#    Wrapper to produce a .sam file for use in pair plotting and general chicanery
#    Copyright (C) 2010 Michael Imelfort
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################


def mkindex(database, algorithm):
    os.system('bwa index -a '+algorithm+' '+database)

def aln(database, readfile, outfile, threads):
    os.system('bwa aln -t '+threads+' '+database+' '+readfile+' > '+outfile)

def sampe(database, sai_1, sai_2, readfile_1, readfile_2, outfile):
    if outfile is None:
        os.system('bwa sampe '+database+' '+sai_1+' '+sai_2+' '+readfile_1+' '+readfile_2)
    else:
        os.system('bwa sampe '+database+' '+sai_1+' '+sai_2+' '+readfile_1+' '+readfile_2+' > '+outfile)

def samse(database, sai_1, readfile_1, outfile):
    if outfile is None:
        os.system('bwa samse '+database+' '+sai_1+' '+readfile_1)
    else:
        os.system('bwa samse '+database+' '+sai_1+' '+readfile_1+' > '+outfile)

def bwasw(database, readfile_1, outfile):
    if outfile is None:
        os.system('bwa bwasw '+database+' '+readfile_1 )
    else:
        os.system('bwa bwasw '+database+' '+readfile_1+' >'+outfile )

def safeRemove(fileName):
    if os.path.isfile(fileName):
        os.system('rm ' + fileName)

# Entry sub. Parse vars and call parseSamBam
#
if __name__ == '__main__':

    # intialise the options parser
    parser = OptionParser("\n\n %prog [options]")
    parser.add_option("-1", "--reads_1", type="string", dest="readfile_1", help="The first data of a paired read file")
    parser.add_option("-2", "--reads_2", type="string", dest="readfile_2", help="The second data of a paired read file")
    parser.add_option("-d", "--database", type="string", dest="database", help="The scaffold, query, database...")
    parser.add_option("-a", "--bwa_algorithm", type="string", dest="algorithm", help="The algorithm bwa uses for indexing 'bwtsw' or 'is' [default: is]")
    parser.add_option("-k", "--keep", action="store_true", dest="keepfiles", help="Keep all the database index files etc after (see also --kept) [default: false]")
    parser.add_option("-K", "--kept", action="store_true", dest="keptfiles", help="Assume the indices already exist, don't re-make them (and don't delete them) (e.g. previously this script was run with -k/--keep [default: false]")
    parser.add_option("-s", "--sam_filename", type="string",
            dest="samfilename", help="The name for the final sam file name [default: STDOUT]")
    parser.add_option("-S", "--single", action="store_true", dest="singleEnd", help="Use this for non-paired reads [default: false]")
    parser.add_option("-L", "--long_reads", action="store_true",dest="longReads", help="The input is long reads (eg. 454), sets the search algorithm to BWA-SW")
    parser.add_option("-t", "--threads", type="string", dest="threads",
            default="1", help="The number of threads to use when aligning")

    # get and check options
    (opts, args) = parser.parse_args()
    if(opts.readfile_2 is None):
        # single ended!
        doSings = True
        if (opts.database is None or opts.readfile_1 is None ):
            print ('You need to specify a multiple fasta file and ONE read file (single ended)')
            parser.print_help()
            sys.exit(1)
    else:
        doSings = False
        if (opts.database is None or opts.readfile_1 is None ):
            print ('You must specify both -1 and -2 and -d for a paired alignment.  For single ended just use -1 and -d')
            parser.print_help()
            sys.exit(1)

    # override defaults
    if(opts.algorithm is None):
        algorithm = "is"
    else:
        algorithm = opts.algorithm

    # create indexes if required
    if(opts.keptfiles is None):
        print('making indices')
        mkindex(opts.database, algorithm)

    # run the actual alignment
    if(opts.longReads):
        bwasw(opts.database, opts.readfile_1, opts.samfilename)
    else:
        aln(opts.database, opts.readfile_1, "out_bwa_sa1.sai", opts.threads)
        if(doSings is False):
            aln(opts.database, opts.readfile_2, "out_bwa_sa2.sai", opts.threads)
            sampe(opts.database, "out_bwa_sa1.sai", "out_bwa_sa2.sai", opts.readfile_1, opts.readfile_2, opts.samfilename)
        else:
            samse(opts.database, "out_bwa_sa1.sai", opts.readfile_1, opts.samfilename)

    # clean up
    if(opts.keepfiles is None and opts.keptfiles is None):
        safeRemove(opts.database+'.amb')
        safeRemove(opts.database+'.ann')
        safeRemove(opts.database+'.bwt')
        safeRemove(opts.database+'.pac')
        safeRemove(opts.database+'.rbwt')
        safeRemove(opts.database+'.rpac')
        safeRemove(opts.database+'.rsa')
        safeRemove(opts.database+'.sa')
        safeRemove('out_bwa_sa1.sai')
        if(doSings is False):
            safeRemove('out_bwa_sa2.sai')

