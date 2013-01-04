#!/usr/bin/python
import subprocess
import os
from optparse import OptionParser
import sys
import tempfile
###############################################################################
#
#    makesam.py
#    Wrapper to produce a .sam file for use in pair plotting and general chicanery
#    Copyright (C) 2010-2012 Michael Imelfort
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
    subprocess.call('bwa index -a '+ algorithm+' '+ database, shell=True)

def aln(database, readfile, outfile, threads):
    subprocess.call('bwa aln -t '+ threads+' '+ database+' '+ readfile+' >'+outfile, shell=True)

def sampe(database, sai_1, sai_2, readfile_1, readfile_2, outfile):
    if outfile is None:
        subprocess.call('bwa sampe ' + database+ ' ' + sai_1+ ' ' + sai_2+ ' ' +  readfile_1+ ' ' +  readfile_2, shell=True)
    else:
        subprocess.call('bwa sampe '+database+' '+sai_1+' '+sai_2+' '+readfile_1+' '+readfile_2+' >'+outfile, shell=True)

def samse(database, sai_1, readfile_1, outfile):
    if outfile is None:
        subprocess.call('bwa samse '+database+' '+sai_1+' '+readfile_1, shell=True)
    else:
        subprocess.call('bwa samse '+database+' '+sai_1+' '+readfile_1+' >'+outfile, shell=True)

def bwasw(database, readfile_1, readfile_2, outfile, threads):
    if outfile is None:
        if readfile_2 is None:
            subprocess.call('bwa bwasw -t '+threads+' '+database+' '+readfile_1, shell=True)
        else:
            subprocess.call('bwa bwasw -t '+threads+' '+database+' '+readfile_1+' '+readfile_2, shell=True)
    else:
        if readfile_2 is None:
            subprocess.call('bwa bwasw -t '+threads+' '+database+' '+readfile_1+' >'+outfile , shell=True)
        else:
            subprocess.call('bwa bwasw -t '+threads+' '+database+' '+readfile_1+' '+readfile_2+' >'+outfile , shell=True)

def bwasw_to_sorted_indexed_bam(database, readfile_1, readfile_2, outfile):
    if readfile_2 is None:
        subprocess.call('bwa bwasw -t '+threads+' '+database+' '+readfile_1+' | samtools view -SubhF 4 - |samtools sort - '+outfile , shell=True)
    else:
        subprocess.call('bwa bwasw -t '+threads+' '+database+' '+readfile_1+' '+readfile_2+' | samtools view -SubhF 4 - |samtools sort - '+outfile , shell=True)
    # samtools index cannot be piped, so a tmpfile is required
    cmd2 = 'samtools index '+outfile+'.bam'
    subprocess.call(cmd2, shell=True)

def sampe_to_sorted_indexed_bam(database, sai_1, sai_2, readfile_1, readfile_2, outfile):
    cmd1 = 'bwa sampe '+database+' '+sai_1+' '+sai_2+' '+readfile_1+' '+readfile_2+' | samtools view -SubhF 4 - |samtools sort - '+outfile
    subprocess.call(cmd1, shell=True)
    # samtools index cannot be piped, so a tmpfile is required
    cmd2 = 'samtools index '+outfile+'.bam'
    subprocess.call(cmd2, shell=True)

def samse_to_sorted_indexed_bam(database, sai_1, readfile_1, outfile):
    cmd1 = 'bwa samse '+database+' '+sai_1+' '+readfile_1+' | samtools view -SubhF 4 - |samtools sort - '+outfile
    subprocess.call(cmd1, shell=True)
    # samtools index cannot be piped, so a tmpfile is required
    cmd2 = 'samtools index '+outfile+'.bam'
    subprocess.call(cmd2, shell=True)

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
    parser.add_option("-b", "--bam_filename", type="string",
            dest="bamfilename", help="Output a sorted indexed bam file, of this name")
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
            sys.stderr.write('You need to specify a multiple fasta file and ONE read file (single ended)'+"\n")
            parser.print_help()
            sys.exit(1)
    else:
        doSings = False
        if (opts.database is None or opts.readfile_1 is None ):
            sys.stderr.write('You must specify both -1 and -2 and -d for a paired alignment.  For single ended just use -1 and -d'+"\n")
            parser.print_help()
            sys.exit(1)

    # override defaults
    if(opts.algorithm is None):
        algorithm = "is"
    else:
        algorithm = opts.algorithm

    # create indexes if required
    if(opts.keptfiles is None):
        sys.stderr.write('making indices'+"\n")
        sys.stderr.flush
        mkindex(opts.database, algorithm)

    # run the actual alignment
    output_file = None
    if opts.samfilename is not None:
        if opts.samfilename.endswith('.sam'):
            output_file = opts.samfilename
        else:
            output_file = opts.samfilename + '.sam'

    bam_output_file = None
    if opts.bamfilename is not None:
        if opts.bamfilename.endswith('.bam'):
            bam_output_file = opts.bamfilename[:-4] #samtools renames the output file with .bam already
        else:
            bam_output_file = opts.bamfilename

    #create tmp files
    sai1 = tempfile.mkstemp(suffix='.sai')
    sai2 = tempfile.mkstemp(suffix='.sai')
    
    if(opts.longReads):
        if bam_output_file is None:
            bwasw(opts.database, opts.readfile_1,opts.readfile_2,
                    output_file, opts.threads)
        else:
            bwasw_to_sorted_indexed_bam(opts.database,
                    opts.readfile_1,opts.readfile_2, bam_output_file,
                    opts.threads)
    else:
        aln(opts.database, opts.readfile_1, sai1[1], opts.threads)
        if(doSings is False):
            aln(opts.database, opts.readfile_2, sai2[1], opts.threads)
            if bam_output_file is None:
                sampe(opts.database, sai1[1], sai2[1], opts.readfile_1, opts.readfile_2,
                      output_file)
            else:
                sampe_to_sorted_indexed_bam(opts.database, sai1[1], sai2[1], opts.readfile_1, opts.readfile_2,
                      bam_output_file)
        else:
            if bam_output_file is None:
                samse(opts.database, sai1[1], opts.readfile_1, output_file)
            else:
                samse_to_sorted_indexed_bam(opts.database, sai1[1], opts.readfile_1, bam_output_file)
                

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
    safeRemove(sai1[1])
    safeRemove(sai2[1])

