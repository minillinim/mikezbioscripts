#!/usr/bin/python
import subprocess
import os
from optparse import OptionParser
import sys
import pysam
###############################################################################
#
#    plotPipe.py
#    Main pipeline for contig stats plotting
#    Copyright (C) Michael Imelfort
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

#
# Entry sub. Parse vars and call parseSamBam
#
if __name__ == '__main__':

    # intialise the options parser
    parser = OptionParser("\n\n %prog -s samfile -c contigfile [options]")
    parser.add_option("-s", "--sam_fileName", type="string", dest="samFileName", help="Specify a sam file to parse")
    parser.add_option("-N", "--number_SAM", type="int", dest="samFileStop", help="Specify how many SAM / BAM records to parse")
    parser.add_option("-b", "--binary", action="store_true", dest="useBinary", help="Set this if you use a BAM file [default: false]")
    parser.add_option("-c", "--contigs", type="string", dest="contigFileName", help="Specify a contig file to parse")
    parser.add_option("-k", "--keep_all", action="store_true", dest="keepAll", help="Set this if you want to keep the output csv files [default: false]")

    print "Plotting stage"

    # get and check options
    (opts, args) = parser.parse_args()
    if (opts.samFileName is None):
        print ('You need to specify a sam file to parse')
        parser.print_help()
        sys.exit(1)
    if (opts.contigFileName is None):
        print ('You need to specify a contigs file to parse')
        parser.print_help()
        sys.exit(1)

    print "Command line parameters look good. Parsing SAM"

    # work out how many to parse
    if (opts.samFileStop is None):
        stops = '';
    else:
        stops = '-N ' + str(opts.samFileStop)

    # parse SAM
    if (opts.useBinary is None):
        parse_cmd = "sam2PairPlotCSV.py -c -s " + opts.samFileName + ' ' + stops
    else:
        parse_cmd = "sam2PairPlotCSV.py -c -b -s " + opts.samFileName + ' ' + stops

    os.system(parse_cmd)

    print "Parsing contig file"

    # parse the contig file
    os.system("makeKmerGCCSV.pl -in "+opts.contigFileName+" -kmers 3,4,8,12")# -silent")

    print "Sedifying the output"

    # merge some of the results:
    os.system("sed -e 's/\"mode_coverage\"//' -e '/^$/d' *mode_coverage.csv > mods.pp.csv")
    os.system('echo \'"contigId","seq_length","GC","seq_length_withN","mode_cov"\' > header.pp')
    os.system("paste -d ',' "+opts.contigFileName+"_GC.csv mods.pp.csv > pasted.pp")
    os.system("cat header.pp pasted.pp > "+opts.contigFileName+"_GC.csv")
    os.system("rm header.pp")
    os.system("rm pasted.pp")
    os.system("rm mods.pp.csv")

    print "Making purdy pictures"

    # make purdy pictures
    if (opts.useBinary):
        try:
            samFile = pysam.Samfile(opts.samFileName, 'rb')
        except:
            print "Unable to open BAM file -- did you supply a SAM file?"
            sys.exit(1)
    else:
        try:
            samFile = pysam.Samfile(opts.samFileName, 'r')
        except:
            print "Unable to open SAM file -- did you supply a BAM file?"
            sys.exit(1)

    num_contigs =  len (samFile.header['SQ'])

    # store all the mappings for each contig in it's own list
    CSVFileName = "map_out"
    for i in range(0, num_contigs):
        con_id = samFile.header['SQ'][i]['SN']
        con_id_pair_csv = con_id.replace(' ','_') + '_' + CSVFileName + '.csv'
        con_id_cov_csv = con_id.replace(' ','_') + '_' + CSVFileName + '_coverage.csv'
        con_id_mode_cov_csv = con_id.replace(' ','_') + '_' + CSVFileName + '_mode_coverage.csv'
        con_id_heat_csv =  con_id.replace(' ','_') + "_KC.csv";
        cmd = "pipelinePlotter.py -p "+con_id_pair_csv+" -c "+con_id_cov_csv+" -m "+con_id_mode_cov_csv+" -H "+con_id_heat_csv+" -l " + str(con_id)
        os.system(cmd)
        if(opts.keepAll is None):
            os.system("rm "+ con_id_pair_csv)
            os.system("rm "+ con_id_cov_csv)
            os.system("rm "+ con_id_mode_cov_csv)
            os.system("rm "+ con_id_heat_csv)
