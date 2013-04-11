#!/usr/bin/env python
import pysam
import subprocess
import os
from optparse import OptionParser
import random
import sys
###############################################################################
#
#    sam2PairPlotCSV.py
#    Main script for producing csv files for sam pair and coverage plots.
#    Copyright (C) 2010 Adam Skarshewski, Michael Imelfort
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
# Directionality
###############################################################################
#  0  DD  <---  --->   if the first read disagrees with ref then the second is ahead and in the reverse orientation (wrt the first)
#  3  DA  <---  <---   if the first read disagrees with ref then the second is ahead and in the same orientation
#  2  AD  --->  <---   if the first read agrees with ref then the second is ahead and in the reverse orientation
#  3  AA  --->  --->   if the first read agrees with ref then the second is ahead and in the same orientation
#  4  Single hit: this read agrees with it's ref
#  5  Single hit: this read disagrees with it's ref
#  6  Different references: this read agrees with it's ref
#  7  Different references: this read disagrees with it's ref
###############################################################################
#
# Add a new entry to the Dict or append it onto the end of an existing entry
#
def addPosRevToContigDictOrParse(dictionary, key, value, parsed_dictionary, rl):
    if key in dictionary:
        # check to make sure that the position for these two guys are the same
        # don't add it twice if this is trhe case
        if(dictionary[key][0][0] != value[0]):
            (start, ins, code) = getMappingCode(dictionary[key], value, rl)
            # kill the old entry
            del dictionary[key]
            # add the parsed one...
            if start in parsed_dictionary:
                parsed_dictionary[start].append((ins, code))
            else:
                # { start_pos : (insert, code) }
                parsed_dictionary[start] = [(ins, code)]
    else:
        dictionary[key] = [value]

#
# place an int in the dictionary or add to the one there...
#
def addOrAdd(dictionary, key, value):
    if key in dictionary:
        dictionary[key] = dictionary[key] + value
    else:
        dictionary[key] = value

#
# Strip the _ or - from the end of a query name
#
def sanitiseQName(qname):
    q_len = len(qname)
    if(qname[q_len-2] == '_' or qname[q_len-2] == '-'):
        return qname[ 0:(q_len-2)]
    else:
        return qname

#
# get a mapping code for the mapped pair
#
def getMappingCode(read_1, read_2, readLen):
    # these guys are not necessarilly sorted
    # do this first
    if(read_1[0][0] < read_2[0]):
        read1 = read_1[0]
        read2 = read_2
    else:
        read2 = read_1[0]
        read1 = read_2

    ins_size = read2[0] + readLen - read1[0]
    if read1[1] and read2[1]:
        return (read1[0], ins_size, 3)
    elif read1[1] == read2[1]:
        return (read1[0], ins_size, 3)
    elif read1[1] == 0:
        return (read1[0], ins_size, 0)
    else:
        return (read1[0], ins_size, 2)

#
# Parse through the sam file and work out orientations etc for each entry
# Fill all the data structs needed for .csv file creation
#
def parseSamBam(samPath, useBinary, makeCoverage, stopAt, CSVFileName):
    # open the SAM/BAM
    if (useBinary):
        try:
            samFile = pysam.Samfile(samPath, 'rb')
        except:
            print "Unable to open BAM file -- did you supply a SAM file instead?"
            sys.exit(1)
    else:
        try:
            samFile = pysam.Samfile(samPath, 'r')
        except:
            print "Unable to open SAM file -- did you supply a BAM file instead?"
            sys.exit(1)

    num_contigs =  len (samFile.header['SQ'])
    total_size = reduce(lambda x,y:x+y, samFile.lengths)

    # take not of if the user set a stop point
    if( 0 == stopAt):
        stopAt = -1

    #
    # samFile.header['SQ'][i]['SN'] IS the contig fasta ID
    #

    # store all the mappings for each contig in it's own list
    contig_mappings = {}
    for i in range(0, num_contigs):
        contig_mappings[samFile.header['SQ'][i]['SN']] = {}

    # store all the parsed mappings in this dictionary
    parsed_mappings = {}
    for i in range(0, num_contigs):
        parsed_mappings[samFile.header['SQ'][i]['SN']] = {}

    # open the .csv files
    csv_files = {}
    for i in range(0, num_contigs):
        csv_file = open(samFile.header['SQ'][i]['SN'].replace(' ','_') + '_' + CSVFileName + '.csv', 'wb')
        csv_file.write('"position","insertsize","direction"\n')
        # Do not close the filehandle now because they are still needed. This might cause a "Too many files open" error though...
        csv_files[samFile.header['SQ'][i]['SN']] = [csv_file]

    # start parsing
    samIter = samFile.fetch()
    try:
        samRecord = samIter.next()
    except:
        print "No reads found"
        return None

    # now get the read length for all the reads
    rl = len(samRecord.seq)
    parsed_lines = 0
    total_parsed = 0

    # first get all the "Good" mappings and add them to a dictionary
    # this is a crude type of sort...
    while (0 != stopAt):
        stopAt = stopAt - 1
        parsed_lines = parsed_lines + 1
        name = sanitiseQName(samRecord.qname)
        pair = (samRecord.pos, samRecord.is_reverse)

        # if there was a mapping. Update the parsed mappings dictionary
        # { start_pos : (insert, code), [(insert, code), ...] }
        if(0 <= samRecord.rname):
            addPosRevToContigDictOrParse(contig_mappings[samFile.getrname(samRecord.rname)], name, pair, parsed_mappings[samFile.getrname(samRecord.rname)], rl)

        # keep the user in the loop
        if(2000000 <= parsed_lines):
            total_parsed = total_parsed + parsed_lines
            print str(total_parsed) + " done..."
            parsed_lines = 0

        # next record!
        try:
            samRecord = samIter.next()
        except:
            break

    print "Parsed: " + str(parsed_lines)

    # if the user has specifed that we'd like to coverage then it's time to take a quick detour
    # and do this now
    if(makeCoverage):
        # for each contig in the file
        for i in range(0, num_contigs):
            # reset these guys for the new contig
            v_height = 0            # The height at any given position
            incre_dict = {}         # Per position height increment
            decre_dict = {}         # per position height decrement
            max_occurance = 0;      # The number of times the most commonly occuring v_height occurs
            mode_coverage = 0;      # The most commonly occuring v_height
            mode_dict = {}          # for keeping track of mode

            # make the .csv file and header
            con_id = samFile.header['SQ'][i]['SN']
            con_length = samFile.header['SQ'][i]['LN']
            cov_file = open(con_id.replace(' ','_') + '_' + CSVFileName + '_coverage.csv', 'wb')
            cov_file.write('"position","coverage"\n')

            # for each position in the contig
            for position in range(1, int(con_length)+1):
                if position in parsed_mappings[con_id]:
                    # there was at least one hit here! Go through all of them
                    for hit_number in range(0, len(parsed_mappings[con_id][position])):
                        insert = parsed_mappings[con_id][position][hit_number][0]
                        # add in for the read
                        addOrAdd(incre_dict, position, 1)
                        addOrAdd(decre_dict, position + rl, 1)
                        # add in for it's pair
                        addOrAdd(incre_dict, position + insert, 1)
                        addOrAdd(decre_dict, position + rl + insert, 1)
                if position in incre_dict:
                    # add to the coverage at this position
                    v_height = v_height + incre_dict[position]
                    del incre_dict[position]
                if position in decre_dict:
                    # add to the coverage at this position
                    v_height = v_height - decre_dict[position]
                    del decre_dict[position]

                cov_file.write('{0},{1}\n'.format(position, v_height))

                # update the "mode coverage"
                if(v_height in mode_dict):
                    mode_dict[v_height] = mode_dict[v_height] + 1
                    if(mode_dict[v_height] > max_occurance):
                        max_occurance = mode_dict[v_height]
                        mode_coverage = v_height
                else:
                    mode_dict[v_height] = 1

            cov_file.close()

            # write the mode
            cov_file = open(con_id.replace(' ','_') + '_' + CSVFileName + '_mode_coverage.csv', 'wb')
            cov_file.write('"mode_coverage"\n')
            cov_file.write('{0}\n'.format(mode_coverage))
            cov_file.close()

    # now we save all the lost souls
    # All the entries in the contig_mapping dictionary will be singly mapped
    # or mapped to two separate references
    holding_dictionary = {}
    for i in range(0, num_contigs):
        con_id = samFile.header['SQ'][i]['SN']
        for key in contig_mappings[con_id].iterkeys():
            if key in holding_dictionary:
                # this guy MUST be in a different scaffold
                # and this is a mess...
                if (holding_dictionary[key][1][0][1]):
                    code = 6
                else:
                    code = 7
                parsed_mappings[holding_dictionary[key][0]][holding_dictionary[key][1][0][0]] = [(1, code)]

                if (contig_mappings[con_id][key][0][1]):
                    code = 6
                else:
                    code = 7
                parsed_mappings[con_id][contig_mappings[con_id][key][0][0]] = [(1, code)]
                del holding_dictionary[key]
            else:
                # { READ_ID : ( contigID , ( pos , reversed ) ) }
                holding_dictionary[key] = (con_id, contig_mappings[con_id][key])

        del contig_mappings[con_id]

    # collect the last of the single mapped
    for key in holding_dictionary.iterkeys():
        if (holding_dictionary[key][1][0][1]):
            code = 4
        else:
            code = 5
        parsed_mappings[holding_dictionary[key][0]][holding_dictionary[key][1][0][0]] = [(1, code)]

    # now we print .csv files and close them
    for i in range(0, num_contigs):
        con_id = samFile.header['SQ'][i]['SN']
        csv_file = csv_files[con_id][0]
        for key in sorted(parsed_mappings[con_id].iterkeys()):
            csv_file.write('{0},{1},{2}\n'.format(key, parsed_mappings[con_id][key][0][0], parsed_mappings[con_id][key][0][1]))
        csv_file.close()

    del holding_dictionary
    del parsed_mappings
    del contig_mappings

    # close sam
    samFile.close()

#
# Entry sub. Parse vars and call parseSamBam
#
if __name__ == '__main__':

    # intialise the options parser
    parser = OptionParser("%prog -s samFileName [-b] [-c] [-o CSV fileName] [-v]\n\t-- parse a sam/bam file and produce csv files")
    parser.add_option("-s", "--sam", type="string", dest="samFileName", help="Give a SAM/BAM file name")
    parser.add_option("-b", "--binary", action="store_true", dest="useBinary", help="Set this if you use a BAM file [default: false]")
    parser.add_option("-c", "--coverage", action="store_true", dest="makeCoverage", help="Set this to output coverage information too [default: false]")
    parser.add_option("-o", "--csv_fileName", type="string", dest="CSVFileName", help="Specify a name for the CSV file [default: map_out.csv]")
    parser.add_option("-N", "--number_SAM", type="int", dest="samFileStop", help="Specify how many SAM / BAM records to parse")

    # get and check options
    (opts, args) = parser.parse_args()
    if (opts.samFileName is None):
        # could be that we're using this as the cleaver
        print("Specify a SAM/BAM filename")
        parser.print_help()
        sys.exit(1)

    # override defaults
    if(opts.CSVFileName is None):
        CSV_file_name = 'map_out'
    else:
        CSV_file_name =  opts.CSVFileName

    # override defaults
    if(opts.makeCoverage is None):
        makeCoverage = False
    else:
        makeCoverage = True

    if(opts.samFileStop is None):
        stopPoint = 0
    else:
        stopPoint = opts.samFileStop

    # do stuff
    parseSamBam(opts.samFileName, opts.useBinary, makeCoverage, stopPoint,  CSV_file_name)
