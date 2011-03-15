#!/usr/bin/perl
################################################################################
#
#    File: cutDown_ILL.pl
#    Version: 0.1
#
#    Cuts down a shuffled ILL type fasta file into shorter reads.
#    The new cutdown fasta file will be called "FILESTEM_cut_xx.extension" or
#    where FILESTEM is the stem name of the input file and xx is the length
#    of the reads in the output file.
#    
#    Thus the input file being cut down is "FILESTEM.extension"
#    
#    Each read is cut into a number of smaller reads. To identify the source of
#    the smaller reads the fasta header is modified. If an input read is called:
#        
#        >bob_read
#        
#    Then the first 4 reads cut from bob will be:
#            
#    >1_bob_read
#    >2_bob_read
#    >3_bob_read
#    >4_bob_read
#    
#    Copyright (C) 2010 2011 Michael Imelfort
#
#    This file is part of the SaSSY assembler project
#
#    SaSSY is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SaSSY is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with SaSSY. If not, see <http://www.gnu.org/licenses/>.
#
#    Description:
#
#    this file reads from the DataManagerDefs.conf file
#    and produces the DataManagerDefs.cpp and DataManagerDefs.h files
#    for use with the DataManager object found in DataManager.cpp/.h
#    You will only need to call this script when you update the
#    DataManagerDefs.conf file
# 
#    Original Authors: Mike Imelfort 2010 2011
#
################################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;
use List::Util qw[min max];

#CPAN modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params and print copyright
printAtStart();
my $options = checkParams();

# set defaults and override if called for
my $offset = 1;
my $extension = "fasta";
my $read_length = $options->{'outrl'};
if( exists $options->{'offset'} )   { $offset = $options->{'offset'}; }
if( exists $options->{'extension'} ) { $extension = $options->{'extension'}; }

# set filenames
my $inputfile = $options->{'filestem'}.".".$extension;
my $outputfile = $options->{'filestem'}."_cut_$read_length.$extension";

# just a holder string
my $seq_string;

# try open the input / ouput files
open my $IN_FILE, "<", $inputfile or die $!;
open my $OUT_FILE, ">", $outputfile or die $!;

# do all the work!
cut_reads();

# clean up all the open files
close $IN_FILE;
close $OUT_FILE;

#################################################################################################################
sub cut_reads {
    my $num_processed = 0;          # total number processed
    my $total = 0;
    
    while (<$IN_FILE>)
    {
        $num_processed++;
        my $forward_head = $_;
        my $forward_seq = <$IN_FILE>;
        my $reverse_head = <$IN_FILE>;
        my $reverse_seq = <$IN_FILE>;

        chomp $forward_head;
        chomp $reverse_head;
        chomp $forward_seq;
        chomp $reverse_seq;
        
        $forward_head =~ s/>//;
        $reverse_head =~ s/>//;
        $forward_seq =~ s/ //g;
        $reverse_seq =~ s/ //g;

        # how many kmers can we cut from each guy?
        my $for_len = length $forward_seq;
        my $rev_len = length $reverse_seq;

        # we need these guys
        my $mnsl = $for_len;
        my $mxsl= $rev_len;
        my $forward_is_smallest = 1;
        if($rev_len < $for_len) { $forward_is_smallest = 0; $mnsl = $rev_len; $mxsl = $for_len;}

        # if the read length is larger than the shortest
        # there's nothing we can do
        next if($read_length > $mnsl);

        # easy to cut kmers
        my $num_kmers = int (($mnsl - $read_length)/$offset) + 1;

        my $for_remaining = $for_len - $read_length;
        my $rev_remaining = $rev_len - $read_length;

        my $off_counter = 0;
        while($off_counter < $num_kmers)
        {
            print $OUT_FILE ">".($off_counter + 1)."_$forward_head\n".(substr ($forward_seq, ($off_counter * $offset), $read_length))."\n";
            print $OUT_FILE ">".($off_counter + 1)."_$reverse_head\n".(substr ($reverse_seq, ($off_counter * $offset), $read_length))."\n";
            $total++;
            $for_remaining -= $offset;
            $rev_remaining -= $offset;
            $off_counter++;
        }

        # fix these guys
        $for_remaining += $offset;
        $rev_remaining += $offset;

        if(exists $options->{'last'})
        {
            # first make sure that at least one of them is at 0.
            if($for_remaining != 0 && $rev_remaining != 0)
            {
                if($for_remaining < $offset)
                {
                    print $OUT_FILE ">".($off_counter + 1)."_L_$forward_head\n".(substr ($forward_seq, $for_len - $read_length, $read_length))."\n";
                    $for_remaining = 0;
                }
                else
                {
                    print $OUT_FILE ">".($off_counter + 1)."_L_$forward_head\n".substr(($forward_seq, ($off_counter * $offset), $read_length))."\n";
                    $for_remaining -= $offset;
                }
                if($rev_remaining < $offset)
                {
                    print $OUT_FILE ">".($off_counter + 1)."_L_$reverse_head\n".(substr ($reverse_seq, $rev_len - $read_length, $read_length))."\n";
                    $rev_remaining = 0;
                }
                else
                {
                    print $OUT_FILE ">".($off_counter + 1)."_L_$reverse_head\n".(substr ($reverse_seq, ($off_counter * $offset), $read_length))."\n";
                    $rev_remaining -= $offset;
                }
                $total++;
                
                $off_counter++;
                
                if(exists $options->{'all'})
                {
                    # now snap up the remainder of the longest
                    # if they were the same length, these would both be at 0 now
                    if($for_remaining != 0 || $rev_remaining != 0)
                    {
                        if($forward_is_smallest)
                        {
                            while($rev_remaining >= $offset)
                            {
                                print $OUT_FILE ">".($off_counter + 1)."_A_$forward_head\n".(substr ($forward_seq, rand($for_len - $read_length), $read_length))."\n";
                                print $OUT_FILE ">".($off_counter + 1)."_A_$reverse_head\n".(substr ($reverse_seq, ($off_counter * $offset), $read_length))."\n";
                                $rev_remaining -= $offset;  
                                $off_counter++;
                                $total++;
                                
                            }
                            if($rev_remaining != 0)
                            {
                                print $OUT_FILE ">".($off_counter + 1)."_A_$forward_head\n".(substr ($forward_seq, rand($for_len - $read_length), $read_length))."\n";
                                print $OUT_FILE ">".($off_counter + 1)."_A_$reverse_head\n".(substr ($reverse_seq, $rev_len - $read_length, $read_length))."\n";
                                $total++;
                                
                            }
                        }
                        else
                        {
                            while($for_remaining >= $offset)
                            {
                                print $OUT_FILE ">".($off_counter + 1)."_A_$forward_head\n".(substr ($forward_seq, ($off_counter * $offset), $read_length))."\n";
                                print $OUT_FILE ">".($off_counter + 1)."_A_$reverse_head\n".(substr ($reverse_seq, rand($rev_len - $read_length), $read_length))."\n";
                                $for_remaining -= $offset;  
                                $off_counter++;
                                $total++;
                                
                            }
                            if($rev_remaining != 0)
                            {
                                print $OUT_FILE ">".($off_counter + 1)."_A_$forward_head\n".(substr ($forward_seq, $for_len - $read_length, $read_length))."\n";
                                print $OUT_FILE ">".($off_counter + 1)."_A_$reverse_head\n".(substr ($reverse_seq, rand($rev_len - $read_length), $read_length))."\n";
                                $total++;
                                
                            }
                        }
                    }
                }
            }
        }
    }
    
    print "----------------------------------------------------------------\n";
    print "  Report: $outputfile\n";
    print "----------------------------------------------------------------\n";
    print "  Processed:\t$num_processed\n";
    print "  Total reads made:\t$total\n";
    print "----------------------------------------------------------------\n";
    
}

#################################################################################################################
sub checkParams {
    my @standard_options = ( "help+", "man+", "outrl:i", "filestem:s", "extension:s", "offset:i", "last+", "all+");
    my %options;

    # Add any other command line options, and the code to handle them
    GetOptions( \%options, @standard_options );

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # If the -man option is set, run perldoc for this
    #
    exec("perldoc $0") if $options{'man'};
    
    # be sure that all the information is correct...
    if (!(exists $options{'outrl'}))
    { 
        print "You need to specify the readlength for the output file\n"; 
        exec("pod2usage $0"); 
    }
    if (!(exists $options{'filestem'}))
    { 
        print "You need to specify the input file\n"; 
        exec("pod2usage $0"); 
    }
    
    # all implies last
    if (exists $options{'all'}) { $options{'last'} = 1; }
    
    return \%options;
}

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2011 Michael Imelfort
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

   cutDown_ILL.pl

=head1 COPYRIGHT

   Copyright (c) 2010 Michael Imelfort
   All rights reserved.

=head1 DESCRIPTION

   Cuts down a shuffled ILL type fasta file into shorter reads.
   The new cutdown fasta file will be called "FILESTEM_cut_xx.extension" or
   where FILESTEM is the stem name of the input file and xx is the length
   of the reads in the output file.

   Thus the input file being cut down is "FILESTEM.extension"

   Each read is cut into a number of smaller reads. To identify the source of the
   smaller reads the fasta header is modified. If an input read is called:

   >bob_read

   Then the first 4 reads cut from bob will be:

   >1_bob_read
   >2_bob_read
   >3_bob_read
   >4_bob_read

=head1 SYNOPSIS

cutDown_ILL.pl -outrl INT -filestem STRING [-offset INT] [-extension STRING] [-last] [-all] [-help] [-man]

    BASIC

      -filestem STRING     The name of the input file (minus extension)
      -outrl INT           The readlength for the output file
      -offset INT          The offset between cut reads [default: 1]
      -extension STRING    The extension of the file - fasta, fna, etc... [default: fasta]
      
    ADVANCED
    
      -last                Cut the "last" kmer from a read regardless of if it matches the offset. 
                           Ie. Given 65bp reads with -outrl = 63 and -offset = 4. Using "-last" would
                           give 2 kmers; starting at the first and third bases. Without "-last" the program
                           would only give the first and the remaining two bases would be lost.
      -all                 Use when the input reads are different lengths. Selects from the longest according
                           to the offset given, but selects random kmers from the shorter to make up the pair.
                           -all implies -last
    MISC
    
      -help              Displays basic usage information
      -man               Displays more detailed information

         
=cut

