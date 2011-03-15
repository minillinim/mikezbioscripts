#!/usr/bin/perl -w
###############################################################################
#
#    barcodeByMers.pl
#
#
#
#
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

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;
use IO::File;
use Data::Dumper;

#CPAN modules

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params and print copyright
my $options = checkParams();
if(!exists $options->{'silent'}) { printAtStart(); }

# global variables

# override defaults
my $kmer_length = 4;
if(exists $options->{'k'}) { $kmer_length = $options->{'k'}; }
my $window_size = 2000;
if(exists $options->{'w'}) { $window_size = $options->{'w'}; }

# abutting kmers and windows by default.
my $kmer_offset = $kmer_length;
my $window_offset = $window_size;
if(exists $options->{'ko'}) { $kmer_offset = $options->{'ko'}; }
if(exists $options->{'wo'}) { $window_offset = $options->{'wo'}; }

# by default, use the entire read
my $use_all_read = 1;
if(exists $options->{'no_lo'}) { $use_all_read = 0; }

# make the kmer array
my @mer_array = makeMers($kmer_length, 1, 1);

# do no averages by default
my $do_aves = 0;
my $aves_fh;
if(exists $options->{'ave'}) {
    $aves_fh = IO::File->new($options->{'ave'}, "w");
    $do_aves = 1;
    printOutHeaderRow($aves_fh);
}

# open the input file
my $in_file = IO::File->new($options->{'in'}, "r");
# open the output file
my $out_file = IO::File->new($options->{'out'}, "w");
printOutHeaderRow($out_file);

my $id = "";
my $sequence = "";
my $in_fasta = 0;

while(my $line = <$in_file>)
{
    chomp($line);
    if($line =~ />/gi)
    {
        if(0 == $in_fasta)
        {
            # first sequence
            $id = $line;
            $in_fasta = 1;
        }
        else
        {
            # start of a new sequence
            # deal with the old one first....
            # do GC
            $id =~ s/>//gi;
            $id =~ s/ /_/gi;
            if(0 != printSeqBarcodes($id, cutWindows($sequence), $out_file, $aves_fh)) { print "Sequence: [$id] is too short\n"; }

            # prepare for the new one
            $id = $line;
            $sequence = "";
        }
    }
    else
    {
        $sequence .= $line;
    }
}

if($in_fasta)
{
    # do the last one...
    $id =~ s/>//gi;
    $id =~ s/ /_/gi;
    if(0 != printSeqBarcodes($id, cutWindows($sequence), $out_file, $aves_fh)) { print "Sequence: [$id] is too short\n"; }
}

# close the files
$in_file->close();
$out_file->close();
if(1 == $do_aves)
{
    $aves_fh->close();
}

###############################################################################
 # SUBS
###############################################################################
sub printOutHeaderRow
{
    #
    # print the header row for the csv files
    #
    my ($fh) = @_;
    print $fh "\"SequenceID\"";
    foreach my $kmer (@mer_array)
    {
        print $fh ",\"$kmer\"";
    }
    print $fh "\n";
}

sub printSeqBarcodes {
    my ($id, $barcodes_ref, $fh, $ave_fh) = @_;
    my $window_counter = 0;
    my @aves_array = ();
    if(0 ==  scalar @{ $barcodes_ref }) { return 1; }
    
    # intialise the averages array
    if(1 == $do_aves)
    {
        foreach my $kmer (@mer_array)
        {
            push @aves_array, 0;
        }
    }
    
    foreach my $window_barcode_ref (@$barcodes_ref)
    {
        print $fh "$id"."__$window_counter";
        my $i = 0;
        foreach my $current_count (@$window_barcode_ref)
        {
            print $fh ",$current_count";
            if( 1 == $do_aves)
            {
                $aves_array[$i] += $current_count;
                $i++;
            }
        }
        print $fh "\n";
        $window_counter++;
    }
    
    if(1 == $do_aves)
    {
        print $ave_fh "$id"."__AVERAGES";
        foreach my $count_total (@aves_array)
        {
            print $ave_fh ",".(int($count_total/$window_counter));
        }
        print $ave_fh "\n";
    }
    return 0;
}

sub cutWindows {
    #
    # tap into them global vars!
    # cut the sequence into windows and pass them along
    # and then make barcodes for ALL!
    #
    my ($sequence) = @_;
    my $seq_length = length $sequence;
    # do nothing if the window if too short
    my @sequence_barcodes = ();
    if($seq_length < $window_size) { return \@sequence_barcodes; }
    my $sub_start =  0;
    while($sub_start + $window_size <= $seq_length)
    {
        push @sequence_barcodes, cutMersNBarcode(substr $sequence, $sub_start, $window_size);
        $sub_start += $window_offset;
    }
    if(1 == $use_all_read)
    {
        $sub_start -= $window_offset;
        if($sub_start + $window_size != $seq_length)
        {
            push @sequence_barcodes, cutMersNBarcode(substr $sequence, ($seq_length - $window_size), $window_size);
        }
    }
    return \@sequence_barcodes;
}

sub cutMersNBarcode {
    #
    # cut each window in kmers and do barcodes
    #
    my ($window) = @_;
    my $window_length = length $window;
    my $sub_start = 0;
    
    # get a fresh map to do counting in
    my $mer_map_ref = merArray2Map();
    
    # cut into kmers and add to the map
    while($sub_start + $kmer_length <= $window_length)
    {
        # look out for non ACGT chars!
        my $this_mer = substr $window, $sub_start, $kmer_length;
        if(exists $$mer_map_ref{$this_mer})
        {
            $$mer_map_ref{$this_mer}++;
        }
        $sub_start += $kmer_offset;
    }
    $sub_start -= $kmer_offset;
    if($sub_start + $kmer_length != $window_length)
    {
        # look out for non ACGT chars!
        my $this_mer = substr $window, $sub_start, $kmer_length;
        if(exists $$mer_map_ref{$this_mer})
        {
            $$mer_map_ref{$this_mer}++;
        }
    }
    
    # return the barcode for this window
    return merMap2Barcode($mer_map_ref);
}

sub merArray2Map {
    #
    # return an map made from the @mer_array
    #
    my %return_map = ();
    foreach my $this_mer (@mer_array)
    {
        $return_map{$this_mer} = 0;
    }
    return \%return_map;
}

sub merMap2Barcode {
    #
    # make an array barcode from a hash of kmer counts
    #
    my ($mer_map_ref) = @_;
    my @barcode = ();
    foreach my $this_mer (@mer_array)
    {
        push @barcode, ($$mer_map_ref{$this_mer});
    }
    return \@barcode;
}

sub makeMers {
    #
    # main wrapper for the working part of the kmer list making algorithm
    # $kmer_length and $sort must be either 1 or 0
    #
    my ($mer_length, $sort, $laurenize) = @_;
    my ($mer_array_ref, $gc_array_ref) = makeMersRec($mer_length - 1);
    my %output_hash = ();

    if((1 == $laurenize) and (0 == ($mer_length % 2)))
    {
        # if this option is set then the user would like to "laurenize" the data
        # remove all but the lexicographically lowest versions of each kmer.
        # already done if the k-mer length is odd
        for my $i (0 .. scalar @{ $mer_array_ref } - 1)
        {
            if(isLexiLowerThan(@$mer_array_ref[$i], revComp(@$mer_array_ref[$i])))
            {
               $output_hash{@$mer_array_ref[$i]} = @$gc_array_ref[$i];
            }
        }
    }
    else
    {
        for my $i (0 .. scalar @{ $mer_array_ref } - 1)
        {
            $output_hash{@$mer_array_ref[$i]} = @$gc_array_ref[$i];
        }        
    }


    # the output_hash now holds a (perhaps laurenized) hash of mers => gc count
    my @output_array = ();
    
    # sort by GC if need be
    if(1 == $sort)
    {
        foreach my $key (sort { $output_hash{$a} <=> $output_hash{$b} }( keys(%output_hash) ) ) {
            push @output_array, $key;
        }
    }
    else
    {
        foreach my $key ( keys(%output_hash) ) {
            push @output_array, $key;
        }
    }
    
    return @output_array;
}

sub makeMersRec {
    #
    # recursive function
    # makes a kmer list of a given length
    # 
    my ($mer_length) = @_;

    my @this_rounds_array = ();
    my @this_rounds_gc = ();
    my @alphabet = ('A', 'C', 'G', 'T');
    my @GC = (0,1,1,0);

    if($mer_length == 0) {
        return (\@alphabet, \@GC);
    }
    else
    {
        # get the lower array
        my ($mer_array_ref, $gc_array_ref) = makeMersRec($mer_length - 1);
        foreach my $alpha_index (0 .. $#alphabet) 
        {
            foreach my $mer_index ( 0 .. (scalar @{ $mer_array_ref } - 1) )
            {
                push @this_rounds_array, ($alphabet[$alpha_index].@$mer_array_ref[$mer_index]);
                push @this_rounds_gc, ($GC[$alpha_index] + @$gc_array_ref[$mer_index]);
            }     
        }
        return (\@this_rounds_array, \@this_rounds_gc);
    }
}

sub revComp {
# reverse compliment a sequence
    my ($input) = @_;
    my $output = "";
    my $current_char;
    my $j = 0; 
    for(my $i = ((length $input) - 1); $i >= 0; $i--)
    { 
        $current_char = substr ($input, $i, 1);
        if($current_char eq "a" || $current_char eq "A") { $output .=  "T"; }
        elsif($current_char eq "c" || $current_char eq "C") { $output .=  "G"; }
        elsif($current_char eq "g" || $current_char eq "G") { $output .=  "C"; }
        elsif($current_char eq "t" || $current_char eq "T") { $output .=  "A"; }
        $j++; 
    }
    return $output;
}

sub isLexiLowerThan {
# return 1 if A is strictly lexicograpically lower than B
# else 0
    my ($read_a, $read_b) = @_;
    for(my $i = 0; $i <= (length $read_a); $i++)
    {
        if((substr ($read_a, $i, 1)) lt (substr ($read_b, $i, 1)))
        {
            return 1;
        }
    }
    return 0;
}

sub checkParams 
{
    my @standard_options = ( "help+", "in:s", "out:s", "ave:s", "k:i", "w:i", "ko:i", "wo:i", "no_lo+", "silent+" );
    my %options;
    
    # Add any other command line options, and the code to handle them
    GetOptions( \%options, @standard_options );
    
    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));
    
    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};
    
    # needs an input file
    #
    exec("pod2usage $0") if ( !(exists $options{'in'}) );
    
    # needs an output file
    #
    exec("pod2usage $0") if ( !(exists $options{'out'}) );
    
    return \%options;
}

sub printAtStart {
print<<"EOF";
----------------------------------------------------------------
 $0
 Copyright (C) 2010 Michael Imelfort

 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
----------------------------------------------------------------
EOF
}

__DATA__

=head1 NAME

    barcodeByMers.pl

=head1 COPYRIGHT

   copyright (C) 2010 Michael Imelfort

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

  barcodeByMers.pl
  
  Given a multiple fasta file, calculate the length of the sequences, and
  also the gc_content. Also calculate k-mer coverages for a given set of kmers.
  Output everything to a CSV file

=head1 SYNOPSIS

    barcodeByMers.pl -in fasta_file -out output_file [-ave averages_file] [-k kmer_length] [-w window_size] [-ko kmer_offset] [-wo window_offset] [-no_lo] [-help] [-silent]

    -in fasta_file            The file to work on
    -out output_file          The file to write to
    -ave averages_file        Filename to write whole contig average barcodes to
    -k kmer_length            The kmer length to use [default: 4]
    -w window_size            Window size to bin kmers in [default: 2000]
    -ko kmer_offset           Offset to move kmers along in window [default: kmer_length]
    -wo window_offset         Offset to move the window along by
    -no_lo                    Discard the end of the contig if it is not a whole window length [default: false]
                              ie. Default -> Make the last window start from the end of the read and work backwards (use all of the read)
    [-silent]                 Output nothing extra to the screen
    [-help]                   Displays basic usage information

=cut


