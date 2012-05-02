#!/usr/bin/env perl
###############################################################################
#
#    barcodeByMers.pl
#    
#    Create kmer barcodes for a set of sequences
#
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
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;
use Carp;
use threads;
use threads::shared;
 
#CPAN modules
use Bio::SeqIO;
use Data::Dumper;

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# edit here to log all external commands 
my $global_log_commands = 0;

# ext command failure levels
use constant {
    IGNORE_FAILURE => 0,
    WARN_ON_FAILURE => 1,
    DIE_ON_FAILURE => 2
};

# get input params and print copyright
my $global_options = checkParams();
if(!exists $global_options->{'silent'}) { printAtStart(); }

######################################################################
# CODE HERE
######################################################################
# override defaults
my $global_kmer_length = overrideDefault(4, 'kmer_length');      # default kmer length of 4
my $global_window_size = overrideDefault(2000, 'window_size');
my $global_cut_off_len = overrideDefault($global_window_size , 'cutoff');

# abutting kmers and windows by default.
my $global_kmer_offset = overrideDefault(1,'kmer_offset');
my $global_window_offset = overrideDefault($global_window_size, 'window_offset');

# by default, use the entire read
my $global_use_all_read = 1;
if(exists $global_options->{'no_lo'}) { $global_use_all_read = 0; }

# threading! -> default of two threads, one for file parsing and one for 
# munging. Any more than 2 and the exess all get used for munging
my $global_num_threads = overrideDefault(2, 'threads');
my $global_working_threads = $global_num_threads-1;

if(!exists $global_options->{'silent'}) {
print<<EOF
    Making barcodes:
----------------------------------------------------------------
    Kmer length:                       $global_kmer_length mers
    Kmer offset:                       $global_kmer_offset bp 
    Window size:                       $global_window_size bp
    Window offset:                     $global_window_offset bp
    Rejecting all seqs shorter than:   $global_cut_off_len bp
    Using:                             $global_num_threads threads
----------------------------------------------------------------
EOF
}

# make the kmer array
my @global_mer_array = makeMers($global_kmer_length, 1);

# do no averages by default
my $global_do_aves = 0;
my $global_aves_fh;
if(exists $global_options->{'ave'}) {
    $global_aves_fh = openWrite($global_options->{'ave'});
    $global_do_aves = 1;
    printOutHeaderRow($global_aves_fh);
}

# a place to store good sequences -> shared across threads!
my @global_good_seqs = ();
my $global_good_seqs_counter = 0;
my $global_still_parsing = 1;

# keep count of where we're at
my $global_update_amount = 1000;

my $global_valid_seqs = 0;
my $global_parsed_seqs = 0;
my $global_complete_seqs = 0;
my $global_rejected_seqs = 0;

share (@global_good_seqs);
share ($global_good_seqs_counter);
share ($global_still_parsing);
share ($global_valid_seqs);

# open the output file
my $global_out_fh = openWrite($global_options->{'out'});
printOutHeaderRow($global_out_fh);

# one thread for parsing through the seqio object
my $seqio = Bio::SeqIO->new(-file => $global_options->{'in'}, '-format' => 'Fasta');
my $parse_thread = threads->new(\&packSeqNIDs, $seqio);
sleep(2);

# IO management
my $update_count = $global_update_amount;
my @super_full = ();
my @super_aves = ();

# do seqs in batches
my $batch_size = 5;
$batch_size *= 2;

# keep parsing as long as we have to
my $keep_parsing = 0;
{
    lock($global_valid_seqs);
    if($global_valid_seqs > $global_parsed_seqs) { $keep_parsing = 1; }
}
while(1 == $keep_parsing)
{
    # start N-1 threads at a time
    my @all_threads = ();
    my $num_started_this_loop = 0;
    foreach my $thread_counter (1..$global_working_threads) {
        {
            # get the next good sequence
            lock(@global_good_seqs);
            lock($global_good_seqs_counter);
            if($global_good_seqs_counter > $#global_good_seqs) {
                my $keep_waiting = 0;
                {
                    lock($global_still_parsing);
                    $keep_waiting = $global_still_parsing;
                }
                
                while(1 == $global_still_parsing) {
                    # give the parser a chance to catch up
                    sleep(1);
                    # check if we have something to do yet
                    last if($global_good_seqs_counter <= $#global_good_seqs);
                    # perhaps we can go through this list again...
                    {
                        lock($global_still_parsing);
                        $keep_waiting = $global_still_parsing;
                    }
                }
                {
                    lock($global_still_parsing);
                    if((0 == $global_still_parsing) and ($global_good_seqs_counter > $#global_good_seqs))
                    {
                        $thread_counter = $global_working_threads;
                        $keep_parsing = 0;
                        last;
                    }
                }
            }
            if(1 == $keep_parsing)
            {
                # got a sequence to parse!
                my $seq_id = $global_good_seqs[$global_good_seqs_counter];
                $global_good_seqs_counter++;
                my $seq_string = $global_good_seqs[$global_good_seqs_counter];
                $global_good_seqs_counter++;

                # cut barcodes on a new thread
                $num_started_this_loop++;
                $global_parsed_seqs++;
                #print "starting: $seq_id $global_parsed_seqs\n";
                
                # DO NOT REMOVE THIS TWO STEP THREAD START!
                # SEE: http://www.perlmonks.org/?node_id=265269
                my ($thr) = threads->new(\&makeSeqBarcodes, ($seq_id, $seq_string));
                $all_threads[$thread_counter] = $thr;
            }
        }
    }
    
    foreach my $thread_counter (1..$num_started_this_loop) {
        # get the next thread and wait for him to complete
        my $this_thread = $all_threads[$thread_counter];
        my ($full_bcode, $aves_bcode) = $this_thread->join;
    
        # Save on IO ops
        push @super_full, $full_bcode;
        push @super_aves, $aves_bcode;
           
        # print print and update the user
        $global_complete_seqs++;
        $update_count--;
        if(0 == $update_count) {
            if(!exists $global_options->{'silent'}) { print "Processed: $global_complete_seqs\n"; }
    
            # print now!
            foreach my $line (@super_full) { print $global_out_fh $line; }
            foreach my $line (@super_aves) { print $global_aves_fh $line; }
                
            $#super_full = -1;
            $#super_aves = -1;
                
            $update_count = $global_update_amount;
        }
    }
    # another round?
    {
        lock($global_valid_seqs);
        if($global_valid_seqs > $global_parsed_seqs) { $keep_parsing = 1; }
    }  
}

# print remaining!
foreach my $line (@super_full) { print $global_out_fh $line; }
foreach my $line (@super_aves) { print $global_aves_fh $line; }

# close this thread off
$parse_thread->join;

# close the files
close($global_out_fh);
if(1 == $global_do_aves)
{
    close($global_aves_fh);
}

if(!exists $global_options->{'silent'})
{
print<<EOF
    Processed: $global_complete_seqs sequences
    Rejected: $global_rejected_seqs sequences
----------------------------------------------------------------
EOF
}

######################################################################
# CUSTOM SUBS
######################################################################
sub packSeqNIDs {
    #-----
    # Go through the fasta file and pack the global
    # arrays with seqs and IDs
    #
    my ($seqio) = @_;
    while(my $seq = $seqio->next_seq)
    {
        my $seq_string = $seq->seq;
        next if ((length $seq_string < $global_window_size) or (length $seq_string < $global_cut_off_len));
        my $seq_id = $seq->id;
        $seq_id =~ s/ /_/gi;
        # push this onto the global list of good sequences
        {
            lock(@global_good_seqs);
            push @global_good_seqs, $seq_id;
            push @global_good_seqs, $seq_string;
            {
                lock($global_valid_seqs);
                $global_valid_seqs++;
                #print "--> $global_valid_seqs\n";
            }
            
        }
    }
    
    # let the other threads know that we're done!
    {
        lock($global_still_parsing);
        $global_still_parsing = 0;
    }
}

sub printOutHeaderRow
{
    #-----
    # print the header row for the csv files
    #
    my ($fh) = @_;
    print $fh "\"SequenceID\"";
    foreach my $kmer (@global_mer_array)
    {
        print $fh ",\"$kmer\"";
    }
    print $fh "\n";
}

sub makeSeqBarcodes {
    my ($id, $seq_string) = @_;
    my $barcodes_ref = cutWindows($seq_string);
    my $window_counter = 0;
    my @aves_array = ();
    
    # return the barcode via these two strings
    my $aves_bcode = "";
    my $full_bcode = "";
    
    if(0 ==  scalar @{ $barcodes_ref }) { return 1; }
    
    # intialise the averages array
    if(1 == $global_do_aves)
    {
        foreach my $kmer (@global_mer_array)
        {
            push @aves_array, 0;
        }
    }
    
    foreach my $window_barcode_ref (@$barcodes_ref)
    {
        $full_bcode .= "$id"."__$window_counter";
        my $i = 0;
        foreach my $current_count (@$window_barcode_ref)
        {
            $full_bcode .= ",$current_count";
            if( 1 == $global_do_aves)
            {
                $aves_array[$i] += $current_count;
                $i++;
            }
        }
        $full_bcode .= "\n";
        $window_counter++;
    }
    
    if(1 == $global_do_aves)
    {
        $aves_bcode .= "$id";
        foreach my $count_total (@aves_array)
        {
            $aves_bcode .= sprintf(",%.6f",($count_total/$window_counter));
        }
        $aves_bcode .= "\n";
    }
    return ($full_bcode, $aves_bcode);
}

sub cutWindows {
    #-----
    # tap into them global vars!
    # cut the sequence into windows and pass them along
    # and then make barcodes for ALL!
    #
    my ($sequence) = @_;
    my $seq_length = length $sequence;

    # do nothing if the window if too short
    my @sequence_barcodes = ();
    if($seq_length < $global_window_size) { return \@sequence_barcodes; }
    my $sub_start =  0;
    while($sub_start + $global_window_size <= $seq_length)
    {
        push @sequence_barcodes, cutMersNBarcode(substr $sequence, $sub_start, $global_window_size);
        $sub_start += $global_window_offset;
    }
    if(1 == $global_use_all_read)
    {
        $sub_start -= $global_window_offset;
        if($sub_start + $global_window_size != $seq_length)
        {
            push @sequence_barcodes, cutMersNBarcode(substr $sequence, ($seq_length - $global_window_size), $global_window_size);
        }
    }
    return \@sequence_barcodes;
}

sub cutMersNBarcode {
    #-----
    # cut each window in kmers and do barcodes
    #
    my ($window) = @_;
    my $window_length = length $window;
    my $sub_start = 0;
    
    # get a fresh map to do counting in
    my $mer_map_ref = merArray2Map();
    
    # cut into kmers and add to the map
    while($sub_start + $global_kmer_length <= $window_length)
    {
        # look out for non ACGT chars!
        my $this_mer = lowlexi(substr $window, $sub_start, $global_kmer_length);
        if(exists $$mer_map_ref{$this_mer})
        {
            $$mer_map_ref{$this_mer}++;
        }
        $sub_start += $global_kmer_offset;
    }
    $sub_start -= $global_kmer_offset;
    if($sub_start + $global_kmer_length != $window_length)
    {
        # look out for non ACGT chars!
        my $this_mer = lowlexi(substr $window, $sub_start, $global_kmer_length);
        if(exists $$mer_map_ref{$this_mer})
        {
            $$mer_map_ref{$this_mer}++;
        }
    }
    
    # return the barcode for this window
    return merMap2Barcode($mer_map_ref);
}

sub merArray2Map {
    #-----
    # return an map made from the @global_mer_array
    #
    my %return_map = ();
    foreach my $this_mer (@global_mer_array)
    {
        $return_map{$this_mer} = 0;
    }
    return \%return_map;
}

sub merMap2Barcode {
    #-----
    # make an array barcode from a hash of kmer counts
    #
    my ($mer_map_ref) = @_;
    my @barcode = ();
    foreach my $this_mer (@global_mer_array)
    {
        push @barcode, ($$mer_map_ref{$this_mer});
    }
    return \@barcode;
}

sub makeMers {
    #
    # main wrapper for the working part of the kmer list making algorithm
    # $global_kmer_length and $sort must be either 1 or 0
    #
    my ($mer_length, $sort) = @_;
    my ($mer_array_ref, $gc_array_ref) = makeMersRec($mer_length - 1);

    # remove all but the lexicographically lowest versions of each kmer.
    my %output_hash = ();
    for my $i (0 .. scalar @{ $mer_array_ref } - 1)
    {
        $output_hash{lowlexi(@$mer_array_ref[$i])} = @$gc_array_ref[$i];
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
    #-----
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

sub revcompl {
    #-----
    # Reverse complement a sequence
    #
    my ($seq) = @_;
    $seq =~ tr/ACGTNacgtn/TGCANtgcan/;
    return scalar reverse $seq;
}

sub lowlexi {
    #-----
    # Return the lowest lexicographical form of a sequence
    #
    my ($seq) = @_;
    my $rev_seq = revcompl($seq);
    if($seq lt $rev_seq) { return $seq; }
    else { return $rev_seq; }
}

######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS
  
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help+", "threads|t:i", "in|i:s", "out|o:s", "ave|a:s", "kmer_length|k:i", "window_size|w:i", "kmer_offset|M:i", "window_offset|D:i", "cutoff|c:i", "no_lo+", "silent+" );
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsory items
    if(!exists $options{'in'} ) { printParamError ("We need an input file to process!"); }
    if(!exists $options{'out'} ) { printParamError ("We need an output file to process!"); }
    #if(!exists $options{''} ) { printParamError (""); }

    return \%options;
}

sub printParamError
{
    #-----
    # What to do if there's something wrong with a parameter
    #  
    my ($error) = @_;  
    print "**ERROR: $0 : $error\n"; exec("pod2usage $0");
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}

######################################################################
# FILE IO

sub openWrite
{
    #-----
    # Open a file for writing
    #
    my ($fn) = @_;
    open my $fh, ">", $fn or croak "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{   
    #-----
    # Open a file for reading
    #
    my ($fn) = @_;
    open my $fh, "<", $fn or croak "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}

######################################################################
# EXTERNAL COMMANDS
#
# checkAndRunCommand("ls", {
#                          -a => ""
#                          }, 
#                          WARN_ON_FAILURE);

sub checkFileExists {
    #-----
    # Does a file exists?
    #
    my ($file) = @_;
    unless(-e $file) {
        croak "**ERROR: $0 : Cannot find:\n$file\n";
    }
}

sub logExternalCommand
{
    #-----
    # Log a command line command to the command line!
    #
    if(1 == $global_log_commands) {
        print $_[0], "\n";
    }
}

sub isCommandInPath
{
    #-----
    # Is this command in the path?
    #
    my ($cmd, $failure_type) = @_;
    if (system("which $cmd |> /dev/null")) {
        handleCommandFailure($cmd, $failure_type);
    }
}

sub runExternalCommand
{
    #-----
    # Run a command line command on the command line!
    #
    my ($cmd) = @_;
    logExternalCommand($cmd);
    system($cmd);
}

sub checkAndRunCommand
{
    #-----
    # Run external commands more sanelier
    # 
    my ($cmd, $params, $failure_type) = @_;
    
    isCommandInPath($cmd, $failure_type);
    
    # join the parameters to the command
    my $param_str = join " ", map { $_ . " " . $params->{$_}} keys %{$params};
    my $cmd_str = $cmd . " " . $param_str;
    
    logExternalCommand($cmd_str);

    # make sure that all went well    
    if (system($cmd_str)) {
         handleCommandFailure($cmd_str, $failure_type)
    }
}

sub handleCommandFailure {
    #-----
    # What to do when all goes bad!
    #
    my ($cmd, $failure_type) = @_;
    if (defined($failure_type)) {
        if ($failure_type == DIE_ON_FAILURE) {
            croak "**ERROR: $0 : " . $! . "\n";
        } elsif ($failure_type == WARN_ON_FAILURE) {
            carp "**WARNING: $0 : "  . $! . "\n";
        }
    }
}

######################################################################
# MISC

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) Michael Imelfort
    
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

   copyright (C) Michael Imelfort

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

  Given a multiple fasta file, calculate the length of the sequences, and
  also the gc_content. Also calculate k-mer coverages for a given set of kmers.
  Output everything to a CSV file
  
=head1 SYNOPSIS

    barcodeByMers.pl -in|i FILE -out|o FILE

    -i -in FILE               The file to work on
    -o -out FILE              The file to write to
    [-a -ave FILE]            Filename to write whole contig average barcodes to [default: NO averages]
    [-t -threads INT]         The number of threads to use [default: 2] 
                              NOTE: For reasons not worth going into here, it's generally
                              worth your while to run this guy with many threads
                              Try using 2 cores less than the capacity of your machine.
                              It won;t hurt nothing!
    [-c -cutoff INT]          Reject all sequences shorter than this [default: WINDOW_SIZE]
    [-k -kmer_length INT]     The kmer length to use [default: 4]
    [-M -kmer_offset INT]     Offset to move kmers along in window [default: 1]
    [-w -window_size INT]     Window size to bin kmers in [default: 2000]
    [-D -window_offset INT]   Offset to move the window along by [default: WINDOW_SIZE]
    [-no_lo]                  Discard the end of the contig if it is not a whole window length [default: false]
                              ie. Default -> Make the last window start from the end of the read and work backwards (use all of the read)
    [-silent]                 Output nothing extra to the screen
    [-help]                   Displays basic usage information
         
=cut
