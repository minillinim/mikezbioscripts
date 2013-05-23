#!/usr/bin/perl

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

my $options = check_params();

# set defaults and override if called for
my $offset = 1;
my $multi_rls = 0;
my $overhang_insert = 0;
my $directory = "";
my $extension = "fasta";
my $read_length = $options->{'outrl'};
if( exists $options->{'offset'} )   { $offset = $options->{'offset'}; }
if( exists $options->{'multirl'} )   { $multi_rls = $options->{'multirl'}; }
if( exists $options->{'ohangins'} )  { $overhang_insert = $options->{'ohangins'}; }
if( exists $options->{'extension'} ) { $extension = $options->{'extension'}; }
if( exists $options->{'directory'})
{
    $directory = $options->{'directory'}."/";
    `mkdir -p $directory`;
}

# set filenames
my $inputfile = $options->{'filestem'}.".".$extension;
my $outputfile = $options->{'filestem'}."_cut_".$options->{'outrl'}.".".$extension;
my $overhangfile = "";
if(0 != $overhang_insert) { $overhangfile = $options->{'filestem'}."_overhang_".$options->{'outrl'}."_".$overhang_insert.".".$extension; }

my $seq_string;

# try open the input / ouput files
open my $IN_FILE, "<", $inputfile or die $!;
open my $OUT_FILE, ">", $outputfile or die $!;
my $OHANG_FILE;

# now we can cut down
if(0 != $overhang_insert) { 
   open $OHANG_FILE, ">", $overhangfile or die $!; 
   cut_reads();
   close $OHANG_FILE;
}
else
{
    cut_reads();
}

# clean up all the open files
close $IN_FILE;
close $OUT_FILE;

#################################################################################################################
sub cut_reads {
    my $num_skipped = 0;
    while (<$IN_FILE>)
    {
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

        my $for_seq_length = length $forward_seq;
        my $rev_seq_length = length $reverse_seq;
        my $min_length = 0;
        my $num_kmers = 0;
        my $seq_used = 0;
        my $overhang_string = "";

        if ($for_seq_length < $read_length) 
        {
            if($rev_seq_length < $read_length) { $num_skipped++; next; }
            else
            {
                # reverse is the longer one, we can make no good reads
                # but if we have an overhang then we can cut it entirely from
                # the reverse strand.
                if($overhang_insert != 0)
                {
                    cut_ohang($reverse_head, $reverse_seq);
                }
            }
        }
        else
        {
            if($rev_seq_length < $read_length)
            {
                # forward is the longer one, we can make no good reads
                # but if we have an overhang then we can cut it entirely from
                # the forward strand.            
                if($overhang_insert != 0)
                {
                    cut_ohang($forward_head, $forward_seq);
                }
            }
            else
            {
                # we can cut at least 1 read from both but we'll have to work out which 
                # is the longer one.
                $min_length = min($for_seq_length, $rev_seq_length);
                $num_kmers = int(($min_length - $read_length)/$offset) + 1;
                $seq_used = ($num_kmers  - 1) * $offset + $read_length;
                if($overhang_insert != 0)
                {
                    if($min_length == $rev_seq_length) 
                    { 
                        if( ( $for_seq_length - $seq_used ) >= ($read_length + $overhang_insert) )
                        {
                            cut_ohang($forward_head, (substr $forward_seq, $seq_used));
                        }
                    }
                    else 
                    {    
                        if( ( $rev_seq_length - $seq_used ) >= ($read_length + $overhang_insert) )
                        {
                            cut_ohang($reverse_head, (substr $reverse_seq, $seq_used));
                        }
                    }
                }
            }
        }
        my $off_counter = 0;
        while($off_counter < $num_kmers)
        {
            print $OUT_FILE ">".($off_counter + 1)."_$forward_head\n".(substr ($forward_seq, ($off_counter * $offset), $read_length))."\n";
            print $OUT_FILE ">".($off_counter + 1)."_$reverse_head\n".(substr ($reverse_seq, ($off_counter * $offset), $read_length))."\n";
            $off_counter++;
        }
    }
    print "Skipped: $num_skipped\n";
}

#################################################################################################################
sub cut_ohang {
    my ($seq_head, $seq_str) = @_;
    my $num_cuts = int( ( (length $seq_str) - (2*$read_length + $overhang_insert) ) / $offset) + 1;
    my $for_off = 0;
    my $rev_off = $read_length + $overhang_insert;
    my $cut_count = 1;
    while($cut_count <= $num_cuts)
    {
        print $OHANG_FILE ">"."OH$cut_count"."_$seq_head"."_1\n".(substr $seq_str, $for_off, $read_length)."\n>"."OH$cut_count"."_$seq_head"."_2\n".(substr $seq_str, $rev_off, $read_length)."\n";
        $for_off += $offset;
        $rev_off += $offset;
        $cut_count++;
    }
}

#################################################################################################################
sub check_params {
    my @standard_options = ( "help+", "man+", "outrl:i", "filestem:s", "extension:s", "offset:i", "multirl+", "ohangins:i", "directory:s");
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
    if((exists $options{'ohangins'}) && (!exists $options{'multirl'}))
    {
        print "Cannot specify an overhang insert size without the multirl flag\n";
        exec("pod2usage $0");
    }
    return \%options;
}

__DATA__

=head1 NAME

   cutDown_454.pl

=head1 COPYRIGHT

   Copyright (c) 2010 Michael Imelfort
   All rights reserved.

=head1 DESCRIPTION

   Cuts down a (454) shuffled fasta file into shorter reads. Produces at most 2 new files.
   The new cutdown fasta file will be called "FILENAME_cut_xx" where FILENAME is the
   name of the input file and xx is the length of the reads in the output file.
   Each read is cut into a number of smaller reads. To identify the source of the 
   smaller reads the fasta header is modified. If an input read is called:

   >bob_read

   Then the first 4 reads cut from bob will be:

   >1_bob_read
   >2_bob_read
   >3_bob_read
   >4_bob_read

   For 454 type data where there are many different read lengths and where the
   read length of the forward strand is not equal to the readlength of the reverse 
   strand a lot of sequence can go to waste. By using -ohangins parameter you can 
   cut new raired reads from the otherwise wasted overhang. These are always written 
   to the file "FILENAME_overhang_xx_yy" where xx is the read length and yy is the 
   insert size. For the same input read as above, we would get:

   >OH1_bob_read
   >OH2_bob_read
   >OH3_bob_read   etc...

=head1 SYNOPSIS

   cutDown_454.pl -outrl INT -filestem STRING [-extension STRING] [-multirl] [-ohangins INT] [-offset INT] [-directory STRING] [-help] [-man]

      -directory         Directory to place output file into (default: running dir)
      -filestem          The name of the input file (minus extension)
      -extension         The extension of the file (fasta, fna etc. default fasta)
      -offset            The offset between cut reads (default: 1)
      -outrl             The readlength for the output file
      -multirl           Set this flag if the input file has many different readlengths (ie, 454 data) 
      -ohangins          The insert size for when we make overhanging reads (must have -multirl flag)

      -help              Displays basic usage information
      -man               Displays more detailed information

         
=cut

