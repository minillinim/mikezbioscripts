#!/usr/bin/perl -w
###############################################################################
#
#    makeKmerGCCSV.pl  -- does just that
#    Copyright (C) 2009-2012 Lauren Bragg, Michael Imelfort
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


# override defaults
my $out_GC_csv_filename = $options->{'in'}."_GC.csv";
if(exists $options->{'out'})
{
    $out_GC_csv_filename = $options->{'out'}."_GC.csv";
}


# global variables
my @unique_kmer_ids = ();           # ids for unique kmers per kmer size chosen
my @kmer_id_maps = ();              # maps to link the kmers to ids for each kmer size
my @id_kmer_maps = ();              # maps to link the ids to kmers for each kmer size
my @id_counts_maps = ();            # global counts of kmers for each kmer size
my @kmer_lists_per_sequence = ();   # ordered lists of kmer_ids for each sequence seen
my $num_sequences = -1;
my @ids = ();
my @seq_lengths = ();

# specifically for sequence GC count plotting
my $tetramers;
my @gc_at_pos_per_sequence = ();    # the gc vs position values per sequence

# open the input file
my $in_file = IO::File->new($options->{'in'}, "r") or die "can't find: $options->{'in'} $!";
# open the output file
my $GC_file = IO::File->new($out_GC_csv_filename, "w");

my @kmers = ();
my $num_kmers = 0;

# we need to know where the 4mer is hidden
my $kmer_4_index = -1;
if(exists $options->{'kmers'})
{
    @kmers = split /,/, $options->{'kmers'};
    $num_kmers = $#kmers + 1;
    for my $current_kmer (0 .. $#kmers)
    {
        $kmers[$current_kmer] = int($kmers[$current_kmer]);
        $unique_kmer_ids[$current_kmer] = 0;

        if(4 == $kmers[$current_kmer])
        {
            $kmer_4_index = $current_kmer;
        }

        my %tmp_kmer_id_map = ();
        $kmer_id_maps[$current_kmer] = \%tmp_kmer_id_map;

        my %tmp_id_kmer_map = ();
        $id_kmer_maps[$current_kmer] = \%tmp_id_kmer_map;

        my %tmp_id_counts_map = ();
        $id_counts_maps[$current_kmer] = \%tmp_id_counts_map;
    }

    if(-1 == $kmer_4_index)
    {
        die("You need to specify a kmer of length 4 for anything good to ever happen to you!\n");
    }

    $tetramers = makeTetramerGCmap();
}

my $id = "";
my $sequence = "";
my $in_fasta = 0;

if(!exists $options->{'silent'}) { print "Start reading file\n"; }

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
            if(!exists $options->{'silent'}) { print "Processing: $id\n"; }
            push @ids, $id;
            my ($seq_length, $seq_length_with_n, $gc_count) = getGC($sequence);
            my $gc_percentage = ($seq_length > 0)? ($gc_count / $seq_length) : 0;
            print $GC_file $id, ",", $seq_length,",",  $gc_percentage, ",", $seq_length_with_n,  "\n";
            push @seq_lengths, ($seq_length_with_n - 1);

            if(0 != $num_kmers)
            {
                cutMers($sequence);
            }

            $num_sequences++;

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
    if(!exists $options->{'silent'}) { print "Processing Last: $id\n"; }
    push @ids, $id;
    my ($seq_length, $seq_length_with_n, $gc_count) = getGC($sequence);
    my $gc_percentage = ($seq_length > 0)? ($gc_count / $seq_length) : 0;
    print $GC_file $id, ",", $seq_length,",",  $gc_percentage, ",", $seq_length_with_n,  "\n";
    push @seq_lengths, ($seq_length_with_n - 1);
    if(0 != $num_kmers)
    {
        cutMers($sequence);
    }

    $num_sequences++;
}

# Now we've amde all the kmer lists
# convert all the IDs in the sequence specific arrays to the actual counts
for my $current_seq_index (0 .. $num_sequences)
{
    for my $current_kmer (0 .. $#kmers)
    {
        my $this_length_kmer_list = 0;

        # get the kmer list for this kmer (so we can get the length)
        my $current_kmer_list_ref = $kmer_lists_per_sequence[$current_seq_index][$current_kmer]; #$this_id_lists_per_kmer[$current_kmer];
        my @current_kmer_list = @$current_kmer_list_ref;

        for my $j (0 .. $#current_kmer_list)
        {
            $kmer_lists_per_sequence[$current_seq_index][$current_kmer][$j] = $id_counts_maps[$current_kmer]{$kmer_lists_per_sequence[$current_seq_index][$current_kmer][$j]};
            $this_length_kmer_list++;
        }

        if($this_length_kmer_list < $seq_lengths[$current_seq_index])
        {
            for my $fixer ($this_length_kmer_list .. $seq_lengths[$current_seq_index])
            {
                $kmer_lists_per_sequence[$current_seq_index][$current_kmer][$fixer] = 0;
                #$current_kmer_list[$fixer] = 0;
            }
        }
    }
}

if(!exists $options->{'silent'}) { print "Printing\n"; }

# finally we can print to file
for my $current_seq_index (0 .. $num_sequences)
{
    # open the kmer_count file
    my $out_KC_csv_filename = $ids[$current_seq_index]."_KC.csv";
    my $KC_file = IO::File->new($out_KC_csv_filename, "w");

    # print the headers and prep the offsets
    my @buffer_offsets = ();
    print $KC_file "\"position\"";
    for my $current_kmer (0 .. $#kmers)
    {
        print $KC_file ",\"K_".($kmers[$current_kmer])."\"";
        push @buffer_offsets, int($kmers[$current_kmer]/2);
    }
    print $KC_file ",\"GC\"\n";

    for my $position (0 .. $seq_lengths[$current_seq_index])
    {
        print $KC_file "$position";
        for my $current_kmer (0 .. $#kmers)
        {
            # get the kmer list for this kmer
            if($position >= $buffer_offsets[$current_kmer])
            {
                #print $KC_file ",".$this_id_lists_per_kmer[$current_kmer][$position - $buffer_offsets[$current_kmer]];
                print $KC_file ",". $kmer_lists_per_sequence[$current_seq_index][$current_kmer][$position - $buffer_offsets[$current_kmer]];
            }
            else
            {
                print $KC_file ",0";
            }
        }
        if(defined $gc_at_pos_per_sequence[$current_seq_index][$position])
        {
            print $KC_file ",".$gc_at_pos_per_sequence[$current_seq_index][$position];
        }
        else
        {
            print $KC_file ",0";
        }
        print $KC_file "\n";
    }
    $KC_file->close();
}

# close the files
$in_file->close();
$GC_file->close();

###############################################################################
 # SUBS
###############################################################################
sub getGC {
    #
    # get the gc count and length for a sequence
    #
    my ($sequence) = @_;
    my $seq_length = 0;
    my $seq_length_with_n = 0;
    my $gc_count = 0;
    foreach my $base (split(//, $sequence))
    {
        if($base =~ /[CG]/gi)
        {
            $gc_count++;
        }
        if($base !~ /[N]/gi)
        {
            $seq_length++;
        }
        $seq_length_with_n++;
    }
    return ($seq_length, $seq_length_with_n, $gc_count);
}

sub makeTetramerGCmap {

    my %tmer_map = ();
    $tmer_map{"AAAA"} = 0; $tmer_map{"AAAC"} = 1; $tmer_map{"AAAG"} = 1; $tmer_map{"AAAT"} = 0; $tmer_map{"AACA"} = 1; $tmer_map{"AACC"} = 2; $tmer_map{"AACG"} = 2; $tmer_map{"AACT"} = 1;
    $tmer_map{"AAGA"} = 1; $tmer_map{"AAGC"} = 2; $tmer_map{"AAGG"} = 2; $tmer_map{"AAGT"} = 1; $tmer_map{"AATA"} = 0; $tmer_map{"AATC"} = 1; $tmer_map{"AATG"} = 1; $tmer_map{"AATT"} = 0;
    $tmer_map{"ACAA"} = 1; $tmer_map{"ACAC"} = 2; $tmer_map{"ACAG"} = 2; $tmer_map{"ACAT"} = 1; $tmer_map{"ACCA"} = 2; $tmer_map{"ACCC"} = 3; $tmer_map{"ACCG"} = 3; $tmer_map{"ACCT"} = 2;
    $tmer_map{"ACGA"} = 2; $tmer_map{"ACGC"} = 3; $tmer_map{"ACGG"} = 3; $tmer_map{"ACGT"} = 2; $tmer_map{"ACTA"} = 1; $tmer_map{"ACTC"} = 2; $tmer_map{"ACTG"} = 2; $tmer_map{"ACTT"} = 1;
    $tmer_map{"AGAA"} = 1; $tmer_map{"AGAC"} = 2; $tmer_map{"AGAG"} = 2; $tmer_map{"AGAT"} = 1; $tmer_map{"AGCA"} = 2; $tmer_map{"AGCC"} = 3; $tmer_map{"AGCG"} = 3; $tmer_map{"AGCT"} = 2;
    $tmer_map{"AGGA"} = 2; $tmer_map{"AGGC"} = 3; $tmer_map{"AGGG"} = 3; $tmer_map{"AGGT"} = 2; $tmer_map{"AGTA"} = 1; $tmer_map{"AGTC"} = 2; $tmer_map{"AGTG"} = 2; $tmer_map{"AGTT"} = 1;
    $tmer_map{"ATAA"} = 0; $tmer_map{"ATAC"} = 1; $tmer_map{"ATAG"} = 1; $tmer_map{"ATAT"} = 0; $tmer_map{"ATCA"} = 1; $tmer_map{"ATCC"} = 2; $tmer_map{"ATCG"} = 2; $tmer_map{"ATCT"} = 1;
    $tmer_map{"ATGA"} = 1; $tmer_map{"ATGC"} = 2; $tmer_map{"ATGG"} = 2; $tmer_map{"ATGT"} = 1; $tmer_map{"ATTA"} = 0; $tmer_map{"ATTC"} = 1; $tmer_map{"ATTG"} = 1; $tmer_map{"ATTT"} = 0;
    $tmer_map{"CAAA"} = 1; $tmer_map{"CAAC"} = 2; $tmer_map{"CAAG"} = 2; $tmer_map{"CAAT"} = 1; $tmer_map{"CACA"} = 2; $tmer_map{"CACC"} = 3; $tmer_map{"CACG"} = 3; $tmer_map{"CACT"} = 2;
    $tmer_map{"CAGA"} = 2; $tmer_map{"CAGC"} = 3; $tmer_map{"CAGG"} = 3; $tmer_map{"CAGT"} = 2; $tmer_map{"CATA"} = 1; $tmer_map{"CATC"} = 2; $tmer_map{"CATG"} = 2; $tmer_map{"CATT"} = 1;
    $tmer_map{"CCAA"} = 2; $tmer_map{"CCAC"} = 3; $tmer_map{"CCAG"} = 3; $tmer_map{"CCAT"} = 2; $tmer_map{"CCCA"} = 3; $tmer_map{"CCCC"} = 4; $tmer_map{"CCCG"} = 4; $tmer_map{"CCCT"} = 3;
    $tmer_map{"CCGA"} = 3; $tmer_map{"CCGC"} = 4; $tmer_map{"CCGG"} = 4; $tmer_map{"CCGT"} = 3; $tmer_map{"CCTA"} = 2; $tmer_map{"CCTC"} = 3; $tmer_map{"CCTG"} = 3; $tmer_map{"CCTT"} = 2;
    $tmer_map{"CGAA"} = 2; $tmer_map{"CGAC"} = 3; $tmer_map{"CGAG"} = 3; $tmer_map{"CGAT"} = 2; $tmer_map{"CGCA"} = 3; $tmer_map{"CGCC"} = 4; $tmer_map{"CGCG"} = 4; $tmer_map{"CGCT"} = 3;
    $tmer_map{"CGGA"} = 3; $tmer_map{"CGGC"} = 4; $tmer_map{"CGGG"} = 4; $tmer_map{"CGGT"} = 3; $tmer_map{"CGTA"} = 2; $tmer_map{"CGTC"} = 3; $tmer_map{"CGTG"} = 3; $tmer_map{"CGTT"} = 2;
    $tmer_map{"CTAA"} = 1; $tmer_map{"CTAC"} = 2; $tmer_map{"CTAG"} = 2; $tmer_map{"CTAT"} = 1; $tmer_map{"CTCA"} = 2; $tmer_map{"CTCC"} = 3; $tmer_map{"CTCG"} = 3; $tmer_map{"CTCT"} = 2;
    $tmer_map{"CTGA"} = 2; $tmer_map{"CTGC"} = 3; $tmer_map{"CTGG"} = 3; $tmer_map{"CTGT"} = 2; $tmer_map{"CTTA"} = 1; $tmer_map{"CTTC"} = 2; $tmer_map{"CTTG"} = 2; $tmer_map{"CTTT"} = 1;
    $tmer_map{"GAAA"} = 1; $tmer_map{"GAAC"} = 2; $tmer_map{"GAAG"} = 2; $tmer_map{"GAAT"} = 1; $tmer_map{"GACA"} = 2; $tmer_map{"GACC"} = 3; $tmer_map{"GACG"} = 3; $tmer_map{"GACT"} = 2;
    $tmer_map{"GAGA"} = 2; $tmer_map{"GAGC"} = 3; $tmer_map{"GAGG"} = 3; $tmer_map{"GAGT"} = 2; $tmer_map{"GATA"} = 1; $tmer_map{"GATC"} = 2; $tmer_map{"GATG"} = 2; $tmer_map{"GATT"} = 1;
    $tmer_map{"GCAA"} = 2; $tmer_map{"GCAC"} = 3; $tmer_map{"GCAG"} = 3; $tmer_map{"GCAT"} = 2; $tmer_map{"GCCA"} = 3; $tmer_map{"GCCC"} = 4; $tmer_map{"GCCG"} = 4; $tmer_map{"GCCT"} = 3;
    $tmer_map{"GCGA"} = 3; $tmer_map{"GCGC"} = 4; $tmer_map{"GCGG"} = 4; $tmer_map{"GCGT"} = 3; $tmer_map{"GCTA"} = 2; $tmer_map{"GCTC"} = 3; $tmer_map{"GCTG"} = 3; $tmer_map{"GCTT"} = 2;
    $tmer_map{"GGAA"} = 2; $tmer_map{"GGAC"} = 3; $tmer_map{"GGAG"} = 3; $tmer_map{"GGAT"} = 2; $tmer_map{"GGCA"} = 3; $tmer_map{"GGCC"} = 4; $tmer_map{"GGCG"} = 4; $tmer_map{"GGCT"} = 3;
    $tmer_map{"GGGA"} = 3; $tmer_map{"GGGC"} = 4; $tmer_map{"GGGG"} = 4; $tmer_map{"GGGT"} = 3; $tmer_map{"GGTA"} = 2; $tmer_map{"GGTC"} = 3; $tmer_map{"GGTG"} = 3; $tmer_map{"GGTT"} = 2;
    $tmer_map{"GTAA"} = 1; $tmer_map{"GTAC"} = 2; $tmer_map{"GTAG"} = 2; $tmer_map{"GTAT"} = 1; $tmer_map{"GTCA"} = 2; $tmer_map{"GTCC"} = 3; $tmer_map{"GTCG"} = 3; $tmer_map{"GTCT"} = 2;
    $tmer_map{"GTGA"} = 2; $tmer_map{"GTGC"} = 3; $tmer_map{"GTGG"} = 3; $tmer_map{"GTGT"} = 2; $tmer_map{"GTTA"} = 1; $tmer_map{"GTTC"} = 2; $tmer_map{"GTTG"} = 2; $tmer_map{"GTTT"} = 1;
    $tmer_map{"TAAA"} = 0; $tmer_map{"TAAC"} = 1; $tmer_map{"TAAG"} = 1; $tmer_map{"TAAT"} = 0; $tmer_map{"TACA"} = 1; $tmer_map{"TACC"} = 2; $tmer_map{"TACG"} = 2; $tmer_map{"TACT"} = 1;
    $tmer_map{"TAGA"} = 1; $tmer_map{"TAGC"} = 2; $tmer_map{"TAGG"} = 2; $tmer_map{"TAGT"} = 1; $tmer_map{"TATA"} = 0; $tmer_map{"TATC"} = 1; $tmer_map{"TATG"} = 1; $tmer_map{"TATT"} = 0;
    $tmer_map{"TCAA"} = 1; $tmer_map{"TCAC"} = 2; $tmer_map{"TCAG"} = 2; $tmer_map{"TCAT"} = 1; $tmer_map{"TCCA"} = 2; $tmer_map{"TCCC"} = 3; $tmer_map{"TCCG"} = 3; $tmer_map{"TCCT"} = 2;
    $tmer_map{"TCGA"} = 2; $tmer_map{"TCGC"} = 3; $tmer_map{"TCGG"} = 3; $tmer_map{"TCGT"} = 2; $tmer_map{"TCTA"} = 1; $tmer_map{"TCTC"} = 2; $tmer_map{"TCTG"} = 2; $tmer_map{"TCTT"} = 1;
    $tmer_map{"TGAA"} = 1; $tmer_map{"TGAC"} = 2; $tmer_map{"TGAG"} = 2; $tmer_map{"TGAT"} = 1; $tmer_map{"TGCA"} = 2; $tmer_map{"TGCC"} = 3; $tmer_map{"TGCG"} = 3; $tmer_map{"TGCT"} = 2;
    $tmer_map{"TGGA"} = 2; $tmer_map{"TGGC"} = 3; $tmer_map{"TGGG"} = 3; $tmer_map{"TGGT"} = 2; $tmer_map{"TGTA"} = 1; $tmer_map{"TGTC"} = 2; $tmer_map{"TGTG"} = 2; $tmer_map{"TGTT"} = 1;
    $tmer_map{"TTAA"} = 0; $tmer_map{"TTAC"} = 1; $tmer_map{"TTAG"} = 1; $tmer_map{"TTAT"} = 0; $tmer_map{"TTCA"} = 1; $tmer_map{"TTCC"} = 2; $tmer_map{"TTCG"} = 2; $tmer_map{"TTCT"} = 1;
    $tmer_map{"TTGA"} = 1; $tmer_map{"TTGC"} = 2; $tmer_map{"TTGG"} = 2; $tmer_map{"TTGT"} = 1; $tmer_map{"TTTA"} = 0; $tmer_map{"TTTC"} = 1; $tmer_map{"TTTG"} = 1; $tmer_map{"TTTT"} = 0;
    return \%tmer_map;
}

sub cutMers {
    #
    # tap into them global vars!
    #
    my ($sequence) = @_;
    my $seq_length = length $sequence;
    my @current_id_list_per_mer = ();
    my $id;
    for my $current_kmer (0 .. $#kmers)
    {
        # vars needed for this loop
        my @current_id_list = ();
        my @gc_list = ();
        my $current_mer = $kmers[$current_kmer];
        my $total_chomps = $seq_length - $current_mer;

        my $start_chomp = 0;
        while($start_chomp < $total_chomps)
        {
            # get the current kmer
            my $this_mer = substr $sequence, $start_chomp, $current_mer;

            if(4 == $current_mer)
            {
                push @gc_list, $tetramers->{$this_mer};
            }

            #check if it's in the global array
            if(exists $kmer_id_maps[$current_kmer]{$this_mer})
            {
                # we've seen it before
                $id = $kmer_id_maps[$current_kmer]{$this_mer};
                $id_counts_maps[$current_kmer]{$id}++;
            }
            else
            {
                # make a new id
                $unique_kmer_ids[$current_kmer]++;
                $id = $unique_kmer_ids[$current_kmer];

                $id_kmer_maps[$current_kmer]{$id} = $this_mer;
                $kmer_id_maps[$current_kmer]{$this_mer} = $id;

                $id_counts_maps[$current_kmer]{$id} = 1;
            }
            push @current_id_list, $id;
            $start_chomp++;
        }
        push @current_id_list_per_mer, \@current_id_list;

        # if the current kmer is 4 push the gc count back on the list
        if(4 == $current_mer)
        {
            push @gc_at_pos_per_sequence, \@gc_list;
        }
    }
    push @kmer_lists_per_sequence, \@current_id_list_per_mer;
}

sub checkParams {
    my @standard_options = ( "help+", "in:s", "out:s", "kmers:s", "silent+" );
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

    return \%options;
}


sub printAtStart {
print<<"EOF";
----------------------------------------------------------------
 $0
 Copyright (C) 2009 - 2012 Lauren Bragg, Michael Imelfort

 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
----------------------------------------------------------------
EOF
}

__DATA__

=head1 NAME

    makeKmerGCCSV.pl

=head1 COPYRIGHT

   copyright (C) 2009 - 2012 Lauren Bragg, Michael Imelfort

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

  makeKmerGCCSV.pl
  Given a multiple fasta file, calculate the length of the sequences, and
  also the gc_content. Also calculate k-mer coverages for a given set of kmers.
  Output everything to a CSV file

=head1 SYNOPSIS

   makeKmerGCCSV.pl -in fasta_file [-kmers kmer1[,kmer2,..]] [-out output_filestem] [-help] [-silent]

    -in fasta_file            The file to work on
    -kmers kmer1[,kmer2,...]  Comma separated list of kmers to work out coverage for
    -out output_file          Output filestem [default: input filename]
    [-help]                   Displays basic usage information
    [-silent]                 Suppress all output to the display

=cut


