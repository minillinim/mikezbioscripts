#!/usr/bin/perl
###############################################################################
#
#    sammy2.pl
#    
#    <one line to give the program's name and a brief idea of what it does.>
#
#    Copyright (C) 2011 Michael Imelfort
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

#CPAN modules

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params and print copyright
printAtStart();
my $options = checkParams();

######################################################################
# CODE HERE
######################################################################

    # open files
    open my $output_fh, ">", $options->{'outfile'} or die $!;
    open my $input_fh, "<", $options->{'infile'} or die $!;
    
    # get input sequence
    my $global_sequence_header = <$input_fh>;
    my @sequence_lines = <$input_fh>;
    my $global_sequence = join "", @sequence_lines;
    $global_sequence =~ s/\n//g;
    my $global_rev_sequence = revcompl($global_sequence);
    my $global_source_length = length $global_sequence;
    
    close $input_fh;
    
    # label these guys globals
    my $global_read_length = $options->{'readlength'};
    my $global_read_counter = 1;
    my $global_orientation = $options->{'orientation_type'}; 
    
    my %global_ins_values = ();
    my $global_ins_upper_rand = 0;
    my $global_prob_resolution = 10000;
    my $global_inv_prob_res = 0.0001;
    
    # for use when we have a database
    my %global_db_reads = ();
    
    # prep the array we'll use to hold inset sizes
    prep_insert();

    my $num_pairs = $options->{'numpairs'};
    if(exists $options->{'database'})
    {
        # parse in the database
        open my $db_fh, "<", $options->{'database'} or die $!;
        while(<$db_fh>)
        {
            chomp $_;
            $global_db_reads{$_} = 1;
        }
        close $db_fh;
        while($num_pairs > 0)
        {
            my $r1 = "NULL";
            my $r2 = "NULL";
            while (1){
                ($r1, $r2) = cut_reads(get_insert());
                if((exists $global_db_reads{lowlexi($r1)}) and (exists $global_db_reads{lowlexi($r2)})) 
                { last; }
            }
            print $output_fh ">Read$global_read_counter"."_1\n$r1\n>Read$global_read_counter"."_2\n$r2\n";
            $num_pairs--;
            $global_read_counter++;
        }
        
    }
    else
    {
        while($num_pairs > 0)
        {
            my ($r1, $r2) = cut_reads(get_insert());
            print ">Read$global_read_counter"."_1\n$r1\n>Read$global_read_counter"."_2\n$r2\n";
            $num_pairs--;
            $global_read_counter++;
        }
    }

    # clean up
    close $output_fh;

######################################################################
# CUSTOM SUBS
######################################################################

sub get_insert {
    # Returns an random integer between [ mean insert - (3 * stdev) --> mean insert + (3 * stdev) ]    
    # if run repeatedly, the integers returned will fit a normal distribution
    # this is not 100% true of short read data but it will do    
    # if you're not happy, write your own version here...
    my $rand_num = 1;
    while($rand_num > $global_ins_upper_rand)
    {  
        $rand_num = int(rand() * $global_prob_resolution) / $global_prob_resolution;
    } 
    return $global_ins_values{$rand_num};
}

sub prep_insert {
# prep the arrays we'll need to use in get insert

    # these are the limits for the insert sizes
    my $upper_insert_limit = $options->{'insert'} + 3 * $options->{'deviation'};
    my $lower_insert_limit = $options->{'insert'} - 3 * $options->{'deviation'};

    # need these guys to make a normal distribution
    my $pi = 4 * atan2(1, 1);
    my $cumulative_area = 0;
    my $multiplier = 1 / (int($options->{'deviation'}) * sqrt(2 * $pi));
    my $two_sig_squared =  int($options->{'deviation'}) * int($options->{'deviation'}) * 2;
    my $step = 0.001;
    
    my $x  = 0;
    my $key = $lower_insert_limit; 
    my $prev_cum_area = 0;
    
# do the first half notch 
    for($x = $lower_insert_limit; $x < $lower_insert_limit + 0.5; $x = $x + $step)
    {
        my $height = $multiplier * exp(-1 * ($x - int($options->{'insert'})) * ($x - int($options->{'insert'})) / $two_sig_squared);
        $cumulative_area = $cumulative_area + $height * $step; 
    }
# update the hash
    my $new_cum_area = int($cumulative_area * $global_prob_resolution);
    for(my $counter = $new_cum_area; $counter >= 0; $counter--)
    {
        $global_ins_values{($counter / $global_prob_resolution)} = $key;
    }
    $prev_cum_area = $new_cum_area;
    $key++; 
    
# do all the middle ones
    $lower_insert_limit = $lower_insert_limit + 0.5;
    while($x < $upper_insert_limit - 0.5)
    {  
        for($x = $lower_insert_limit; $x < $lower_insert_limit + 1; $x = $x + $step)
        {
            my $height = $multiplier * exp(-1 * ($x - int($options->{'insert'})) * ($x - int($options->{'insert'})) / $two_sig_squared);
            $cumulative_area = $cumulative_area + $height * $step; 
        }
        $lower_insert_limit++;
        
        $new_cum_area = int($cumulative_area * $global_prob_resolution);
        for(my $counter = $new_cum_area; $counter > $prev_cum_area; $counter--)
        {
            $global_ins_values{($counter / $global_prob_resolution)} = $key;
        }
        $prev_cum_area = $new_cum_area;
        $key++;
        $x = $lower_insert_limit; 
    }
# do the final half notch 
    for($x = $lower_insert_limit; $x < $upper_insert_limit; $x = $x + $step)
    {
        my $height = $multiplier * exp(-1 * ($x - int($options->{'insert'})) * ($x - int($options->{'insert'})) / $two_sig_squared);
        $cumulative_area = $cumulative_area + $height * $step; 
    }
    $new_cum_area = int($cumulative_area * $global_prob_resolution);
    for(my $counter = $new_cum_area; $counter >= $prev_cum_area; $counter--)
    {
         $global_ins_values{($counter / $global_prob_resolution)} = $upper_insert_limit;
    }
    $global_ins_upper_rand = $cumulative_area;
}

sub cut_reads {
    # actually cut the read (and add errors) based on template
    # information
    # return reads
    my($insert_size) = @_;
    my $read_1 = "";
    my $read_2 = "";
    
    # choose a random start site
    my $start = int(rand($global_source_length - $insert_size - $global_read_length));
    
    # choose a random strand
    my $strand = int(rand(1.99999));
    
    if($global_orientation == 0) {
        # <-- -->
        if($strand == 0) { $read_1 = revcompl(substr($global_sequence, $start, $global_read_length));  $read_2 = substr($global_sequence, $start + $insert_size - $global_read_length, $global_read_length); }
        else { $read_1 = revcompl(substr($global_rev_sequence, $start, $global_read_length));  $read_2 = substr($global_rev_sequence, $start + $insert_size - $global_read_length, $global_read_length); }
    }
    elsif($global_orientation == 2) {
        # --> <--
        if($strand == 0) { $read_1 = substr($global_sequence, $start, $global_read_length);  $read_2 = revcompl(substr($global_sequence, $start + $insert_size - $global_read_length, $global_read_length)); }
        else { $read_1 = substr($global_rev_sequence, $start, $global_read_length);  $read_2 = revcompl(substr($global_rev_sequence, $start + $insert_size - $global_read_length, $global_read_length)); }
    }
    else { # if($global_orientation == 1) || ($global_orientation == 3)
        # <-- <-- or --> -->
        if($strand == 0) { $read_1 = revcompl(substr($global_sequence, $start, $global_read_length));  $read_2 = revcompl(substr($global_sequence, $start + $insert_size - $global_read_length, $global_read_length)); }
        else { $read_1 = revcompl(substr($global_rev_sequence, $start, $global_read_length));  $read_2 = revcompl(substr($global_rev_sequence, $start + $insert_size - $global_read_length, $global_read_length)); }
    }
    return ($read_1, $read_2);
}

sub revcompl {
    my ($seq) = @_;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return reverse $seq;
}

sub lowlexi {
    my ($seq) = @_;
    my $rev_seq = revcompl($seq);
    if($seq lt $rev_seq) { return $seq; }
    else { return $rev_seq; }
}

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "verbose|v+", "readlength|l:i", "insert|i:i", "deviation|d:i", "numpairs|n:i", "infile|f:s", "outfile|o:s", "orientation_type|t:i", "database|b:s" );
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

    if (!exists $options{'insert'}) 
    {
        print "***ERROR: Missing insert\n";
        exec("pod2usage $0");        
    }
    
    if (!exists $options{'deviation'})
    {
        print "***ERROR: Missing deviation\n";
        exec("pod2usage $0");        
    }
    
    if (!exists $options{'infile'}) 
    {
        print "***ERROR: Missing infile\n";
        exec("pod2usage $0");        
    }

    if (!exists $options{'outfile'}) 
    {
        print "***ERROR: Missing outfile\n";
        exec("pod2usage $0");        
    }

    if (!exists $options{'readlength'})
    {
        print "***ERROR: Missing readlength\n";
        exec("pod2usage $0");        
    }
    
    if (!exists $options{'orientation_type'})
    {
        print "***ERROR: Missing orientation_type\n";
        exec("pod2usage $0");        
    }
    else
    {
        if(!($options{'orientation_type'} == 0 || $options{'orientation_type'} == 1 || $options{'orientation_type'} == 2 || $options{'orientation_type'} == 3))
        {
            print "***ERROR: Invalid orientation_type. Need [0,1,2,3]\n";
            exec("pod2usage $0");        
        }
    }
    
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

    sammy2.pl

=head1 COPYRIGHT

   copyright (C) 2011 Michael Imelfort

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

   Insert detailed description here

=head1 SYNOPSIS

    sammy2.pl readlength|l INT insert|i INT deviation|d INT numpairs|n INT infile|f STRING outfile|o STRING orientation_type|t INT database|b STRING [-verbose|v] [-help|h]

        --readlength -l        INT         The length of the reads
        --insert -i            INT         The insert size
        --deviation -d         INT         The standard deviation
        --numpairs -n          INT         The number of pairs
        --infile -f            STRING      Input file to cut reads from
        --outfile -o           STRING      Output file
        --orientation_type -t  INT         Orientation type [0,1,2,3]
        [--database -b]        STRING      Name of database to use for read checking
        [--verbose -v]                     Be verbose...
        [--help -h]                        Displays basic usage information
         
=cut

