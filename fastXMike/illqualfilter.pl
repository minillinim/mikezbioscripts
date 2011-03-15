#!/usr/bin/perl
################################################################################
#
#    File: illqualfilter.pl
#    Version: 0.1
#
#    Take raw ill files (dual fastq) and produce a filtered 
#    shuffled fasta file as output
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

# override defaults
my $trim  = 0;
my $reject_CO = 0;
my $trim_chars = "\@ABCDEFGHIJ";
if(exists $options->{'trimmers'}) { $trim_chars = $options->{'trimmers'}; }
print "\n Removing all reads with \"Ns\" and ";
if(exists $options->{'trim'}) { $trim = 1; print "trimming "; $reject_CO = $options->{'trim'}; }
print "any reads with quality\n strings containing the following: $trim_chars\n\n";
if($reject_CO < 0) { $reject_CO = 0; }
# 4 stats
my $num_removed_qual = 0;       # removed due to low quality score
my $num_removed_seq = 0;        # removed becuase it contained an N
my $num_changed = 0;            # number changed but not removed
my $num_processed = 0;          # total number processed

# do stuff
# open all files
my $in1fh;
my $in2fh;
my $outfh;
open $in1fh, "<", $options->{'in1'} or die "could not open file \"$options->{'in1'}\" $!\n";
my $use_2 = 1;
if(!($options->{'in2'} eq "0")) {
    open $in2fh, "<", $options->{'in2'} or die "could not open file \"$options->{'in2'}\" $!\n";
}
else
{
    $use_2 = 0;
}
open $outfh, ">", $options->{'out'} or die "could not open file \"$options->{'out'}\" $!\n";

while(<$in1fh>)
{
    my $f1_l1 = $_;
    my $f1_l2 = <$in1fh>;
    my $f1_l3 = <$in1fh>;
    my $f1_l4 = <$in1fh>;

    my $f2_l1;
    my $f2_l2;
    my $f2_l3;
    my $f2_l4;

    if(1 == $use_2)
    {
        $f2_l1 = <$in2fh>;
        $f2_l2 = <$in2fh>;
        $f2_l3 = <$in2fh>;
        $f2_l4 = <$in2fh>;
    }
    else
    {
        $f2_l1 = <$in1fh>;
        $f2_l2 = <$in1fh>;
        $f2_l3 = <$in1fh>;
        $f2_l4 = <$in1fh>;
    }
    $num_processed++;
    
    # check out if the sequence contains an 'N' or if the quality score contains one of the trimmer chars
    if($f1_l2 =~ /N/ || $f2_l2 =~ /N/ ) {
        $num_removed_seq++; 
        next;
    }
    
    if($f1_l2 =~ /^[C]*$/ || $f2_l2 =~ /^[C]*$/ ) {
        $num_removed_seq++; 
        next;
    }
    
    my $trim_res = 0;
    if($f1_l4 =~ /[$trim_chars]/ || $f2_l4 =~ /[$trim_chars]/) {
        # check to see if we should just trim this guy and not remove
        if(1 == $trim) { 
            $trim_res = trim_smart(\$f1_l2, $f1_l4, \$f2_l2, $f2_l4, $reject_CO); 
            if(0 == $trim_res) { $num_changed++; }
            else { $num_removed_qual++; next; }
        }
        else { $num_removed_qual++;  next;}
    }
    # once we get here we know all is good, but check the value of trim_res.
    $f1_l1 =~ s/@/>/;
    $f2_l1 =~ s/@/>/;
    print $outfh $f1_l1;
    print $outfh $f1_l2;
    print $outfh $f2_l1;
    print $outfh $f2_l2;
}

print "----------------------------------------------------------------\n";
print "  Report: $options->{'out'}\n";
print "----------------------------------------------------------------\n"; 
print "  Processed:\t$num_processed\n";
print "  Removed N:\t$num_removed_seq\t(".($num_removed_seq/$num_processed*100)."%)\n";
print "  Removed Q:\t$num_removed_qual\t(".($num_removed_qual/$num_processed*100)."%)\n";
print "  Trimmed Q:\t$num_changed\t(".($num_changed/$num_processed*100)."%)\n";
print "  Remaining:\t".($num_processed-$num_removed_seq-$num_removed_qual)."\t(".(($num_processed-$num_removed_seq-$num_removed_qual)/$num_processed*100)."%)\n";
print "----------------------------------------------------------------\n";

close $in1fh;
if(1 == $use_2 ) { close $in2fh; }
close $outfh;
###################################################################################
 # SUBS
###################################################################################
sub trim_smart {
    # find the first instance of a trimchar and trim the sequences accordingly.
    my ($r1_rf, $q1, $r2_rf, $q2, $reject) = @_;
    my $sub_pos = 0;
    foreach my $qual (split //, $q1)
    {
        if($qual =~ /[$trim_chars]/)
        {
            last;
        }
        $sub_pos++;
    }
    if($sub_pos > $reject)
    {
        chomp $$r1_rf;
        $$r1_rf = substr $$r1_rf, 0, $sub_pos;
        $$r1_rf .=" \n";
    }
    else { return -1; }
        
    $sub_pos = 0;
    foreach my $qual (split //, $q2)
    {
        if($qual =~ /[$trim_chars]/)
        {
            last;
        }
        $sub_pos++;
    }
    if($sub_pos > $reject)
    {
        chomp $$r2_rf;
        $$r2_rf = substr $$r2_rf, 0, $sub_pos;
        $$r2_rf .=" \n";
    }
    else { return -1; }
    
    return 0;
}

sub checkParams {
    my @standard_options = ( "help+", "in1:s", "in2:s", "out:s", "trimmers:s", "trim:i" );
    my %options;

    # Add any other command line options, and the code to handle them
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    if(!exists $options{'in1'} || !exists $options{'in2'} || !exists $options{'out'})
    {
        exec("pod2usage $0"); 
        die;
    }

    return \%options;
}


sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2010, 2011 Michael Imelfort
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    illqualfilter.pl

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

   Take raw ill files (dual fastq) and produce a filtered shuffled fasta 
   file as output. very basic

=head1 SYNOPSIS

    illqualfilter.pl -in1 FILE1 -in2 FILE2 -out SHUFF [-trimmers String] [-trim REJECT] [-help]
 
    -in1 FILE1          First fastq file
    -in2 FILE2          Second fastq file or '0' if the first file is a shuffled fastq
    -out SHUFF          Output shuffled fasta file
    [-trimmers String]  List of quality scores to trim away [default : "@ABCDEFGHIJ"]
    [-trim REJECT]      Remove bad bases - not the whole read [default : remove whole] 
                        REJECT any trimmed reads shorter than this length, 0 for ignore
    [-help]             Displays basic usage information
         
=cut
