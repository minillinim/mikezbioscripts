#!/usr/bin/perl
###############################################################################
#
#    remove_exact_dupes.pl
#    Version: 0.1
#    
#    Looks at a shuffled FASTA file and removes any exact duplicate pairs
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

######################################################################
# CODE HERE
######################################################################
my $num_removed = 0;            # removed pairs
my $num_processed = 0;          # total number processed
my %read_keys = ();             # map used to determine dupes...

# open all files
my $infh;
my $outfh;
open $infh, "<", $options->{'in'} or die "could not open file \"$options->{'in'}\" $!\n";
open $outfh, ">", $options->{'out'} or die "could not open file \"$options->{'out'}\" $!\n";

while(<$infh>)
{
    my $f1_l1 = $_;
    my $f1_l2 = <$infh>;
    my $f1_l3 = <$infh>;
    my $f1_l4 = <$infh>;
    
    chomp $f1_l2;
    chomp $f1_l4;
    
    $num_processed++;
    
    my $key_string = makeKey($f1_l2, $f1_l4);
    
    if(!exists $read_keys{$key_string})
    {
        $read_keys{$key_string} = 1;
        print $outfh $f1_l1;
        print $outfh $f1_l2."\n";
        print $outfh $f1_l3;
        print $outfh $f1_l4."\n";
    }
    else
    {
        $num_removed++;
    }
}

print "----------------------------------------------------------------\n";
print "  Report: $options->{'out'}\n";
print "----------------------------------------------------------------\n";
print "  Processed:\t$num_processed\n";
print "  Removed:\t$num_removed\t(".(($num_processed-$num_removed)/$num_processed*100).")\n";
print "  Remaining:\t".($num_processed-$num_removed)."\n";
print "----------------------------------------------------------------\n";

close $infh;
close $outfh;


######################################################################
# CUSTOM SUBS
######################################################################
sub makeKey {
    my ($seq_1, $seq_2) = @_;
    my $key_str;
    
    my $llc_seq_1 = llc($seq_1);
    my $llc_seq_2 = llc($seq_2);
    
    if($llc_seq_1 lt $llc_seq_2) { $key_str = $seq_1.$seq_2; }
    else { $key_str =  revcompl($seq_2). revcompl($seq_1); }
    return $key_str;
}

sub llc {
    my ($query) = @_;
    my $rev_q = revcompl($query);
    if($rev_q lt $query) { return $rev_q; }
    return $query;
}

sub revcompl {
    my ($rev) = @_;
    $rev =~ tr/ACGTacgt/TGCAtgca/; 
    return reverse $rev;
}
    
######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help+", "in:s", "out:s" );
    my %options;

    # Add any other command line options, and the code to handle them
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    if(!exists $options{'in'} || !exists $options{'out'})
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
 Copyright (C) 2011 Michael Imelfort
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    remove_exact_dupes.pl

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

   Looks at a shuffled FASTA file and removes any exact duplicate pairs

=head1 SYNOPSIS

    remove_exact_dupes.pl -in INFILE -out OUTFILE [-help]

      -in INFILE         Shuffled fasta file to search
      -out OUTFILE       File to output unique reads to
      [-help]            Displays basic usage information
         
=cut

