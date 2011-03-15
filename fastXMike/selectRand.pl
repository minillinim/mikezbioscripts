#!/usr/bin/perl
################################################################################
#
#    selectRand.pl
#    Version: 0.1
#
#    Select illumina reads at random!!!
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

open my $in_fh, "<", $options->{'in'};
my $out_fn = $options->{'in'}.$options->{'num'}.".rand";
open my $out_fh, ">", $out_fn;

# get the total number of reads..
my $exe_string = "grep \">\" ".$options->{'in'}." | wc -l";
my $totes = int(`$exe_string`) / 2;

# work out the cutoff
my $cut = $options->{'num'} / $totes;

my $count = 0;

while(<$in_fh>)
{
  my $header = $_;
  my $seq = <$in_fh>;
  my $header_2 = <$in_fh>;
  my $seq_2 = <$in_fh>;

  if(rand() <= $cut)
  {
    print $out_fh $header;
    print $out_fh $seq;
    print $out_fh $header_2;
    print $out_fh $seq_2;
    $count++;
  }
  last if $count >= $options->{'num'};
}
close $in_fh;
close $out_fh;

######################################################################
# CUSTOM SUBS
######################################################################


######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "in|i:s", "num|n:s" );
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
    
    if(!exists($options{'num'}))
    {
        print "You need to specify the number of reads\n";
        exec("pod2usage $0")
    }

    if(!exists($options{'in'}))
    {
        print "You need to specify an input file\n";
        exec("pod2usage $0")
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

    selectRand.pl

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

   Select reads at random from a shuffled fasta file.

=head1 SYNOPSIS

    selectRand.pl  [-help|h]

      [-in -i]                     Input file (must be shuffled fasta)
      [-num -n]                    Number of reads to select
      [-help -h]                   Displays basic usage information
      
=cut

