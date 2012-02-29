#!/usr/bin/perl
###############################################################################
#
#    velvet_fixer.pl
#    Reverses the second read in each pair so that velvet doesn't freak out
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

# open all files
my $infh;
my $outfh;
open $infh, "<", $options->{'in'} or die "could not open file \"$options->{'in'}\" $!\n";
open $outfh, ">", $options->{'out'} or die "could not open file \"$options->{'out'}\" $!\n";

while(<$infh>)
{
    my $tmp = $_;
    print $outfh $tmp;
    $tmp = <$infh>;
    print $outfh $tmp;
    $tmp = <$infh>;
    print $outfh $tmp;
    $tmp = <$infh>;
    chomp $tmp;
    print $outfh revcompl($tmp)."\n";
}

close $infh;
close $outfh;


######################################################################
# CUSTOM SUBS
######################################################################

sub revcompl {
    my ($rev) = @_;
    $rev =~ tr/ACGTacgt/TGCAtgca/; 
    return reverse scalar $rev;
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

    velvet_fixer.pl

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

   Reverses the second read in each pair so that velvet doesn't freak out

=head1 SYNOPSIS

    velvet_fixer.pl -in INFILE -out OUTFILE [-help]

      -in INFILE         Shuffled fasta file to search
      -out OUTFILE       File to output unique reads to
      [-help]            Displays basic usage information
         
=cut

