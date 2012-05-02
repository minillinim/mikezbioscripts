#!/usr/bin/env perl
###############################################################################
#
#    getGT.pl
#    
#    Select information and sequences from files based on sequence lengths
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
use Data::Dumper;

#CPAN modules
use Bio::SeqIO;

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
printAtStart();
my $global_options = checkParams();

######################################################################
# CODE HERE
######################################################################
my $global_prefix_tag = overrideDefault("GT_".$global_options->{'cutoff'}, 'tag');

# for user stats
my $global_parsed = 0;
my $global_kept = 0;

# see if there are any aux files to open
my @aux_in_fhs = ();
my @aux_out_fhs = ();
my @aux_headers = ();
my @aux_seps = ();

# to hold all of the data!
my @aux_data = ();

if(exists $global_options->{'aux'})
{
    my @aux_fns = split(/,/, $global_options->{'aux'});
    foreach my $aux_fn (@aux_fns)
    {
        push @aux_in_fhs, openRead($aux_fn);
        my $aux_out_fn = $global_prefix_tag."_".$aux_fn;
        push @aux_out_fhs, openWrite($aux_out_fn);
        
        # make a hash to hold the data
        my %tmp_hash = ();
        push @aux_data, \%tmp_hash;
    }
    
    # check and chomp headers!
    if(exists $global_options->{'auxh'})
    {
        my $aux_counter = 0;
        my @auxhs = split(/,/, $global_options->{'auxh'});
        if($#auxhs != $#aux_fns) { croak "**ERROR: $0 : AUX files and AUX headers must be the same length!\n";} 
        foreach my $auxh (@auxhs)
        {
            # if it has a header, then chomp it off and print it into the file
            if(1 == $auxh)
            {
                my $in_fh = $aux_in_fhs[$aux_counter];
                my $out_fh = $aux_out_fhs[$aux_counter];
                my $header = <$in_fh>;
                print $out_fh $header;
            }
            $aux_counter++;
        }
    }
    
    # get the separators in order
    if(exists $global_options->{'auxs'})
    {
        my @auxss = split(/,/, $global_options->{'auxs'});
        foreach my $auxs (@auxss)
        {
            if('c' eq $auxs) { push @aux_seps, ","; }
            elsif ('t' eq $auxs) { push @aux_seps, "\t"; }
            else { croak "**ERROR: $0 : Unknown separator \"$auxs\"\n"; }
        }
    }
    else
    {
        # all commas
        foreach my $aux_fn (@aux_fns)
        {
            push @aux_seps, ",";
        }
    }
    
    # now parse in all the aux info
    my $aux_counter = 0;
    foreach my $aux_fh (@aux_in_fhs)
    {
        my $aux_data_ref = $aux_data[$aux_counter];
        while(<$aux_fh>)
        {
            # make the "key => value" be contig_header => csv_line  
            my @fields = split(/$aux_seps[$aux_counter]/, $_);
            ${$aux_data_ref}{$fields[0]} = $_;
        }
        $aux_counter++;
    }    

    # close the aux in fhs
    foreach my $aux_fh (@aux_in_fhs)
    {
        close($aux_fh);
    }
}

# parse the fasta file
my $seqio_in = Bio::SeqIO->new(-file => $global_options->{'fasta'}, '-format' => 'Fasta');
my $seqio_out = Bio::SeqIO->new(-file => ">$global_prefix_tag"."_".$global_options->{'fasta'}, '-format' => 'Fasta');
while(my $seq = $seqio_in->next_seq)
{ 
    $global_parsed++;
    my $seq_string = $seq->seq;
    my $seq_id = $seq->id;
    
    next if(length $seq_string < $global_options->{'cutoff'});
    $global_kept++;
    
    # this sequence passed the test!
    $seqio_out->write_seq($seq);
    
    # now check out the aux data
    my $aux_counter = 0;
    while($aux_counter <= $#aux_data)
    {
        my $aux_data_ref = $aux_data[$aux_counter];
        if(exists ${$aux_data_ref}{$seq_id})
        {
            # we have an entry for this contig!
            # get the file handle and print!
            my $out_fh = $aux_out_fhs[$aux_counter];
            print $out_fh ${$aux_data_ref}{$seq_id};
        }
        $aux_counter++;
    }   
}

# close the aux out fhs
foreach my $aux_fh (@aux_out_fhs)
{
    close($aux_fh);
}

print "Parsed: $global_parsed\n";
print "Retained: $global_kept\n";

######################################################################
# CUSTOM SUBS
######################################################################

######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "cutoff|c:i", "fasta|f:s", "aux|a:s", "auxh|h:s", "auxs|s:s", "tag|t:s", "help|h+");
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
    if(!exists $options{'cutoff'} ) { printParamError ("Need a cut off!"); }
    if(!exists $options{'fasta'} ) { printParamError ("Need a file to parse!"); }
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

    getGT.pl

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

   Select information and sequences from files based on sequence lengths

=head1 SYNOPSIS

    getGT.pl -cutoff|c INT -fasta|f FILE
    
      -cutoff -c INT               Cutoff to select above (inclusive)
      -fasta -f FILE               Fasta file to parse
      [-tag -t STRING]             Tag to attach to processed files [default: GT_cutoff]
      [-aux -a FILE[,FILE...]]     Auxilary files to parse
      [-auxh -h INT[,INT...]]      Header flags for aux files
      [-auxs -s CHAR[,CHAR]...]]   The separators used in the aux files [default ',']
                                   Use 'c' == COMMA, 't' = TAB
      [-help -h]                   Displays basic usage information
         
         
    EX:
       
       getGT.pl -c 1000 -f bob.fa -a cov1.csv,gt.csv -h 0,1 -s c,t 
       
       This will select all seqs from bob.fa which are AT LEAST 1000bp long
       AND it will select all the relevant rows from cov1.csv and gc.csv. 
       cov1.csv has no header line and is comma separated.
       gc.csv has a header and is tab separated. 
       
       If you leave out auxh, we'll assume no headers!
       If you leave out auxs, we'll assum commas!
       
       This script assumes that the first column in the csv files contains the
       contig header. There should be no holes.
       
=cut
