#!/usr/bin/perl
use strict;
 # HWI-EAS68_8259_FC30BPV_PE:1:1:1407:416#CTCACGC/1:AAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGC:hhhhhhh[hhhhhhhhhhhhhhhhhhhhh[hNhhhhhhhhhhhhhAhhAh
my ($working_file_name, $suffix_read, $suffix_file) = @ARGV;
my @file_root = split /\./, $working_file_name;
my $out_fn = $file_root[0].$suffix_file.".fasta";

open my $in_fh, "<", $working_file_name;
open my $out_fh, ">", $out_fn;

while(<$in_fh>)
{
  my $line = $_;
  chomp $line;

  my @line_fields = split /:/, $line;
  my @final_co_ord = split /#/, $line_fields[4];
  print $out_fh ">$line_fields[1]:$line_fields[2]:$line_fields[3]:$final_co_ord[0]$suffix_read\n$line_fields[5]\n";
}


close $in_fh;
close $out_fh;
