#!/usr/bin/perl
use strict;

my ($working_file_name) = @ARGV;
my $out_A_fn = $working_file_name."_A";
my $out_B_fn = $working_file_name."_B";

open my $in_fh, "<", $working_file_name;
open my $out_A_fh, ">", $out_A_fn;
open my $out_B_fh, ">", $out_B_fn;

while(<$in_fh>)
{
  chomp $_;
  my $line_1 = $_;
  my $line_2 = <$in_fh>;
  my $line_3 = <$in_fh>;
  my $line_4 = <$in_fh>;
  my $line_5 = <$in_fh>;
  my $line_6 = <$in_fh>;
  my $line_7 = <$in_fh>;
  my $line_8 = <$in_fh>;

  my @co_ords = split /:/, $line_1;
  my @final_co_ord = split /#/, $co_ords[4];
  print $out_A_fh ">$co_ords[1]:$co_ords[2]:$co_ords[3]:$final_co_ord[0]-1\n$line_2";
  print $out_B_fh ">$co_ords[1]:$co_ords[2]:$co_ords[3]:$final_co_ord[0]-2\n$line_6";
}


close $in_fh;
close $out_A_fh;
close $out_B_fh;
