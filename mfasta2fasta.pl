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
  my $line_1 = $_;
  my $line_2 = <$in_fh>;
  my $line_3 = <$in_fh>;
  my $line_4 = <$in_fh>;

  print $out_A_fh "$line_1$line_2";
  print $out_B_fh "$line_3$line_4";
}

close $in_fh;
close $out_A_fh;
close $out_B_fh;
