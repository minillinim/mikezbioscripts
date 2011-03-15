#!/usr/bin/perl
use strict;

my ($in_fn) = @ARGV;

open my $in_fh, "<", $in_fn;
my $out_fn_1 = $in_fn."_A";
my $out_fn_2 = $in_fn."_B";
open my $out_fh_1, ">", $out_fn_1;
open my $out_fh_2, ">", $out_fn_2;

while(<$in_fh>)
{
  my $line_1 = $_;
  my $line_2 = <$in_fh>;
  my $line_3 = <$in_fh>;
  my $line_4 = <$in_fh>;

  my $line_5 = <$in_fh>;
  my $line_6 = <$in_fh>;
  my $line_7 = <$in_fh>;
  my $line_8 = <$in_fh>;

  next if ($line_2 =~ /N/);
  next if ($line_6 =~ /N/);

  $line_1 =~ s/@/>/;
  $line_5 =~ s/@/>/;

  print $out_fh_1 "$line_1$line_2";
  print $out_fh_2 "$line_5$line_6";

}

close $out_fh_1;
close $out_fh_2;
close $in_fh;
