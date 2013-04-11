#!/usr/bin/perl
use strict;

my ($in_1_fn, $in_2_fn, $out_fn) = @ARGV;

open my $in_fh_1, "<", $in_1_fn;
open my $in_fh_2, "<", $in_2_fn;
open my $out_fh, ">", $out_fn;

my $num_processed = 0;
my $num_removed = 0;

while(<$in_fh_1>)
{
    $num_processed++;
  my $line_1 = $_;
  my $line_2 = <$in_fh_1>;
  my $line_3 = <$in_fh_1>;
  my $line_4 = <$in_fh_1>;

  my $line_5 = <$in_fh_2>;
  my $line_6 = <$in_fh_2>;
  my $line_7 = <$in_fh_2>;
  my $line_8 = <$in_fh_2>;

  if (($line_2 =~ /N/) || ($line_6 =~ /N/)) {$num_removed++; next;}

  $line_1 =~ s/@/>/;
  $line_5 =~ s/@/>/;

#  $line_1 =~ s/ length=36/_1/;
#  $line_5 =~ s/ length=36/_2/;

  print $out_fh "$line_1$line_2$line_5$line_6";
}

print "---------------------------------------------\n  Report: $out_fn\n---------------------------------------------\n";
print "  Number Processed:\t$num_processed\n  Number Removed:\t$num_removed\t(".($num_removed/$num_processed*100)."%)\n";
print "---------------------------------------------\n";
close $in_fh_1;
close $in_fh_2;
close $out_fh;
