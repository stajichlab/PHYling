#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min);

my $adjust = 0.90;

GetOptions(
    'a|adjust:f' => \$adjust,
);    


my @scores = ();

while(<>) {
    chomp;
    push @scores, $_;
}

my $lowest = min(@scores);
printf( "%.2f\n",$lowest * $adjust);
