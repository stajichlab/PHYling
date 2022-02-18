#!/usr/bin/env perl
use strict;
use warnings;
use Statistics::Descriptive;
my $stat = Statistics::Descriptive::Full->new();

while(<>) {
    my $dat = $_;
    chomp($dat);
    $stat->add_data($dat);
}
print join("\t",0,$stat->standard_deviation*2,$stat->mean),"\n";
