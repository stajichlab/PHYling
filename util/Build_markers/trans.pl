#!/usr/bin/env perl
#
use strict;
use warnings;
use Bio::SeqIO;

my $d = shift || ".";
my $odir = shift || 'pep';
opendir(DIR,$d);
mkdir($odir) if -d $odir;
foreach my $file ( readdir(DIR) ) {
    if ( $file =~ /(\S+)\.txt$/ ) {
	my $base = $1;
	my $in = Bio::SeqIO->new(-format => 'fasta', -file => $file);
	my $out = Bio::SeqIO->new(-format => 'fasta', -file => ">$odir/$base.pep");
	while ( my $seq = $in->next_seq ) {
	    my $keep = undef;
	    for my $frame ( 0..2 ) {
		my $t = $seq->translate(-frame => $frame);
		if ( $t->seq !~ /\*/ ) {
		    $keep = $t;
		    last;
		}
	    }
	    if ( $keep ) {
		$out->write_seq($keep);
	    }
	}
    }
}
