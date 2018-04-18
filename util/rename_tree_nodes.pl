#!env perl
use strict;
use Bio::TreeIO;

# Arguments
# Tree file (with prefix names as the node names)
# Prefix table (from 1KFG/genomes/scripts/make_prefixes.pl)
#  which is tab delimited with 2 columns: 1 is prefix 2 is long name

my $in = Bio::TreeIO->new(-format => 'newick', -file => shift);

my $prefixes = shift;

open(my $fh => $prefixes) || die "$prefixes: $!";
my %map;
while(<$fh>) {
    next if /^Pref/;
    my ($pref,$name) = split;
    $map{$pref} = $name;
}

while( my $tree = $in->next_tree ) {
    for my $node ( grep { $_->is_Leaf } $tree->get_nodes ) {
	my $id = $node->id;
	if( my $lookup = $map{$id} ) {
	    $node->id($lookup);
	} else {
	    warn("no $id in prefix table\n");
	}
    }
    Bio::TreeIO->new(-format => 'newick')->write_tree($tree);
}
