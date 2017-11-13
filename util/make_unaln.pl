#!env perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use Bio::SeqIO;
my $dir = 'search';
my $dbdir = 'pep';
my $ext = "fa";
my $outdir = 'aln';
my $idxfile = 'allseq';
my $MIN_COUNT = 4;
my $debug = 0;
GetOptions(
    'v|debug!'  => \$debug,
    'd|dir:s'   => \$dir,
    'e|ext:s'   => \$ext,
    'I|idx:s'   => \$idxfile,
    'db:s'      => \$dbdir,
    'o|out:s'   => \$outdir,
    'm|min:i'   => \$MIN_COUNT, # minimum number of sp in a gene cluster to use it, otherwise skip it
    );
$idxfile = File::Spec->catfile($dbdir,$idxfile);
if( ! -f $idxfile ) {
 `cat $dbdir/*.$ext | perl -p -e 'if( /^>/ ) { s/>(\\S+).+/>\$1/ } else { s/\\*//g }' > $idxfile`;
 if ( ! -f $idxfile ) {
   die "could not create $idxfile for indexing all seqs";
 }
 `cdbfasta $idxfile`;
}
opendir(BEST,$dir) || die "cannot open dir: $dir, $!";

my %by_gene;
for my $file ( readdir(BEST) ) {
    next unless $file =~ /(\S+)\.best/;
    my $stem = $1;
    open(my $fh => "$dir/$file") || die "$dir/$file: $!";
    while(<$fh>) {
	my ($gene,$name) = split;
	$by_gene{$gene}->{$stem} = $name;
    }
}
for my $gene ( keys %by_gene ) {
 my $ct = scalar keys %{$by_gene{$gene}};
    if ( $ct < $MIN_COUNT ) {
     warn("skipping $gene has only $ct genes\n");
     next;
    }
    open(CDBYANK,"| cdbyank $idxfile.cidx -o $outdir/$gene.$ext") || die $!;    
    my $expected = 0;
    while( my ($sp,$seqname) = each %{$by_gene{$gene}} ) {
	print CDBYANK $seqname,"\n";
        $expected++;
    }
    close(CDBYANK);
    my $nct = `grep -c '^>' $outdir/$gene.$ext`;
    if( $nct != $expected ) { 
     warn("num seqs is $nct while expected $expected for $gene.$ext\n");
     my $seqcheck = Bio::SeqIO->new(-format => 'fasta', -file => "$outdir/$gene.$ext");
     my %seen;
     while (my $s = $seqcheck->next_seq ) {
      $seen{$s->display_id}++; 
      warn("stored ",$s->display_id,"\n");
     }
     for my $sp ( keys %{$by_gene{$gene}} ) {
       my $name = $by_gene{$gene}->{$sp};
         if( ! exists $seen{$name} ) {	
  	  warn("cannot find '$name' in $gene file\n");
	 }
    }
#	exit;
    }
}
