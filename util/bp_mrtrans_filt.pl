#!/usr/bin/env perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
use strict;

# Author:      Jason Stajich <jason-at-bioperl-dot-org>
# Description: Perl implementation of Bill Pearson's mrtrans
#              to project protein alignment back into cDNA coordinates
#

=head1 NAME

bp_mrtrans - implement a transformer of alignments from protein to mrna coordinates

=head1 SYNOPSIS

Usage:
  bp_mrtrans -i inputfile -o outputfile [-if input format] [-of output format] [-s cDNA sequence database]  [-sf cDNA sequence format] [-h]

=head1 DESCRIPTION

This script will convert a protein alignment back into a cDNA.  Loosely
based on Bill Pearson's mrtrans.

The options are:

   -o filename          - the output filename [default STDOUT]
   -of format           - output sequence format
                          (multiple sequence alignment)
                          [default phylip]
   -i filename          - the input alignment [required] (if clipkit filtered, provide clipkit logfile to --filter)
   -if format           - input sequence format
                          (multiple sequence alignment)
                          [ default clustalw]
   -s --seqdb filename  - the cDNA sequence database file
   -sf --seqformat      - the cDNA seq db format (flatfile sequence format)
   --filter             - log file from clipkit processing of the protein alignment
   -h                   - this help menu

=head1 AUTHOR

Jason Stajich, jason-at-bioperl-dot-org

=cut

use strict;
use warnings;
use Bio::AlignIO;
use Bio::SeqIO;
use Getopt::Long;
my %STOP_CODONS = ( 'TAA' => 1, 'TAG' => 1, 'TGA' => 1 );

use constant CODONSIZE => 3;
our $GAP       = '-';
our $CODONGAP  = $GAP x CODONSIZE;

# add our own version of aa_to_dna_aln to support filter
# use Bio::Align::Utilities qw(aa_to_dna_aln);

sub aa_to_dna_aln {
    my ( $aln, $dnaseqs, $filtercols ) = @_;
    unless ( defined $aln
        && ref($aln)
        && $aln->isa('Bio::Align::AlignI') ) {
        croak(
'Must provide a valid Bio::Align::AlignI object as the first argument to aa_to_dna_aln, see the documentation for proper usage and the method signature'
        );
    }
    my $alnlen   = $aln->length;
    my $dnaalign = Bio::SimpleAlign->new();
    $aln->map_chars( '\.', $GAP );
    my @pepseqs;
    foreach my $seq ( $aln->each_seq ) {
        my $aa_seqstr = $seq->seq();
        for ( my $i = 0 ; $i < $alnlen ; $i++ ) {
          my $char = substr( $aa_seqstr, $i, 1 );
          if ( $char eq '*') {
            warn("stop codon found at $i\n");
            $filtercols->[$i] = 1
          }
        }
        push @pepseqs, [$seq->display_id,$aa_seqstr];
    }
    foreach my $s ( @pepseqs) {
        my ($pepid,$pepseq) = @$s;
        my $dnaseq = $dnaseqs->{$pepid} || die "cannot find $pepid\n";
        $dnaseq    = uc $dnaseq->seq();
        my $dnalen = $dnaseqs->{$pepid}->length;
        my $dnaid  = $dnaseqs->{$pepid}->display_id || $pepid; # try to use DNAseq obj ID (issue #137)
        my $nt_seqstr;
        my $j = 0; # this is the CDS sequence counter
        # this is 0-based in the alignment space

        for ( my $i = 0 ; $i < $alnlen ; $i++ ) {
            my $char = substr( $pepseq, $i, 1 );

            if ( $char eq $GAP || $j >= $dnalen ) {
                if ( ! $filtercols || ! $filtercols->[$i] ) {
                  # if filtercols was not given OR this is NOT listed as a filtered column then add the gap
                  $nt_seqstr .= $CODONGAP;
                }
            } else {
              # if filtercols was not given OR this is NOT listed as a filtered column then add the gap
                if( ! $filtercols || ! $filtercols->[$i] ) {
		    my $codon = substr( $dnaseq, $j, CODONSIZE );
		    if ( ! exists $STOP_CODONS{$codon}) {
			$nt_seqstr .= $codon;
		    }
                }
                # always advance the CDS counter since we are skipping this codon
                $j += CODONSIZE;
            }
        }
        $nt_seqstr .= $GAP x ( ( $alnlen * 3 ) - length($nt_seqstr) );

        my $newdna = Bio::LocatableSeq->new(
            -display_id => $dnaid,
            -alphabet   => 'dna',
            -strand     => 1,
            -seq        => $nt_seqstr
        );
        $dnaalign->add_seq($newdna);
    }
    return $dnaalign;
}

# TODO - finish documentation,
#      - add support for extra options in output alignment formats
#        such as idnewline in phylip out to support Molphy input files

my ($iformat,$seqformat,$oformat,$seqdb,$input,$output,$clipkitfilter) = ('clustalw','fasta',
									  'phylip');
my ($help,$usage);

$usage = "usage: bp_mrtrans.pl -i prot_alignment -if align_format -o out_dna_align -of output_format -s cDNA_seqdb -sf fasta --filter clipkit_proteinaln.log\n".
"defaults: -if clustalw
          -of phylip
          -sf fasta\n";

GetOptions(
	   'if|iformat:s'  => \$iformat,
	   'i|input:s'     => \$input,
	   'o|output:s'    => \$output,
	   'of|outformat:s'=> \$oformat,
	   's|seqdb|db:s'  => \$seqdb,
           'filter:s'      => \$clipkitfilter,
	   'sf|seqformat:s'=> \$seqformat,
	   'h|help'        => sub{ exec('perldoc',$0);
				   exit(0)
				   },
	   );

$input ||= shift;
$seqdb ||= shift;
$output ||= shift;
if( ! defined $seqdb ) {
    die("cannot proceed without a valid seqdb\n$usage");
}
if( ! defined $input ) {
    die("cannot proceed without an input file\n$usage");
}

my $indb = new Bio::SeqIO(-file => $seqdb,
			  -format=>$seqformat);
my %seqs;
while( my $seq = $indb->next_seq ) {
    $seqs{$seq->id} = $seq;
}

my @filtersites;
if ( $clipkitfilter ) {
    open(my $fh => $clipkitfilter) || die "cannot open $clipkitfilter";
    # this is 1-based indexing, subtract 1 from colID to get array position
    while(<$fh>) {
	     my ($col,$status,$parsimonyInformative,$gappiness) = split;
       $filtersites[$col-1] = $status eq 'trim' ? 1 : 0;
       #[$status,$parsimonyInformative,$gappiness];
    }
}
my $in = new Bio::AlignIO(-format => $iformat,
			  -file   => $input);
my $out = new Bio::AlignIO(-format => $oformat,
			   -idlength => 22,
			   -interleaved => 0,
			   defined $output ? (-file   => ">$output") : () );

while( my $aln = $in->next_aln ) {
    my $dnaaln = &aa_to_dna_aln($aln,\%seqs,\@filtersites);
    $dnaaln->set_displayname_flat(1);
    $out->write_aln($dnaaln);
}

__END__
