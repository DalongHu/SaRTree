#!/usr/bin/perl -w
use Bio::AlignIO;
my $in = shift;
my $out = shift;
my $fasta = Bio::AlignIO->new(-file=> "$in", -format=>'fasta');
my $phylip = Bio::AlignIO->new(-file=> ">$out",-format=>'phylip');
while(my $seq = $fasta->next_aln){
	$phylip->write_aln($seq);
}
