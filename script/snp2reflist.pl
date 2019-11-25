#!/usr/bin/perl -w
my $in =shift;
my $out = shift;
unless (defined $in and defined $out){
	print "Usage:perl $0 <map(folder)> <out>\n";
	exit;
}
my @maps = glob"$in/*.snplist";
my %ref;
open OUT,">$out";
foreach my $map (@maps){
	open IN,"$map";
	while(<IN>){
		my @t = split /\t/;
		$ref{$t[0]}= $t[1];
	}
}

my @seqs = sort {$a <=> $b} keys %ref;
print OUT "Location\/Strain\tref\n";
my $refseq = '';
foreach my $loc (@seqs){
	print OUT "$loc\t$ref{$loc}\n";
}

