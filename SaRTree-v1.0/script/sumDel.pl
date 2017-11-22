#!/usr/bin/perl -w
#use Bio:SeqIO;
my $del = shift;
my $incompleteSNP =shift;
my $sum = shift;
unless (defined $del and defined $incompleteSNP and defined $sum){
	print "Usage: perl $0 <deletion list> <out(incomplete SNP)> <out(sum of del)>\n";
	exit;
}
open DEL,$del;
open SNP,">$incompleteSNP";
open SUM,">$sum";
my $f = <DEL>;
print SNP "$f";
$f=~s/^[^\t]+\t//;
#print SUM "Start\tEnd\t$f";
my $start = 0;
my $end = 0;
my @seq;
my $tmppattern;
while(<DEL>){
	chomp;
	my @tmp = split /\t/;
	my @base = @tmp[1..(@tmp-1)];
	my @n = (0..(@base-1));
	foreach my $tmp (@tmp[2..(@tmp-1)]){
		if ($tmp ne $tmp[1] and $tmp ne '-'){
			print SNP "$_\n";
			last;
		}
	}
	my $pattern = join '',@tmp[1..(@tmp-1)];
	$pattern =~s/\w/A/g;
	$pattern =~s/-/B/g;
	unless($start){
		$start = $tmp[0];
		$end = $tmp[0];
		@seq = @base;
		$tmppattern = $pattern;
		next;
	}
	if ($tmp[0]-$end == 1 and $tmppattern eq $pattern){
		$end = $tmp[0];
		$seq[$_].=$base[$_] for @n;
	}
	else{
		my $seq = join "\t",@seq;
		$seq=~s/-+/-/g;
		print SUM "$start\t$end\t$seq\n";
		$start = $tmp[0];
		$end = $tmp[0];
		$tmppattern = $pattern;
		@seq = @base;
	}
}
my $seq = join "\t",@seq;
$seq=~s/-+/-/g;
print SUM "$start\t$end\t$seq\n";
