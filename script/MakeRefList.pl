#!/usr/bin/perl -w
my $in = shift;
my $out = shift;
my $namefile = shift;
unless(defined($in) and defined($out) and defined($namefile)){
	print "Usage: perl $0 <IN(ref_fna)> <OUTDIR(dir)> <namelistfile>\n";
	exit;
}
open IN,"$in";
$in=~/([^\/]+)\.[^\.]+/;
my $name = $1;
open OUT ,">$out/ref.list";
open NAME,">>$namefile";
print OUT "Location/Strains\tref";
print NAME "ref\t$name\n";
my $firstline = <IN>;
local $/ = undef;
my $seq = <IN>;
$seq=~s/\n//g;
local $/ = "\n";
my @seq = split //,$seq;
for (my $i = 1;$i<=@seq;$i++){
	print OUT "\n$i\t$seq[$i-1]";
}
