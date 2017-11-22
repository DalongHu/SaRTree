#!/usr/bin/perl -w
my $in =shift;
my $out =shift;
unless(defined($in) and defined($out)){
	print "Usage: perl $0 <IN(list)> <OUT(fasta)>\n";
	exit;
}
open IN,"$in";
open OUT,">$out";
my $firstline = <IN>;
close IN;
chomp $firstline;
my @title = split /\t/,$firstline;
for (my $i = 2;$i<=@title;$i++){
	my $tmp = `cut -f $i $in`;
	$tmp=~s/^.*\n//;
	$tmp=~s/\n//g;
	print OUT ">$title[$i-1]\n$tmp\n";
}
