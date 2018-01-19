#!/usr/bin/perl -w
my $in = shift;
my $list = shift;
my $out = shift;
unless(defined $in and defined $list and defined $out){
	print "Usage: perl $0 <IN> <LIST> <OUT>\n";
	exit;
}

open IN,"$in";
open LIST,"$list";
open OUT,">$out";
my %hash;
while(<LIST>){
	chomp;
	my @tmp = split /\s+/;
	$hash{$tmp[0]} = $tmp[1];
}
my $f = <IN>;
foreach my $key(keys %hash){
	$f=~s/$key/$hash{$key}/;
}
print OUT "$f";

while(<IN>){
	print OUT "$_";
}
close IN;