#!/usr/bin/perl -w
my $in = shift;
my $out = shift;
open IN,"$in";
open OUT,">$out";
my %hash;
my $f = <IN>;
while (<IN>){
	chomp;
	my @tmp = split /\t/,$_;
	$hash{$tmp[-1]}++;
}
foreach my $key(keys %hash){
	print OUT "$key\t$hash{$key}\n";
}
