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
	$hash{$tmp[1]} = $tmp[0];
}
while(<IN>){
	my $tmp = $_;
	foreach my $key(keys %hash){
		$tmp=~s/^$key\t/$hash{$key}\t/;
	} 
	print OUT "$tmp";
}
