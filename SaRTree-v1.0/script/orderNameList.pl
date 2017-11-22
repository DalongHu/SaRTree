#!/usr/bin/perl -w
my $in =shift;
my $list = shift;
my $out = shift;
open IN,$in;
my %n;
while(<IN>){
	chomp;
	my @t = split /\t/;
	$n{$t[0]} = $t[1];
}
close IN;
open LIST, $list;
my @order;
my $f = <LIST>;
chomp $f;
@order = split /\t/,$f;
shift @order;
close LIST;

open OUT,">$out";
foreach my $id(@order){
	my $name = $n{$id};
	print OUT "$id\t$name\n";
}
