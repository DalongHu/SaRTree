#!/usr/bin/perl -w
use Bio::TreeIO;
my $in = shift;
my $list = shift;
my $out = shift;
unless(defined $in and defined $list and defined $out){
	print "Usage: perl $0 <IN> <LIST> <OUT>\n";
	exit;
}
my $file = Bio::TreeIO->new(-file=>"$in",-format=>'newick');
my $inT = $file->next_tree;
my @leaves = $inT->get_leaf_nodes;
open LIST,"$list";
my $outT = Bio::TreeIO->new(-file=>">$out",-format=>'newick');
my %hash;
while(<LIST>){
	chomp;
	my @tmp = split /\s+/;
	$hash{$tmp[1]} = $tmp[0];
}
close LIST;
foreach my $leaf(@leaves){
	my $tmp = $leaf->id;
	$tmp = $hash{$tmp};
	$leaf->id($tmp);
}
$outT->write_tree($inT);
