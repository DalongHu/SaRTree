#!/usr/bin/perl -w
use Bio::TreeIO;
my $in = shift;
my $out =shift;
unless (defined $in and defined $out){
	print "Usage:perl $0 <In(nwk)> <Out(phyloxml)>\n";
	exit;
	
}
my $file = Bio::TreeIO->new(-file=>$in,-format=>'newick');
my $outT = Bio::TreeIO->new(-file=>">$out",-format=>'phyloxml');
my $tree = $file->next_tree;
my @nodes = $tree->get_nodes(-order=>'d');
my $b = 0;
shift @nodes;
foreach my $node(@nodes){
	unless($node->is_Leaf){
		my $id = 'branch_'.$b;
		$node->id($id);
	}
	$b++;
}
$outT->write_tree($tree);

