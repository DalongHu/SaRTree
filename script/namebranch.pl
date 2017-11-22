#!/usr/bin/perl -w
use Bio::TreeIO;
my $in = shift;
my $out = shift;
my $list = shift;
my $loc = shift;
my $name = shift;
unless (defined($in) and defined($out) and defined($list) and defined($rec) and defined $indel and defined $loc and defined($name)){
	print "Usage: perl $0 <in(nwk)> <out(tag)> <branch list> <location file(from StrainLocater)> <namelist>\n";
	exit;
}

open LOC,$loc;
my $rootLine = <LOC>;
$rootLine = <LOC>;
my $rootPattern;
if(defined $rootLine){
        chomp $rootLine;
        my @tmpRoot = split /\t/,$rootLine;
        $rootPattern = $tmpRoot[1];
}
close LOC;
$rootPattern=~tr/AB/BA/;


my %hash;
my %name;
my %value;
open LIST,$list;
open NAME,$name;
my $i=0;
while (<NAME>){
	chomp;
	my @tmp = split /\t/,$_;
	$name{$tmp[0]}=$i;
	$i++;
}
close NAME;
while(<LIST>){
	chomp;
	my @tmp = split /\t/,$_;
	my $branch_no = 'branch_'.$tmp[0];
	$hash{$branch_no}{"p"} = $tmp[1];
	$hash{$branch_no}{"t"} = $tmp[3];
	$hash{$branch_no}{"d"} = $tmp[-1];
}
close LIST;
my $outT = Bio::TreeIO->new(-file=>">$out.nex",-format=>'nexus');
my $outT2 = Bio::TreeIO->new(-file=>">$out.xml",-format=>'phyloxml');
my $inT = Bio::TreeIO->new(-file=>"$in",-format=>'newick');
my $tree = $inT->next_tree;
my @nodes = $tree->get_nodes;
my @leaf = $tree->get_leaf_nodes;
my $l=@leaf;
my $n = @nodes;
foreach my $node(@nodes){
	my @type;
	if ($node->is_Leaf){
		my $id = $node->id;
		$type[$name{$id}]='A';
	}
	else{
		foreach my $leaf ($node->get_all_Descendents){
			if ($leaf->is_Leaf){
				my $id = $leaf->id;
				$type[$name{$id}]='A';
			}
		}
	}
	foreach my $leaf (@leaf){
		my $id = $leaf->id;
		unless (defined($type[$name{$id}])){
			$type[$name{$id}]='B';
		}
	}
	my $type = join '',@type;
	if ($type[0] eq 'B' and $type ne $rootPattern){
		$type=~tr/AB/BA/;
	}
	if ($type!~/B/){
		next;
	}
	if ($node->is_Leaf){
		my $branch_length= $node->branch_length;
		my $value = "$branch_length".'[&events=__'."$hash{$type}".'__'."$rec{$type}".'__'."$indel{$type}".'__]';
		$node->branch_length($value);
	}
	else{
		my $value = '&events='."__$hash{$type}__$rec{$type}__$indel{$type}__";
		$node->bootstrap($value);
	}
}
$outT->write_tree($tree);
$outT2->write_tree($tree);

