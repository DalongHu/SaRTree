#!/usr/bin/perl -w
use Bio::TreeIO;
my $in = shift;
my $out = shift;
my $list = shift;
my $rec = shift;
my $loc = shift;
my $name = shift;
unless (defined($in) and defined($out) and defined($list) and defined($rec) and defined $loc and defined($name)){
	print "Usage: perl $0 <in(nwk)> <out(nex)> <snpsum> <recsum> <location file(from StrainLocater)> <namelist>\n2017-06-30\n2017-06-30\n";
	exit;
}

open LOC,$loc;
my $rootLine = <LOC>;
$rootLine = <LOC>;
my $rootPattern;
if(defined $rootLine){
        chomp $rootLine;
        my @tmpRoot = split /\t/,$rootLine;
        $rootPattern = $tmpRoot[2];
}
close LOC;
$rootPattern=~tr/AB/BA/;

my %hash;
my %rec;
my %name;
my %value;
open LIST,$list;
open REC,$rec;
open NAME,$name;
open TYPE,">$out.type.list";
my $i=0;
while (<NAME>){
	chomp;
	my @tmp = split /\t/,$_;
	$name{$tmp[0]}=$i;
	$i++;
}
while(<LIST>){
	chomp;
	my @tmp = split /\t/,$_;
	$hash{$tmp[0]} = $tmp[1];
}
while(<REC>){
	chomp;
	my @tmp = split /\t/,$_;
	$rec{$tmp[0]} = $tmp[1];

}
my $outT = Bio::TreeIO->new(-file=>">$out",-format=>'nexus');
my $inT = Bio::TreeIO->new(-file=>"$in",-format=>'newick');
my $tree = $inT->next_tree;
my @nodes = $tree->get_nodes(-order=>'d');
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
	print TYPE "$type\n";
	if ($node->is_Leaf){
		my $branch_length= $node->branch_length;
		unless (defined($rec{$type})){
			$rec{$type} = 0;
		}
		unless (defined($hash{$type})){
			$hash{$type} = 0;
		}

		my $value = "$branch_length".'[&events=____'."$hash{$type}".'_____'."$rec{$type}".'___'.',branchID=branch_'.$branchID.']';
		$node->branch_length($value);
#		$value{$id} = $value;
	}
	else{
		unless (defined($rec{$type})){
			$rec{$type} = 0;
		}
		unless (defined($hash{$type})){
			$hash{$type} = 0;
		}
		my $value = '&events='."____$hash{$type}_____$rec{$type}___".',branchID=branch_'.$branchID;
		$node->bootstrap($value);
		my $tmp_branchID = 'branch_'.$branchID;
		$node->id($tmp_branchID);
	}
	$branchID++;
}
# foreach my $leaf(@leaf){
# 	my $id=$leaf->id;
#	my $new = $value{$id};
#	$leaf->id($new);
# }
$outT->write_tree($tree);

