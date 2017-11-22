#!/usr/bin/perl -w
use Bio::TreeIO;
my $raw = shift;
my $out = shift;
my $snp = shift;
my $rec = shift;
my $loc = shift;
my $name = shift;
unless (defined $raw and defined $out and defined $snp and defined $rec and defined $loc and defined $name){
	print "Usage: perl $0 <rawTree(nwk)> <out(branch sum)> <snpList> <recList> <location file(from StrainLocater)> <namelist>\n2017-06-30\n";
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

my %name;
my %data;
my %branch;
open SNP,$snp;
open REC,$rec;
open NAME,$name;
open TYPE,">$out.pattern2branch.txt";
open SUM,">$out";
my $i = 0;
while (<NAME>){
	chomp;
	my @tmp = split;
	$name{$tmp[1]}=$i;
	$i++;
}
close NAME;
my $f = <SNP>;
while(<SNP>){
	chomp;
	my @tmp = split;
	$data{$tmp[-1]}{"snp"}++;
}
close SNP;
$f = <REC>;
while(<REC>){
	chomp;
	my @tmp = split;
	$data{$tmp[-1]}{"rec"}++;
}
close REC;
my $rawT = Bio::TreeIO->new(-file=>"$raw",-format=>'newick');
my $rawTree = $rawT->next_tree;
my @rawNodes = $rawTree->get_nodes;
my @rawLeaf = $rawTree->get_leaf_nodes;
my $branchID = 0;
foreach my $node(@rawNodes){
	my @type;
	my $descendants;
	if ($node->is_Leaf){
		my $id = $node->id;
		$type[$name{$id}]='A';
		$descendants = $id;
	}
	else{
		foreach my $leaf ($node->get_all_Descendents){
			if ($leaf->is_Leaf){
				my $id = $leaf->id;
				$type[$name{$id}]='A';
				$descendants .="$id,";
			}
		}
		$descendants=~s/,$//;
	}
	foreach my $leaf (@rawLeaf){
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
	$data{$type}{"descendants"}=$descendants;
	$branch{$type} = $branchID;
	print TYPE "$type\t$branchID\n";
	unless (defined $data{$type}{"snp"}){
		$data{$type}{"snp"}  = 0;
	}
	unless (defined $data{$type}{"rec"}){
		$data{$type}{"rec"} = 0;
	}
	if ($node->is_Leaf){
		my $id = $node->id;
		$data{$type}{"strain"} = $id;
		$data{$type}{"type"} = 'outer';
		$data{$type}{"boot"} = '-';
	}
	else{
		$data{$type}{"strain"} = '-';
		$data{$type}{"type"} = 'inner';
		my $boot = $node->id;
		if(defined $boot){
			$data{$type}{"boot"} = $boot;
		}
		else{
			$data{$type}{"boot"} = '-';
		}

	}
	$branchID++;
}
print SUM "Branch\tPattern\tStrain\tType\tBootstrap(%)\tSNPs\tRecombinations\tDescentants\n";

foreach my $type(sort {$branch{$a} <=> $branch{$b}} keys %branch){
	my $oBranch = $branch{$type};
	my $oPattern = $type;
	my $oStrain = $data{$type}{"strain"};
	my $oType = $data{$type}{"type"};
	my $oBoot = $data{$type}{"boot"};
	my $oSnp = $data{$type}{"snp"};
	my $oRec = $data{$type}{"rec"};
	my $oDescendants = $data{$type}{"descendants"};
	print SUM "$oBranch\t$oPattern\t$oStrain\t$oType\t$oBoot\t$oSnp\t$oRec\t$oDescendants\n";
}
