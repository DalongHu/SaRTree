#!/usr/bin/perl -w
use Bio::TreeIO;
my $raw = shift;
my $in = shift;
my $out = shift;
my $snp = shift;
my $rec = shift;
my $indel = shift;
my $loc = shift;
my $name = shift;
my $date = shift;
unless (defined $raw and defined $in and defined $out and defined $snp and defined $rec and defined $indel and defined $loc and defined $name and defined $date){
	print "Usage: perl $0 <rawTree(nwk)> <in(nex)> <out(branch sum)> <snpList> <recList> <indelList> <location file(from StrainLocater)> <namelist> <date>\n2017-06-30\n";
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
my %index;
my %data;
my %branch;
my %date;
open SNP,$snp;
open REC,$rec;
open INDEL,$indel;
open NAME,$name;
open DATE,$date;
open TYPE,">$out.pattern2branch.txt";
open SUM,">$out";
my $i = 0;
while (<NAME>){
	chomp;
	my @tmp = split;
	$index{$tmp[0]}=$tmp[1];
	$name{$tmp[1]}=$i;
	$i++;
}
close NAME;
my $last;
while(<DATE>){
	chomp;
	my @tmp = split;
	$date{$tmp[0]} = $tmp[1];
	unless(defined $last){
		$last = $tmp[1];
	}
	elsif($tmp[1] > $last){
		$last = $tmp[1];
	}
}
close DATE;
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
$f = <INDEL>;
while(<INDEL>){
	chomp;
	my @tmp = split;
	$data{$tmp[-1]}{"indel"}++;
}
close INDEL;
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
				$descendants .= "$id,";
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
	$branch{$type} = $branchID;
	$data{$type}{"descendants"}=$descendants;
	print TYPE "$type\t$branchID\n";
	unless (defined $data{$type}{"snp"}){
		$data{$type}{"snp"}  = 0;
	}
	unless (defined $data{$type}{"rec"}){
		$data{$type}{"rec"} = 0;
	}
	unless (defined $data{$type}{"indel"}){
		$data{$type}{"indel"} = 0;
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
my $inT = Bio::TreeIO->new(-file=>"$in",-format=>'nexus');
my $tree = $inT->next_tree;
my @nodes = $tree->get_nodes;
my @leaf = $tree->get_leaf_nodes;
my $root = $tree->get_root_node;
my $height = $root->height;
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
	if ($type[0] eq 'B'){
		$type=~tr/AB/BA/;
	}
	if ($type!~/B/){
		next;
	}
	my $length = $node->branch_length;
	my $depth = $node->depth;
	my $end = $height-$depth;
	$data{$type}{"start"} = int($last-$end-$length);
	if ($node->is_Leaf){
		my $id = $node->id;
		$data{$type}{"end"} = $date{$id};
	}
	else{
		$data{$type}{"end"} = int($last-$end);
	}
	$data{$type}{"year"} = sprintf "%.2f",$length;
	$data{$type}{"rate"} = sprintf "%.3f",$data{$type}{"snp"}/$length;
}
print SUM "Branch\tPattern\tStrain\tStart Date\tEnd Date\tLength(years)\tEvoluationary Rate(SNPs/year)\tType\tBootstrap(%)\tSNPs\tRecombinations\tIndels\n";

foreach my $type(sort {$branch{$a} <=> $branch{$b}} keys %branch){
	my $oBranch = $branch{$type};
	my $oPattern = $type;
	my $oStrain = $data{$type}{"strain"};
	my $oStart = $data{$type}{"start"};
	my $oEnd = $data{$type}{"end"};
	my $oLength = $data{$type}{"year"};
	my $oRate = $data{$type}{"rate"};
	my $oType = $data{$type}{"type"};
	my $oBoot = $data{$type}{"boot"};
	my $oSnp = $data{$type}{"snp"};
	my $oRec = $data{$type}{"rec"};
	my $oIndel = $data{$type}{"indel"};
	my $oDescendants = $data{$type}{"descendants"};	
	print SUM "$oBranch\t$oPattern\t$oStrain\t$oStart\t$oEnd\t$oLength\t$oRate\t$oType\t$oBoot\t$oSnp\t$oRec\t$oIndel\t$oDescendants\n";
}
