#!/usr/bin/perl -w
use Bio::TreeIO;
use Getopt::Long;
my $list;#normal list
my $snp; # the snp.real file or annotated snp list
# my $rec; # next to include rec information
my $tree; # next to include support of nex format
my $out; # out tag
my $ext; # extention of the result to more detailed branches, but maybe wrong
my $low; # low accuracy for "out of tree" detecting
GetOptions (
	"list|l=s"=>\$list,
	"snp|s=s"=>\$snp,
	"tree|t=s"=>\$tree,
	"out|o=s"=>\$out,
	"noext|ne!"=>\$noext,
	"low|w!"=>\$low
);
unless (defined $list and defined $snp and defined $tree and defined $out){
	print "Usage: perl $0 -l <list(sartree output)> -s <snp(list for target tree)> -t <tree(nwk)> -o <out(tag)> [-ne] [-w]\n-list|-l\tinput snp list of queries\n-snp|-s\tsnp list of target tree\n-tree|-t\ttarget tree file(nwk format)\n-out|-o\ta prefix for output files(could include path)\n-noext|-ne\tno extention algorithm, not recommonded\n-low|-l\tlow accuracy algorithm, but could locate more strains with low quality\nWrited by Dalong Hu, new version fitting for normal input and for tree rooting output 2017-08-09\n";
	exit;
}
unless(defined $noext){
	$noext=0;
}

unless(defined $low){
	$low=0;
}
sub getPattern{
	my $input = $_[0];
	my @input = split /\t/,$input;
	my $i=1;
	my @base;
	for($j=0;$j<@input;$j++){
		if($input[$j]!~/\d/){
			push @base,$input[$j];
			$input=~s/$input[$j]/$i/g;
			$i++;
			@input = split /\t/,$input;
		}
	}   
	$input=~s/\t//g;
	$input=~tr/1234/ABCD/;
	return ($input,\@base);
}
my $inT = Bio::TreeIO->new(-file=>"$tree",-format=>'newick');
my $t = $inT->next_tree;
my @nodes = $t->get_nodes;
my @leaf = $t->get_leaf_nodes;

open SNP,$snp;
my %snp;
my %pattern2snp;
my $f = <SNP>;
chomp $f;
my @names =  split /\t/,$f;
@names = @names[1..@leaf];
my %names;
my $i = 0;
foreach my $name(@names){
	$names{$name}=$i;
	$i++;
}
while(<SNP>){
	chomp;
	my @tmp = split;
	my $location = $tmp[0];
	@tmp = @tmp[1..@leaf];
	my $tmp = join "\t",@tmp;
	my ($pattern,$class) = getPattern($tmp);
	$snp{$location}{"pattern"} = $pattern;
	$snp{$location}{"ref"} = shift @$class;
	$snp{$location}{"bases"} = $class; # maybe wrong here
	push @{$pattern2snp{$pattern}},$location;
}
close SNP;
my $j = 0;
my %nodes;
my %hash;
my %id;
my %id2pattern;
my %lnames;
my %side;
my %error;
my @printNodes;
foreach my $node (@nodes){
	my @type;
	my @lnames;
	if ($node->is_Leaf){
		my $id = $node->id;
		unless (defined $names{$id}){
			$error{$id}++;
		}
		$type[$names{$id}]='A';
	}
	else{
		foreach my $leaf ($node->get_all_Descendents){
			if ($leaf->is_Leaf){
				my $id = $leaf->id;
				unless (defined $names{$id}){
					$error{$id}++;
				}
				$type[$names{$id}]='A';
			}
		}
	}
	foreach my $leaf (@leaf){
		my $id = $leaf->id;
		unless (defined $names{$id}){
			$error{$id}++;
		}
		unless (defined $type[$names{$id}]){
			$type[$names{$id}]='B';
		}
		else{
			push @lnames,$id;
		}
	}
	my $type = join '',@type;
	my $r = 'A';
	if ($type[0] eq 'B'){
		$type =~tr/AB/BA/;
		$r = 'B';
	}
	if ($type!~/B/){
		$j++;
		next; #root
	}
	$nodes{$type}{$r}{"id"} = $j; # for the two nodes which are direct descendents of the root
	my $intID = $node->internal_id;
	$hash{$intID} = $type; # maybe error
	$id{$intID} = $j;
	$id2pattern{$j} = $type;
	my $lnames = join ',',@lnames; # ',' maybe not good
	$lnames{$j} = $lnames;
	$side{$j} = $r;
	$printNodes[$j-1] = "$type\t$r\t$lnames";
	$j++;
}
foreach my $error (keys %error){
	print "$error\n";
}

open IN,$list;
my $f2 = <IN>;
chomp $f2;
my @f2 = split /\t/,$f2;
my $str_num = @f2;
my %in;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	my $tmp_loc = $tmp[0];
	unless(defined $snp{$tmp_loc}){
		next;
	}
	for (my $w = 2;$w<$str_num;$w++){
		my $tmp_name = $f2[$w];
		$in{$tmp_name}{$tmp_loc} = $tmp[$w];
	}
}
close IN;
my %result;
my %record;
foreach my $name(sort {$a cmp $b} keys %in){
	foreach my $l (keys %snp){
		unless (defined $nodes{$snp{$l}{"pattern"}}){
			next; # possible pattern unsupport the tree
		}
		my $scoreA = 0;
		my $scoreB = 0;
		if (defined $in{$name}{$l}){
			if ($in{$name}{$l} eq $snp{$l}{"ref"}){ # waste codes now
				if (defined $nodes{$snp{$l}{"pattern"}}{'A'}){
					$scoreA++; # for unkonwn information of the snps of two nodes around the root
					push @{$record{$name}{$nodes{$snp{$l}{"pattern"}}{'A'}{"id"}}{"in"}},$l;
					# next to support rooted tree
				}
				else{
					$scoreB--;
					push @{$record{$name}{$nodes{$snp{$l}{"pattern"}}{'B'}{"id"}}{"out"}},$l;
					
				}
			}
			else{
				foreach my $base(@{$snp{$l}{"bases"}}){ #maybe wrong here
					if ($in{$name}{$l} eq $base){
						if (defined $nodes{$snp{$l}{"pattern"}}{'B'}){
							$scoreB++;
							push @{$record{$name}{$nodes{$snp{$l}{"pattern"}}{'B'}{"id"}}{"in"}},$l;
						}   
						else{
							$scoreA--;
							push @{$record{$name}{$nodes{$snp{$l}{"pattern"}}{'A'}{"id"}}{"out"}},$l;
						}
					} # no else, "else" here equal to new base, 0 score, will add this to support better rooting method
				}
			}
		}
		else{ # regard as same to ref now, next to support gap
			if (defined $nodes{$snp{$l}{"pattern"}}{'A'}){
				$scoreA++;
				push @{$record{$name}{$nodes{$snp{$l}{"pattern"}}{'A'}{"id"}}{"in"}},$l;
			}   
			else{
				$scoreB--;
				push @{$record{$name}{$nodes{$snp{$l}{"pattern"}}{'B'}{"id"}}{"out"}},$l;
			}
		}
		if ($scoreA!=0){
			unless (defined $result{$name}{$nodes{$snp{$l}{"pattern"}}{'A'}{"id"}}){
				$result{$name}{$nodes{$snp{$l}{"pattern"}}{'A'}{"id"}} = $scoreA;
			}
			else{
				$result{$name}{$nodes{$snp{$l}{"pattern"}}{'A'}{"id"}} += $scoreA;
			}
		}
		if ($scoreB!=0){
			unless (defined $result{$name}{$nodes{$snp{$l}{"pattern"}}{'B'}{"id"}}){
				$result{$name}{$nodes{$snp{$l}{"pattern"}}{'B'}{"id"}} = $scoreB;
			}   
			else{
				$result{$name}{$nodes{$snp{$l}{"pattern"}}{'B'}{"id"}} += $scoreB;
			}   
		}  
	}
}
my %max;
my %pattern2located;
open OUT,">$out.outOfTree.txt";
open LOG,">$out.onTree.txt";
open LOG2,">$out.onTree.simple.txt";
print LOG "Strain\tbrID\tbrPattern\trefSide\tScore\tDescendents\tLocations on\tLocations out\tLocations missing\n";
print LOG2 "Strain\tbrID\tNumbers on\tNumbers out\tNumbers missing\n";
foreach my $strain (keys %result){
	my $max = -9999999;
	my $maxp = '';
	my $maxk;
	foreach my $k (keys %{$result{$strain}}){
		my $tmp = $nodes[$k];
		$tmpID = $tmp->internal_id;
		if (defined $nodes{$hash{$tmpID}}){ #possible pattern unsupport the tree
			my $score = $result{$strain}{$k}; # algorithm here is bad
			my $anc = $nodes[$k];
			while ($anc = $anc->ancestor){ # maybe error when it is the root
				my $id = $anc->internal_id;
				if (defined $id{$id}){
					if (defined $result{$strain}{$id{$id}}){ # possible no snp supporting a node
						$score+=$result{$strain}{$id{$id}};
					}
				}

			}
			if ($score > $max){ # next to support equal scores
				$max = $score;
				$maxp = $hash{$tmpID};
				$maxk = $k;
			}
		}
	}
	if ($max <= 0) { # could be other value depending on diff accuracies	
		my $tmpNoRecord;
		unless(defined $record{$strain}{$maxk}{"in"}){
			$tmpNoRecord = 0;
		}
		else{
			$tmpNoRecord = @{$record{$strain}{$maxk}{"in"}};
		}
		if($low){
			if($tmpNoRecord == 0 and $max < 0){
				print OUT "$strain\n";
				next;
			}
		}
		else{
			print OUT "$strain\n";
			next;
		}
	}

	unless($noext){ # for unroot tree
		my @descendents =  $nodes[$maxk]->each_Descendent;
		my %tmpScore;
		my $tmpNoRecord;
		unless(defined $record{$strain}{$maxk}{"out"}){
			$tmpNoRecord = 0;	
		}
		else{
			$tmpNoRecord = @{$record{$strain}{$maxk}{"out"}};
		}
		$tmpScore{$maxk}{"ori"}=$tmpNoRecord;
		foreach my $descendent (@descendents){
			my $desID = $descendent->internal_id;
			unless(defined $record{$strain}{$id{$desID}}{"in"}){
				$tmpNoRecord = 0;
			}
			else{
				$tmpNoRecord = @{$record{$strain}{$id{$desID}}{"in"}};
			}
			if($tmpNoRecord == 0){
				next;# at least one "in" snp on descendent branch, or it couldn't be extent
			}
			unless(defined $result{$strain}{$id{$desID}}){
				$result{$strain}{$id{$desID}} = 0;
			}
			$tmpScore{$id{$desID}}{"ori"} = $tmpNoRecord;
		}
		my $tmpSumScore = 0;
		foreach my $tmpBranchID(keys %tmpScore){
			$tmpSumScore += $tmpScore{$tmpBranchID}{"ori"};
		}
		foreach my $tmpBranchID(keys %tmpScore){
			$tmpScore{$tmpBranchID}{"sum"} = $tmpSumScore - $tmpScore{$tmpBranchID}{"ori"};
		}
		my $tmpMink;
		my $tmpMinSumScore = 999999999;
		my $testEqual = 0;
		foreach my $tmpBranchID(keys %tmpScore){
			if($tmpScore{$tmpBranchID}{"sum"}<$tmpMinSumScore){
				$tmpMinSumScore = $tmpScore{$tmpBranchID}{"sum"};
				$tmpMink = $tmpBranchID; 
				$testEqual = 0;
			}
			elsif($tmpScore{$tmpBranchID}{"sum"} == $tmpMinSumScore){
				my $testTotalSNP = $tmpScore{$tmpBranchID}{"ori"}+abs($result{$strain}{$tmpBranchID}) - $tmpScore{$tmpMink}{"ori"}+abs($result{$strain}{$tmpMink});
				if($testTotalSNP <0){ # when paradox snps are equal in two branches, it is more possible that paradox happens in a branch with more total snps 
					$tmpMinSumScore = $tmpScore{$tmpBranchID}{"sum"};
					$tmpMink = $tmpBranchID;
					$testEqual = 0;
				}
				elsif($testTotalSNP == 0){
					$testEqual = 1;
				}
			}
		}

		if($maxk != $tmpMink and $testEqual == 0){ # when the extreme situation (equal branches could be located) happens, will give up extention
			$max += $result{$strain}{$tmpMink};
			$maxk = $tmpMink;
			$maxp = $id2pattern{$maxk};
		}
	}

	unless (defined $max{$maxk}){
		$max{$maxk} = "$strain";
	}
	else{
		$max{$maxk} .= "\|$strain"; #maybe | is not good
	}
	my $lnames = $lnames{$maxk};
	my $side = $side{$maxk};
	my $l_on = '-';
	my $n_on = 0;
	my $l_out = '-';
	my $n_out = 0;
	if(defined $record{$strain}{$maxk}{"in"}){
		$l_on = join ',',@{$record{$strain}{$maxk}{"in"}};
		$n_on = @{$record{$strain}{$maxk}{"in"}};
	}
	else{
		@{$record{$strain}{$maxk}{"in"}} = ();
	}
	if(defined $record{$strain}{$maxk}{"out"}){
		$l_out = join ',',@{$record{$strain}{$maxk}{"out"}};
		$n_out = @{$record{$strain}{$maxk}{"out"}};
	}
	else{
		@{$record{$strain}{$maxk}{"out"}} = ();
	}
	my $branchRealID = $maxk-1;
	my @checkedSNP = (@{$record{$strain}{$maxk}{"in"}},@{$record{$strain}{$maxk}{"out"}});
	my %checkedSNP;
	foreach my $checkedSNP(@checkedSNP){
		$checkedSNP{$checkedSNP}++;
	}
	my $l_missing = '-';
	my @l_missing;
	my $n_missing = 0;
	foreach my $missing(@{$pattern2snp{$maxp}}){
		unless(defined $checkedSNP{$missing}){
			push @l_missing,$missing;
		}
	}
	if(defined $l_missing[0]){
		$l_missing = join ',',@l_missing;
		$n_missing = @l_missing;
	}
	print LOG "$strain\t$branchRealID\t$maxp\t$side\t$max\t$lnames\t$l_on\t$l_out\t$l_missing\n";
	print LOG2 "$strain\t$branchRealID\t$n_on\t$n_out\t$n_missing\n";
}

open NOD,">$out.nodes.txt";
print NOD "ID\tNumberSNP\tNumberLocated\tLocated\tbrPattern\trefSide\tDescendents\n";
open NOD2,">$out.nodes.simple.txt";
print NOD2 "ID\tNumberSNP\tNumberLocated\tLocated\t\n";
my $i_printnode = 0;
foreach my $printnode(@printNodes){
	my $innernodeID = $i_printnode+1;
	my $number_snp;
	if(defined $pattern2snp{$id2pattern{$innernodeID}}){ 
		$number_snp = @{$pattern2snp{$id2pattern{$innernodeID}}};
	}
	else{
		$number_snp = 0;
	}
	my $number_located = 0;
	my $strain_located = '-';
	if (defined $max{$innernodeID}){
		$strain_located = $max{$innernodeID};
		my @nlnames = split /\|/,$strain_located;
		$number_located = @nlnames;
		print NOD2 "$i_printnode\t$number_snp\t$number_located\t$strain_located\n";
	}
	print NOD "$i_printnode\t$number_snp\t$number_located\t$strain_located\t$printnode\n";
	$i_printnode++;
}

close NOD;
close NOD2;
my $outT = Bio::TreeIO->new(-file=>">$out.tree.nex",-format=>'nexus');
my $m = 0;
foreach my $node(@nodes){
	my $branchRealID = $m - 1;
	if (defined $max{$m}){
		my $lnames = $max{$m}; # lnames is a wrong name for that
		my @nlnames = split /\|/,$lnames;
		my $nlnames = @nlnames;
			if ($node->is_Leaf){
				my $branch_length= $node->branch_length;
				my $value = "$branch_length".'[&located='."$lnames".',numberlocated='.$nlnames.',branchID='.$branchRealID.']';
				$node->branch_length($value);
			}
			else{
				my $value = '&located='."$lnames".',numberlocated='.$nlnames.',branchID='.$branchRealID;
				$node->bootstrap($value);
			}
	}
	else{
		if ($node->is_Leaf){
			my $branch_length= $node->branch_length;
			my $value = "$branch_length".'[&branchID='.$branchRealID.']';
			$node->branch_length($value);
		}
		else{
			my $value = '&branchID='.$branchRealID;
			$node->bootstrap($value);
		}
	}
	$m++;
}

$outT->write_tree($t);






