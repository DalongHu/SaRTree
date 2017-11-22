#!/usr/bin/perl -w
use Bio::TreeIO;
use Getopt::Long;
my $in;
my $top;
my $help;
my $rular;
my $auto;
my $boot;
my $edge;
GetOptions (
	"in|i=s"=>\$in,
	"top|t!"=>\$top,
	"rul|r:f"=>\$rular,
	"auto|a:s"=>\$auto,
	"boot|b!"=>\$boot,
	"edge|e:i"=>\$edge,
	"help|h!"=>\$help
);
my $usage="-------------------treeview v1.2-------------------\n
A commond line tool to view a tree in text interface.\nUsage:\n
perl $0 -i <tree file(nwk)> [Options]\n
==============Options==============
-top/-t\ttopology only, default close\n
-rul/-r [number]\tlength represented by one '-', default 0.5\n
-auto/-a max|min\tauto scale, use \"max\" or \"min\" as tag for diff scale types, availible when -r not indicated, default close\n
-boot/-b\t show bootstrap value, when indicated, the rular will be reset autimatically to fit the scale, default close\n
-edge/-e [integer]\tshould be lareger than 10, show branch length to short the showing '-'s which are longer than indicated value, default 50\n
-help/-h\tshow this message\n
Copyright: Writed by Dalong 12/Dec/2016 ~TT~\n
";
if($help){
	print $usage;
	exit;
}
unless (defined $in){
	print "-------------------treeview v1.2-------------------\nA commond line tool to view a tree in text interface.\nUsage: perl $0 -i <in(nwk)> [-t] [-b] [-r number] [-a max|min] [-e integer]\nCopyright: Writed by Dalong 12/Dec/2016 ~TT~\n";
	exit;
}
unless (defined $auto){
	$auto = '';
}
unless (defined $rular){
	$rular = 0.5;
}
elsif ($rular <=0){
	$rular = 0.001;
	$auto = '';
}
else{
	$auto = '';
}
unless (defined $edge){
	$edge = 50;
}
if ($edge < 10){
	$edge = 10;
}
my $nwk = Bio::TreeIO->new(-format=>'newick',-file=>"$in");
my $tree = $nwk->next_tree;
my $root = $tree->get_root_node;
my @unroot = $root->each_Descendent;
if (@unroot == 3){
	print "unrooted tree, auto root randomly\n";
	$tree->reroot_at_midpoint($unroot[2],"root");
	$root = $tree->get_root_node;
}
unless($tree->is_binary){
	print "Not a binary tree!\n";
}

my @nodes = $tree->get_nodes;
my $minL=100000;
my $maxL=0;
foreach my $no (@nodes){
	my $tmp = $no->branch_length;
	unless (defined $tmp){
		$tmp = 0;
	}
	if($tmp >0){
		if($tmp < $minL){
			$minL = $tmp;
		}
		if($tmp > $maxL){
			$maxL = $tmp;
		}
	}
}
if ($top){
	$maxL = 1;
	$minL = 1;
}
if($auto ne ''){
	if ($auto eq 'max'){
		$rular = $maxL/20;
	}
	if ($auto eq 'min'){
		$rular = $minL;
	}
}
if($boot){
	if($minL/$rular < 3){
		$rular = $minL/3;
	}
}
my @fig;
my %xy;
my $tagr = $root->internal_id;
$xy{$tagr}[0] = 0; #order
$xy{$tagr}[1] = 0; #depth
$xy{$tagr}[2] = '@'; #id
$xy{$tagr}[3] = 0; #up (-1) or down (1)
$xy{$tagr}[4] = 1; #showed branch length
$xy{$tagr}[5] = undef; #output branch length
$xy{$tagr}[6] = 0; #bootsrtap
preorder($root);
my @lines;
my $min=100000;
foreach my $key(keys %xy){
	my $y = $xy{$key}[0];
	if($y<$min){
		$min = $y;
	}
}
foreach my $key(keys %xy){
	my $y = $xy{$key}[0]-$min;
	my $x = $xy{$key}[1];
	my $name = $xy{$key}[2];
	my $len = $xy{$key}[4];
	$lines[$y] = $xy{$key}[3];
	for(my $i=0;$i<$x;$i++){
		push @{$fig[$y]},' ';
	}
	push @{$fig[$y]},'+';
	if (defined $xy{$key}[5]){
		for(my $i =0; $i<int(($len-length($xy{$key}[5])-2)/2); $i++){
			push @{$fig[$y]},'-';
		}
		push @{$fig[$y]},'/';
		my @tmp = split //,$xy{$key}[5];
		foreach my $tmp(@tmp){
			push @{$fig[$y]},$tmp;
		}
		push @{$fig[$y]},'/';
		for(my $i=int(($len-length($xy{$key}[5])-2)/2)+length($xy{$key}[5])+2;$i<$len;$i++){
			push @{$fig[$y]},'-';
		}
	}
	else{
		for(my $i =0; $i<$len; $i++){
			push @{$fig[$y]},'-';
		}
	}
	push @{$fig[$y]},$name;
	push @{$fig[$y]},"\n";
}
if($boot){
	foreach my $key(keys %xy){
		if (defined $xy{$key}[6] and $xy{$key}[6] >0){
			my $y = $xy{$key}[0]-$min;
			my $mark = $xy{$key}[1]+1;
			my @tmp = split //,$xy{$key}[6];
			for(my $i=0;$i<@tmp;$i++){
				$fig[$y-1][$mark+$i] = $tmp[$i];
			}
		}
	}
	
}
for(my $y=0;$y<@fig;$y++){
	for(my $x=0;$x<@{$fig[$y]};$x++){
		if($fig[$y][$x] eq '+'){
			if($lines[$y] > 0){
				my $z = $y-1;
				while($fig[$z][$x] ne '@'){
					$fig[$z][$x] = '|';
					$z--;
				}
			}
			if($lines[$y] < 0){
				my $z = $y+1;
				while($fig[$z][$x] ne '@'){
					$fig[$z][$x] = '|';
					$z++;
				}
			}
		}
	}
}


foreach my $y(@fig){
	my $line = join '',@{$y};
	print "$line";
}
my $ru;
if ($rular < 0.001 or $rular >1000){
	$ru = sprintf("%.2e",$rular);
}
else{
	$ru = sprintf("%.3f",$rular);
}
print "\n|-|\n$ru\n";
sub preorder{
	my $node = $_[0];
	if($node->is_Leaf){
		return;
	}
	my $anc = $node->internal_id;
	my @two = $node->each_Descendent;
	my $up = $two[0];
	my $tag = $up->internal_id;
	my $length = $up->branch_length;
	if($top){
		$length = 1;
	}
	foreach my $tmp (keys %xy){
		if($xy{$tmp}[0] < $xy{$anc}[0]){
			$xy{$tmp}[0]--;
		}
	}
	my $tmp = int(($length)/$rular);
	if ((($length/$rular) - $tmp) > 0.5){
		$tmp++;
	}
	if ($tmp > $edge){
		$tmp = $edge;
		my $le;
		if ($length>1000 or $length <0.001){
			$le = sprintf("%.2e",$length);
		}
		else{
			$le = sprintf("%.3f",$length);
		}
		$xy{$tag}[5] = $le;
	}
	if ($boot){
		unless($up->is_Leaf){
			$xy{$tag}[6] = $up->id;
		}
	}
	$xy{$tag}[4] = $tmp;
	$xy{$tag}[0] = $xy{$anc}[0]-1;
	$xy{$tag}[1] = $xy{$anc}[1]+$xy{$anc}[4]+1;
	if ($up->is_Leaf){
		$xy{$tag}[2] = $up->id;
	}
	else{
		$xy{$tag}[2] = '@';
	}
	$xy{$tag}[3] = -1;
	
	my $down = $two[1];
	$tag = $down->internal_id;
	$length = $down->branch_length;
	if($top){
		$length = 1;
	}
	foreach my $tmp (keys %xy){
		if($xy{$tmp}[0] > $xy{$anc}[0]){
			$xy{$tmp}[0]++;
		}
	}
	$tmp = int(($length)/$rular);
	if ((($length/$rular) - $tmp) > 0.5){
		$tmp++;
	}
	if ($tmp > $edge){
		$tmp = $edge;
		my $le;
		if ($length>1000 or $length <0.001){
			$le = sprintf("%.2e",$length);
		}
		else{
			 $le = sprintf("%.3f",$length);
		}
		$xy{$tag}[5] = $le;
	}
	if ($boot){
		unless($down->is_Leaf){
			$xy{$tag}[6] = $down->id;
		}
	}
	$xy{$tag}[4] = $tmp;
	$xy{$tag}[0] = $xy{$anc}[0]+1;
	$xy{$tag}[1] = $xy{$anc}[1]+$xy{$anc}[4]+1;
	if ($down->is_Leaf){
		$xy{$tag}[2] = $down->id;
	}
	else{
		$xy{$tag}[2] = '@';
	}
	$xy{$tag}[3] = 1;

	preorder($up);
	preorder($down);
	return;
}
