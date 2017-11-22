#!/usr/bin/perl -w
use Bio::TreeIO;
my $in = shift;
my $location = shift;
my $out = shift; #tag
unless(defined $in and defined $location and defined $out){
	print "The script to root a tree.\nUsage: perl $0 <in(nwk)> <location file(from StrainLocater)> <out(tag)>\nCopyright: Wrote by Dalong on 30-Jun-2017\n";
	exit;
}
my $infile = Bio::TreeIO->new(-format=>'newick',-file=>"$in");
my $intree = $infile->next_tree;
my @leaves = $intree->get_leaf_nodes;
my %name;
my %node;
my $i=0;
foreach my $leaf (@leaves){
	my $id= $leaf->id;
	push @name,$id;
	$name{$id}=$i;
	$i++;
}
sub linknodes{
	%node = ();
	my @nodes = $_[0]->get_nodes;
	foreach my $node(@nodes){
		my @tag;
		my @descendents = $node->get_all_Descendents;
		my $number = @descendents;
		unless($number){
			push @descendents,$node;
		}
		foreach my $descendent(@descendents){
			if($descendent->is_Leaf){
				my $id = $descendent->id;
				$tag[$name{$id}] = 'A';
			}
		}
		for(my $j=0;$j<@name;$j++){
			unless(defined $tag[$j]){
				$tag[$j] = 'B';
			}
		}
		my $tag = join '',@tag;
		$node{$tag} = $node;
	}
}
linknodes($intree);
open LOC,$location;
my $line = <LOC>;
$line = <LOC>;
close LOC;
if(defined $line){
	chomp $line;
	my @tmp = split /\t/,$line;
	my $side = $tmp[3];
	my $des = $tmp[5];
	my $on = $tmp[6];
	my $off = $tmp[7];
	my @des = split /,/,$des;
	my @on=();
	if(defined $on) {
		@on = split /,/,$on;
	}
	my @out=();
	if(defined $off) {
		@out = split /,/,$off;
	}
	my @tag;
	foreach $des(@des){
		$tag[$name{$des}] = 'A';
	}
	for(my $j=0;$j<@name;$j++){
		unless(defined $tag[$j]){
			$tag[$j] = 'B';
		}
	}
	my $tag = join '',@tag;
	print "$tag\n";
	foreach my $key(keys %node){
		print "$key\n";
	}
	my $root = $node{$tag};
	$intree->reroot_at_midpoint($root);
	my $tmp = Bio::TreeIO->new(-file=>">$out.tmp.nwk",-format=>'newick');
	$tmp->write_tree($intree);
	$tmp = undef;
	open TMP,"$out.tmp.nwk";
	open TOUT,">$out.tmp.out.nwk";
	my $line = <TMP>;
	close TMP;
	while($line=~/\(([^\(\),]+):([^\(\),]+)\):([^\(\),]+)/){
		my $p1 = $1;
		my $p2 = $2;
		my $p3 = $3;
		my $s = $p2+$p3;
		$line=~s/\($p1:$p2\):$p3/$p1:$s/;
	}
	print TOUT $line;
	close TOUT;
	$infile = Bio::TreeIO->new(-format=>'newick',-file=>"$out.tmp.out.nwk");
	$intree = $infile->next_tree;
	linknodes($intree);
	$root = $node{$tag};
	$tag=~tr/AB/BA/;
	my $new = $node{$tag};
	unless(@on == 0 and @out ==0){
		print @on+1;
		print @out+1;
		my $tmpb = $root->branch_length;
		$tmpb = $tmpb*2*(@on/(@on+@out));
		$root->branch_length($tmpb);
		$tmpb = $new->branch_length;
		$tmpb = $tmpb*2*(@out/(@on+@out));
		$new->branch_length($tmpb);
	}
}
else{
	print "fail\n";
	open LOG,">$out.tmp.out.log";
	print LOG "fail to locate outgroup strain,should use other way to root the tree\n";
}
my $outtree = Bio::TreeIO->new(-file=>">$out.nwk",-format=>'newick');
$outtree->write_tree($intree);
