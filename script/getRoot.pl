#!/usr/bin/perl -w
use Bio::TreeIO;
my $in = shift;
my $snp = shift;
my $out = shift;
unless(defined $in and $snp and defined $out){
	print "Usage: perl $0 <IN(tree)> <snp(list)> <OUT(location)>\nWirted by Dalong Hu, 2017/06/30\n";
	exit;
}

my $file = Bio::TreeIO->new(-file=>"$in",format=>'newick');
my $tree = $file->next_tree;
my $root = $tree->get_root_node;
my @leaf = $tree->get_leaf_nodes;
my @des = $root->each_Descendent;
open SNP,$snp;
my $f = <SNP>;
chomp $f;
my @names =  split /\t/,$f;
@names = @names[1..@leaf];
my $ref= $names[0];
my %names;
my $i = 0;
foreach my $name(@names){
        $names{$name}=$i;
        $i++;
}
close SNP;


open OUT,">$out";
print OUT "Strain\tbrID\tbrPattern\trefSide\tScore\tDescendents\tLocations on\tLocations out\tLocations missing\n";
foreach my $node (@des){
	my @type;
	if ($node->is_Leaf){
                my $id = $node->id;
                $type[$names{$id}]='A';
        }
        else{
                foreach my $leaf ($node->get_all_Descendents){
                        if ($leaf->is_Leaf){
                                my $id = $leaf->id;
				unless(defined $names{$id}){
					print "$id\n";
				}
                                $type[$names{$id}]='A';
                        }
                }
        }
        foreach my $leaf (@leaf){
                my $id = $leaf->id;
		unless (defined $names{$id}){
			print "$id\n";
		}
                unless (defined $type[$names{$id}]){
                        $type[$names{$id}]='B';
                }
        }
        my $type = join '',@type;
        if ($type[0] eq 'A'){
        	print OUT "-\t-\t$type\tB\t-\t-\t-\t-\t-\n";
	}
	
}
