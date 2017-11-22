#!/usr/bin/perl -w
my $in = shift;
my $list = shift;
my $name =shift;
my $out = shift;
unless(defined $in and defined $list and defined $name and defined $out){
	print "Usage: $0 <in> <list> <name> <out>\n";
	exit;
}
open LST,"$list";
my %hash;
while(<LST>){
	chomp;
	my @tmp = split;
	$hash{$tmp[0]} = $tmp[1];
}
open NAME,"$name";
my %name;
while(<NAME>){
	chomp;
	my @tmp = split;
	$name{$tmp[0]}= $tmp[1];
}
open IN,"$in";
open OUT,">$out";
my $f = <IN>;
chomp $f;
foreach my $id(keys %name){
	$f =~ s/$id/$name{$id}/g;
}
print OUT "$f\tBranch\n";
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	my $branch;
	if(defined $hash{$tmp[-1]}){
		$branch =$hash{$tmp[-1]};
	}
	else{
		$branch = '-';
	}
	print OUT "$_\t$branch\n";
}

