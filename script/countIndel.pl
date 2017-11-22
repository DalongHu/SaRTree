#!/usr/bin/perl -w
my $in = shift;
my $out = shift;
unless(defined $in and defined $out){
	print "Usage: perl $0 <IN(indel.seq.xls)> <OUT(indel.xls)>\n";
	exit;
}
open IN,"$in";
open OUT,">$out";
my $f = <IN>;
print OUT "$f";
my %hash;
while(<IN>){
	chomp;
	my @tmp = split;
	for(my $i = 2;$i<@tmp;$i++){
		if ($tmp[$i] ne '-'){
			$tmp[$i] = length($tmp[$i]);
		}
	}
	my $tmp = join "\t",@tmp;
	$hash{$tmp[0]} = $tmp;
}
foreach my $key(sort {$a <=> $b} keys %hash){
	print OUT "$hash{$key}\n";
}
