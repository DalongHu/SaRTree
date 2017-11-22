#!/usr/bin/perl -w
my $dir = shift;
my $out = shift;
unless (defined $dir and defined $out){
	print "Usage: perl $0 <DIR> <OUT>\n";
	exit;
}
my @files = glob "$dir/*.list";
open OUT,">$out";
foreach my $file (@files){
	open IN,$file;
	my $id= <IN>;
	chomp $id;
	$file=~/([^\/]+)\.list/;
	my $name = $1;
	my $i = 0;
	my $q = 0;
	while(<IN>){
		chomp;
		if($_ eq '-' or $_ eq '='){
			$q--;
		}
		$i++;
	}
	close IN;
	my $c = $q*100/$i;
	print OUT "$id\t$name\t$c\n";
}
