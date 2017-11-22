#!/usr/bin/perl -w
my $in = shift;
my $out = shift;
unless(defined($in) and defined($out)){
	    print "Usage:perl $0 <in> <out>\nCopyright: v1.0, writed by Dalong Hu, 2013 05 17 ~~\n";
		    exit;
}

open IN ,"$in";
open OUT ,">$out";
my $firstline = <IN>;
print OUT "$firstline";
while(<IN>){
	my $entry = $_;
	chomp;
	my @entry = split /\t/,$_;
	my $loaction = $entry[0];
	my $i=1;
	my $tmp = '+';
	foreach my $base (@entry){
		if ($i){
			$i=0;
			next;
		}
		if ($tmp eq '+'){
			$tmp = $base;
		}
		elsif($base ne $tmp and $base =~ /[ATGC]/){
			print OUT "$entry";
			last;
		}
	}
}
