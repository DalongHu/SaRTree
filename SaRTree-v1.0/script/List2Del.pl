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
	my @entry = split /\t/,$_;
	if ($_=~/-/ or $_=~/=/){
		s/=/-/g;
		print OUT "$_";
	}
}
