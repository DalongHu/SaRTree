#!/usr/bin/perl -w
my $in = shift;
my $out = shift;
unless(defined($in) and defined($out)){
	    print "Usage:perl $0 <in> <out>\nCopyright: v1.0, writed by Dalong Hu, 2013 05 17 ~~\n";
		    exit;
}
open IN ,"$in";
open OUT ,">$out";
my $f = <IN>;
print OUT "$f";
while(<IN>){
	tr/\=/\-/;
	unless(/-/ or /N/){
		print OUT "$_";
	}
}
