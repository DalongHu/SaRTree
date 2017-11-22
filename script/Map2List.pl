#!/usr/bin/perl -w
my $in = shift;
my $ref = shift;
my $out = shift;
my $outname =shift;
my $namefile= shift;
unless(defined($in) and defined($ref) and defined($out)and defined($outname) and defined($namefile)){
	print "Usage:perl $0 <IN(map)> <REF length> <OUT(list)> <OUT_NAME> <namelist file>\nCopyright: v1.0, writed by Dalong Hu, 2013 05 17 ~~\n";
	exit;
}
open OUT,">$out";
open NAME,">>$namefile";
my @list;
$in=~/([^\/]+)\.map/;
my $name = $1;
print OUT "$outname";
print NAME "$outname\t$name\n";
open MAP,"$in";
while(<MAP>){
	chomp;
	my @entry = split /\t/, $_;
	my $location = $entry[0];
	my $base = $entry[1];
	$list[$location] = $base;
}
for(my $i=1;$i<=$ref;$i++){
	if (defined $list[$i]){
		print OUT "\n$list[$i]";
	}
	else{
		print OUT "\n=";
	}
}
