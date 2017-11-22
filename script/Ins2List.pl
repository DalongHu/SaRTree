#!/usr/bin/perl -w
my $indir = shift;
my $ref = shift;
my $namelist = shift;
my $out = shift;
unless(defined $indir and defined $namelist and defined $ref and defined $out){
	print "Usage:perl $0 <IN(dir)> <REF(fasta)> <Name(lsit)> <OUT(list)>\nCopyright: v1.0, writed by Dalong Hu, 2017 04 5 ~~\n";
	exit;
}

my @file = glob "$indir/*.map.ins";
open NAME,"$namelist";
my %n;
my @id;
while(<NAME>){
	chomp;
	my @t = split /\t/;
	$n{$t[1]} = $t[0];
	push @id,$t[0];
}
shift @id;
close NAME;
open OUT,">$out";
print OUT "Start\tEnd\tref";
my %list;
my %loc;
foreach my $file (@file){
	$file=~/([^\/]+)\.map/;
	my $name = $1;
	my $id = $n{$name};
	open MAP,"$file";
	while(<MAP>){
		chomp;
		my @entry = split /\t/;
		my $location = $entry[0];
		my $base = $entry[1];
		$list{$id}{$location} .= $base;
		$loc{$location}++;
	}
}
foreach my $id (@id){
	print OUT "\t$id";
}
print OUT "\n";
open REF,"$ref";
my $firstline = <REF>;
local $/ = undef;
my $refseq = <REF>;
$refseq =~ s/\n//g;
local $/ = "\n";
my @refseq = split //,$refseq;
my $i = 1;
my %h;
foreach my $base (@refseq){
	my $k = 0;
	if (defined $loc{$i}){
		$k = 1;
	}
	if ($k){
		print OUT "$i\tN/A\t-";
		foreach my $tmp (@id){
			if (defined $list{$tmp}{$i}){
				print OUT "\t$list{$tmp}{$i}";
			}
			else{
				print OUT "\t-";
			}
		}
		print OUT "\n";
	}
	$i++;
}
