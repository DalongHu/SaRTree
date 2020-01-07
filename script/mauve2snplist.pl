#!/usr/bin/perl -w
sub recov{
	my $tmp = $_[0];
	$tmp=~tr/atgc/TACG/;
	$tmp=~tr/ATGC/TACG/;
	$tmp = reverse $tmp;
	return $tmp;
}
my $mauve= shift;
my $snippy = shift;
my $length = shift;
my $map = shift;
unless(defined $mauve and defined $snippy and defined $length and defined $map){
	print "Usage:perl $0 <mauve> <snippy(tab)> <length(ref)> <snplist>\nCopyright: v1.4, writed by Dalong Hu, 2020 Jan 06\n output snp list and gap list ~~\n";
	exit;
}
my %snippy;
open SNP,$snippy;
my $f = <SNP>;
while(<SNP>){
	my @t = split /\t/;
	if($t[2] eq 'snp'){
		$snippy{$t[1]} = $t[4];
	}
}
close SNP;

open MAUVE,"$mauve";
open OUT,">$map";
open GAP,">$map.gap";
local $/ = undef;
my $in = <MAUVE>;
my @result = split /\=/,$in;
my @gap;
my %region;
foreach my $result (@result){
	$result=~s/\#[^\n]+\n//g;
	if (($result=~/\> 1/) and ($result=~/\> 2/)){
		$result=~/\> 1\:(\d+)\-(\d+) ([\-\+]) ([^\n]+)\n([^\>]+)/;
		my $start = $1;
		my $end = $2;
		my $strand = $3;
		my $refname = $4;
		my $ref = $5;
		$region{$start} = $end;
		$refname=~s/.*\///;
		$ref=~s/\n//g;
		$ref=~tr/atgc/ATGC/;
		$result=~/\> 2\:(\d+)\-(\d+) [\-\+] ([^\n]+)\n([^\>]+)/;
		my $seqst = $1;
		my $seqen = $2;
		my $seqname = $3;
		my $seq = $4;
		$seqname=~s/.*\///;
		$seq=~s/\n//g;
		$seq=~tr/atgc/ATGC/;
		if ($strand eq '-'){
			$ref = recov $ref;
			$seq = recov $seq;
		}
		my @ref = split //,$ref;
		my @seq = split //,$seq;
		my $i = 0;
		my $j = $start;
		foreach my $base (@ref){
			if ($base ne '-'){
				if($seq[$i] ne '-'){
					if($base ne $seq[$i] and defined $snippy{$j}){
						if($snippy{$j} eq $seq[$i]){ #redundant to make sure this is checked. 
							print OUT "$j\t$base\t$seq[$i]\n";
						}
						else{
							print OUT "$j\t$base\t$snippy{$j}\n";;
						}
					}
				}
				else{
					push @gap,$j;
				}
				$j++;
			}
			$i++;
		}
		if ($j-1 != $end){
			print "error index!\n";
		}
	}
}

my $tmps = 1;
foreach my $start (sort {$a <=> $b} keys %region){
	if($tmps < $start){
		for(my $i = $tmps;$i<$start;$i++){
			push @gap,$i;
		}
		$tmps = $region{$start}+1;
	}
	else{
		$tmps = $region{$start}+1;
	}
}
if($tmps <=$length){
	for(my $i = $tmps;$i<$length;$i++){
		push @gap,$i;
	}
}
my $gaploc = join "\t",sort {$a <=> $b} @gap;
print GAP "$gaploc";
