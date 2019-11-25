#!/usr/bin/perl -w
sub recov{
	my $tmp = $_[0];
	$tmp=~tr/atgc/TACG/;
	$tmp=~tr/ATGC/TACG/;
	$tmp = reverse $tmp;
	return $tmp;
}
my $mauve= shift;
my $map = shift;
unless(defined($mauve) and defined($map)){
	print "Usage:perl $0 <mauve> <snplist>\nCopyright: v1.2, writed by Dalong Hu, 2019 11 11\n output snp list and gap list ~~\n";
	exit;
}
open MAUVE,"$mauve";
open OUT,">$map";
open GAP,">$map.gap";
print GAP "gap";
local $/ = undef;
my $in = <MAUVE>;
my @result = split /\=/,$in;
foreach my $result (@result){
	$result=~s/\#[^\n]+\n//g;
	if (($result=~/\> 1/) and ($result=~/\> 2/)){
		$result=~/\> 1\:(\d+)\-(\d+) ([\-\+]) ([^\n]+)\n([^\>]+)/;
		my $start = $1;
		my $end = $2;
		my $strand = $3;
		my $refname = $4;
		my $ref = $5;
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
					if($base ne $seq[$i]){
						print OUT "$j\t$base\t$seq[$i]\n";
					}
				}
				else{
					print GAP "\t$j";
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
