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
	print "Usage:perl $0 <mauve> <map>\nCopyright: v1.0, writed by Dalong Hu, 2013 05 17 ~~\n";
	exit;
}
open MAUVE,"$mauve";
open OUT,">$map";
open INS,">$map.ins";
local $/ = undef;
my $in = <MAUVE>;
my @result = split /\=/,$in;
foreach my $result (@result){
	$result=~s/\#[^\n]+\n//g;
	if (($result=~/\> 1/) and ($result=~/\> 2/)){
		$result=~/\> 1\:(\d+)\-(\d+) ([\-\+]) .+\n([^\>]+)/;
		my $start = $1;
		my $end = $2;
		my $strand = $3;
		my $ref = $4;
		$ref=~s/\n//g;
		$result=~/\> 2\:\d+\-\d+ [\-\+] .+\n([^\>]+)/;
		my $seq = $1;
		$seq=~s/\n//g;
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
				print OUT "$j\t$seq[$i]\n";
				$j++;
			}
			else{
				print INS "$j\t$seq[$i]\n";
			}
			$i++;
		}
		if ($j-1 != $end){
			print "error index!\n";
		}
	}
}
