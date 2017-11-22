#!/usr/bin/perl -w
use Statistics::R;
use Getopt::Long;
my $snpList;
my $seqLength;
my $out;
my $help;
my $strict = 0;
my $minCut = 0;
my $maxCut = 20000;
my $stepCut = 10;
my $alpha = 0.05;
GetOptions (
	"snp|s=s"=>\$snpList,
	"length|l=i"=>\$seqLength,
	"out|o=s"=>\$out,
	"strict|t!"=>\$strict,
	"min|n:i"=>\$minCut,
	"max|x:i"=>\$maxCut,
	"step|p:i"=>\$stepCut,
	"alpha|a:f"=>\$alpha,
	"help|h!"=>\$help
);
my $usage = "Usage: perl $0 -s|-snp <snpList> -l|-length <seqLength> -o|-out <out(tag)> [options]\nOptions:\n-s|-snp <string> mandatory, snp list file\n-l|-length <integer> mandatory, reference genome length\n-o|-out <string> mandatory, output tag\n-t|-strict optional, run strict module(default off)\n-n|-min [integer] optional, the minimum cut threshold(default 0)\n-x|-max [integer] optional, the maximum cut threshold(default 20000)\n-p|-step [integer] optional, the adding step of cut threshold(default 10)\n-a|-alpha [float] optional, the significant threshold for ks test(default 0.05)\n-h|-help optional, show this usage\nCopyright: v4.0 Writen by Dalong Hu at 22 Sep 2016 TT\n";
if($help){
	print $usage;
	exit;
}
unless(defined $snpList and defined $seqLength and defined $out){
	print $usage;
	exit;
}

open LIST,"$snpList";
open INFO,">$out.info";
open CUT,">$out.cut";
open REC,">$out.rec";
open LOC,">$out.rec.location";
open REAL,">$out.real";

sub getPattern{
	my $input=$_[0];
	my @input=split/\t/,$input;
	my $i=1;
	for(my $j=0;$j<@input;$j++){
		if($input[$j]!~/\d/){
			$input=~s/$input[$j]/$i/g;
			$i++;
			@input=split/\t/,$input;
		}
	}
	$input=~s/\t//g;
	$input=~tr/1234/ABCD/;
	return $input;
}

my $firstline = <LIST>;
my %list;
my %entry;
print + "RecDetect: Reading ...\n";
while(<LIST>){
	chomp;
	my @entry=split;
	my $location=$entry[0];
	s/\d+\t//;
	my $pattern=getPattern($_);
	$entry{$location} = $_;
	push @{$list{$pattern}{"locations"}},$location;
}
print + "RecDetect: Setting up R ...\n";
my $R = Statistics::R->new();
my %recLocation;
my $totRun = keys %list;
my $haveRun = 0;
my $havePre = 0;
foreach my $key(sort {$a cmp $b} keys %list){
	if (100*$haveRun/$totRun >= $havePre){
		print + "RecDetect: Running $havePre% ...\n";
		$havePre+=10;
	}
	$haveRun++;
	print INFO "$key\t";
	print INFO join "\t",@{$list{$key}{"locations"}};
	print INFO "\n";
	if(@{$list{$key}{"locations"}}<4){
		$list{$key}{"cut"} = $minCut;
		print CUT "$key\ttoo less\n";
		next;
	}
	for(my $i=1;$i<@{$list{$key}{"locations"}};$i++){
		push @{$list{$key}{"distances"}},$list{$key}{"locations"}[$i]-$list{$key}{"locations"}[$i-1];
	}
	push @{$list{$key}{"distances"}},$list{$key}{"locations"}[0]+$seqLength-$list{$key}{"locations"}[-1];

	my @distence=sort{$a<=>$b} @{$list{$key}{"distances"}};
	my %pValue;
	my $count = 0;
	for(my $cut=$minCut;$cut<$maxCut;$cut+=$stepCut){
		my $sum = 0;
		my @sample;
		foreach my $dis(@distence){
			if($dis<$cut){next;}
			push @sample,$dis;
			$sum++;
		}
		if($sum==0){
			next;
		}
		elsif($count == $sum){
			next;
		}
		else{
			$count = $sum;
		}

		$R->set('x',\@sample);
		$R->run(q`t <- rep(0,length(x))`,q`y <- x+jitter(t)`,q`z <- ks.test(y,'pexp',1/mean(y))`);
		my $p = $R->get('z$p.value');
		$pValue{$cut} = $p;
#		print "$key\t$cut\t$p\n";
	}
	my $maxPValue = -1;
	my $maxCutTmp = 0;
	if($strict){
		foreach my $cut (keys %pValue){
			if($pValue{$cut}>$maxPValue){
				$maxPValue = $pValue{$cut};
				$maxCutTmp = $cut;
			}
		}
	}
	else{
		foreach my $cut (sort {$a<=>$b}keys %pValue){
			if($pValue{$cut} > $alpha){
				$maxPValue = $pValue{$cut};
				$maxCutTmp = $cut;
				last;
			}
		}
	}
	if($maxPValue > $alpha){
		$list{$key}{"cut"} = $maxCutTmp;
		print CUT "$key\t$maxCutTmp\n";
	}
	else{
		$list{$key}{"cut"} = $maxCut;
		print CUT "$key\t$maxCut\n"; # this situation show that maxCut is not big enough or something special existing
	}
	
	my @tmpDis = ();
	my @tmpIndex = ();
	my @recEvent;
	my $recEventIndex = 0;
	for(my $i=0;$i<@{$list{$key}{"distances"}};$i++){   
		if($list{$key}{"distances"}[$i]<$list{$key}{"cut"}){   
			push @tmpDis,$list{$key}{"distances"}[$i];
			push @tmpIndex,$i;
		}   
		else{   
			if(@tmpIndex>0){   
				push @tmpIndex,$i;
				foreach my $recLocation(@tmpIndex){
					$recLocation{$list{$key}{"locations"}[$recLocation]}++;
				}
				my $start = $list{$key}{"locations"}[$tmpIndex[0]];
				my $end = $list{$key}{"locations"}[$tmpIndex[-1]];
				$recEvent[$recEventIndex]{"id"} = $tmpIndex[0];
				$recEvent[$recEventIndex]{"endlocation"} = $end;
				my $avg = 0;
				$avg+=$_ for @tmpDis;
				$avg/=@tmpDis;
				$start -= int $avg/2;
				$end += int $avg/2;
				if($start <=0){
					$start+=$seqLength;
				}
				if($end <= 0){
					$end+=$seqLength;
				}
				if($start > $seqLength){
					$start-=$seqLength;
				}
				if($end > $seqLength){
					$end-=$seqLength;
				}
				$recEvent[$recEventIndex]{"start"}=$start;
				$recEvent[$recEventIndex]{"end"}=$end;
				$recEvent[$recEventIndex]{"num"}=@tmpIndex;
				$recEvent[$recEventIndex]{"avg"}=$avg;
				
				#print REC "$key\t$start\t$end\t$num\n";
				$recEventIndex++;
				@tmpDis = ();
				@tmpIndex = ();
			}   
		}
	}
	if(@tmpIndex>0){
		foreach my $recLocation(@tmpIndex){
			$recLocation{$list{$key}{"locations"}[$recLocation]}++;
		}
		$recLocation{$list{$key}{"locations"}[0]}++;
		my $avg = 0;
		$avg+=$_ for @tmpDis;
		my $start;
		my $end;
		if(defined $recEvent[0] and $recEvent[0]{"id"} == 0){
			$avg = ($recEvent[0]{"avg"}*($recEvent[0]{"num"}-1)+$avg)/($recEvent[0]{"num"}-1+@tmpDis);
			$start = $list{$key}{"locations"}[$tmpIndex[0]] - int($avg/2);
			$end = $recEvent[0]{"endlocation"} + int($avg/2);
			if($start <=0){
				$start+=$seqLength;
			}
			if($end <= 0){
				$end+=$seqLength;
			}
			if($start > $seqLength){
				$start-=$seqLength;
			}
			if($end > $seqLength){
				$end-=$seqLength;
			}
			$recEvent[0]{"start"} = $start;
			$recEvent[0]{"end"} = $end;
			$recEvent[0]{"num"} = $recEvent[0]{"num"}+@tmpDis;
		#	print "link\t$key\t$start\t$end\t$recEvent[0]{\"num\"}\n";
		}
		else{
			$start = $list{$key}{"locations"}[$tmpIndex[0]];
			$end = $list{$key}{"locations"}[0];
			$avg/=@tmpDis;
			$start -= int $avg/2;
			$end += int $avg/2;
			if($start <=0){
				$start+=$seqLength;
			}
			if($end <= 0){
				$end+=$seqLength;
			}
			if($start > $seqLength){
				$start-=$seqLength;
			}
			if($end > $seqLength){
				$end-=$seqLength;
			}
			$recEvent[$recEventIndex]{"start"}=$start;
			$recEvent[$recEventIndex]{"end"}=$end;
			$recEvent[$recEventIndex]{"num"}=@tmpDis+1;
		#	print "last\t$key\t$start\t$end\t$recEvent[$recEventIndex]{\"num\"}\n";
		}
	}
	foreach my $recEvent(@recEvent){
		my $tmpStart = ${$recEvent}{"start"};
		my $tmpEnd = ${$recEvent}{"end"};
		my $tmpNum = ${$recEvent}{"num"};
		print REC "$key\t$tmpStart\t$tmpEnd\t$tmpNum\n";	
	}
}
$R->stop();
print + "RecDetect: R stopping ...\n";
print + "RecDetect: Writing ...\n";
foreach my $recLocation(sort {$a<=>$b} keys %recLocation){
	delete $entry{$recLocation};
	print LOC "$recLocation\n";
}
print REAL "$firstline";
foreach my $key (sort {$a<=>$b} keys %entry){
	print REAL "$key\t$entry{$key}\n";
}
print + "RecDetect: finished!\n";
