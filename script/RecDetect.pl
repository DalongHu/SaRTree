#!/usr/bin/perl -w
use Getopt::Long;
my $snpList;
my $seqLength;
my $out;
my $fast;
my $help;
my $strict = 0;
my $minCut = 0;
my $maxCut = 2000;
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
	"fast|f!"=>\$fast,
	"help|h!"=>\$help
);
my $usage = "Usage: perl $0 -s|-snp <snpList> -l|-length <seqLength> -o|-out <out(tag)> [options]\nOptions:\n-s|-snp <string> mandatory, snp list file\n-l|-length <integer> mandatory, reference genome length\n-o|-out <string> mandatory, output tag\n-t|-strict optional, run strict model to filter more recombination regions for strains with high recombination rate(default off)\n-f|-fast optional, run fast model for studies with large number of strains and SNPs, but may be less accurate(default off)\n-n|-min [integer] optional, the minimum cut threshold(default 0)\n-x|-max [integer] optional, the maximum cut threshold(default 20000)\n-p|-step [integer] optional, the adding step of cut threshold(default 10)\n-a|-alpha [float] optional, the significant threshold for ks test(default 0.05)\n-h|-help optional, show this usage\nCopyright: v6.0 Writen by Dalong Hu 23 Jan 2018 TT\n";
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

sub average {
	my $sum_tmp = 0;
	my $num = @_;
	foreach my $tmp (@_){
		$sum_tmp += $tmp;
	}
	$sum_tmp /= $num;
}
sub m_multiply{
	my $x1;
	my $x2;
	my $mm;
	my @x3;
	($x1,$x2,$mm) = @_; 
	for(my $i=0;$i<$mm;$i++){
		for(my $j=0;$j<$mm;$j++){
			my $sum = 0;
			for(my $g=0;$g<$mm;$g++){
				$sum+=${$x1}[$i*$mm+$g]*${$x2}[$g*$mm+$j];
			}
			$x3[$i*$mm+$j] = $sum;
		}
	}
	return @x3;
}
sub m_power{
	my $x1;
	my $ex1;
	my @x2;
	my $ex2;
	my $mm;
	my $nn;
	($x1,$ex1,$mm,$nn) = @_;
	@x2=@{$x1};
	$ex2=$ex1;
	my @order;
	while($nn > 1){
		push @order,$nn;
		$nn = int($nn/2);
	}
	@order = reverse @order;
	foreach my $order (@order){
		my @x3;
		my $ex3;
		@x3=m_multiply(\@x2,\@x2,$mm);
		$ex3 =2*$ex2;
		if(($order%2) ==0){
			@x2 = @x3;
			$ex2 = $ex3;
		}
		else{
			@x2 = m_multiply($x1,\@x3,$mm);
			$ex2 = $ex1+$ex3;
		}
		if($x2[int($mm/2)*$mm+int($mm/2)]>1e140){
			for(my $i=0;$i<$mm*$mm;$i++){
				$x2[$i] = $x2[$i]*1e-140;
			}
			$ex2 +=140;
		}
	}
	return (\@x2,$ex2);
}
sub Knd{
	my $n = $_[0];
	my $d = $_[1];
	my $s = $d*$d*$n;
	if($s>7.24 || ($s>3.76 && $n >99)){
		return 1-2*exp(-(2.00071+0.331/sqrt($n)+1.409/$n)*$s);
	}
	my $kk = int($n*$d)+1;
	my $m = 2*$kk-1;
	my $h = $kk -$n*$d;
	my @hh;
	for(my $i=0;$i<$m;$i++){
		for(my $j=0;$j<$m;$j++){
			if($i - $j +1 <0){
				$hh[$i*$m+$j] = 0;
			}
			else{
				$hh[$i*$m+$j] = 1;
			}
		}
	}
	for(my $i = 0;$i<$m;$i++){
		$hh[$i*$m] -= $h**($i+1);
		$hh[($m-1)*$m+$i] -=$h**($m-$i);
	}
	$hh[($m-1)*$m] += ((2*$h-1>0)?(2*$h-1)**$m:0);
	for(my $i=0;$i<$m;$i++){
		for(my $j=0;$j<$m;$j++){
			if($i-$j+1>0){
				for(my $g = 1;$g<=$i-$j+1;$g++){
					$hh[$i*$m+$j] /= $g;
				}
			}
		}
	}
	my $eH =0;
	my $qq;
	my $eq;
	($qq,$eq) = m_power(\@hh,$eH,$m,$n);
	$s = ${$qq}[($kk-1)*$m+$kk-1];
	for(my $i=1;$i<=$n;$i++){
		$s = $s*$i/$n;
		if($s<1e-140){
			$s*=1e140;
			$eq -= 140;
		}
	}
	$s *=10**$eq;
	return $s;
}


sub ksTest{
	my $num= @_;
	my $mean = average(@_);
	my @CDF;
	foreach my $dis (sort {$a<=>$b} @_) {
		push @CDF, (1-exp(-$dis/$mean));
	}
	my $D = 0;
	my $f0 = 0;
	for (my $i = 0; $i < $num; $i++) {
		my $f1 = ($i+1) / $num;
		my $f2 = $CDF[$i];
		my $delta1 = abs($f0 - $f2);
		my $delta2 = abs($f1 - $f2);
		my $Dt = $delta1 > $delta2 ? $delta1:$delta2;
		if ($Dt > $D) { 
			$D = $Dt;
		}
		$f0 = $f1;
	}
	return 1-Knd($num,$D);
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
print + "RecDetect: Settings...\n";
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
	for(my $cut=$minCut;$cut<=$maxCut;$cut+=$stepCut){
		my $sum = 0;
		my @sample;
		if($fast){
			foreach my $dis(@distence){
				if($dis<$cut){next;}
				push @sample,$dis;
				$sum++;
			}
		}
		else{
			my @loadingSNP;
			if($list{$key}{"distances"}[-1] >= $cut and $list{$key}{"distances"}[0] >= $cut ){
				push @loadingSNP,$list{$key}{"locations"}[0];
			}
			for(my $i=1;$i<@{$list{$key}{"locations"}};$i++){
				if($list{$key}{"distances"}[$i-1] >= $cut and $list{$key}{"distances"}[$i] >= $cut){
					push @loadingSNP,$list{$key}{"locations"}[$i];
					if(@loadingSNP > 1 ){
						my $dis_new = $loadingSNP[-1] - $loadingSNP[-2];
						push @sample,$dis_new;
						$sum++;
					}
				}
			}
			if(@loadingSNP > 1){
				push @sample,$loadingSNP[0]+$seqLength-$loadingSNP[-1];
				$sum++;
			}
		}
		if($sum==0){
			last;
		}
		elsif($count == $sum){
			next;
		}
		else{
			$count = $sum;
		}
		@sample = sort{$a<=>$b} @sample;
		my @data;
		my %rand;
		foreach my $sample(@sample){
			my $jitter;
			do{
				srand();
				$jitter = (rand()-0.5)/500;
				$rand{$jitter}++;
			}while(!(defined $rand{$jitter}));
			my $data = $sample + $jitter;
			push @data,$data;
		}
		my $p = ksTest(@data);
		$pValue{$cut} = $p;
	}
	my $maxPValue = -1;
	my $maxCutTmp = 0;
	my $minFitCut = $maxCut+1;
	foreach my $cut (keys %pValue){
			if($pValue{$cut} > $alpha and $cut < $minFitCut){
				$minFitCut = $cut;
			}
			if($pValue{$cut}>$maxPValue){
				$maxPValue = $pValue{$cut};
				$maxCutTmp = $cut;
			}
	}
	if($strict){
		if($maxPValue > $alpha){
			$list{$key}{"cut"} = $maxCutTmp;
			print CUT "$key\t$maxCutTmp\n";
		}
		else{
			$list{$key}{"cut"} = $maxCut; # here, the empty cut is possible
			print CUT "$key\t$maxCut\n"; # this situation show that maxCut is not big enough or something special existing
		}
	}
	else{
		if($maxPValue > $alpha){
			$list{$key}{"cut"} = $minFitCut;
			print CUT "$key\t$minFitCut\n";
		}
		else{
			$list{$key}{"cut"} = $maxCutTmp;
			print CUT "$key\t$maxCutTmp\n";
		}
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
