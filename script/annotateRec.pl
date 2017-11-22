#!/usr/bin/perl -w
use Bio::SeqIO;
my $rec = shift;
my $gff = shift;
my $ref = shift;
my $loc = shift;
my $out = shift;
unless (defined ($rec) and defined ($gff) and defined $ref and $loc and defined ($out)){
	print "Usage: perl $0 <REC(list)> <GFF> <REF(length)> <location file(from StrainLocater)> <OUT>\n";
	print "Copyright: $0 v1.0, wirted by Dalong Hu at 2017/06/30\n~~~~~~~~~~~~~~~~\n";
	exit;
}

sub chASCII{
	my $tmp = $_[0];
	my $i=1;
	while($i){
		if($tmp=~/\%(..)/){
			my $asc = $1;
			my $char = chr(hex("0x$asc"));
			$tmp=~s/\%$asc/$char/;
		}
		else{
			$i=0;
		}
	}
	return $tmp;
}

open LOC,$loc;
my $rootLine = <LOC>;
$rootLine = <LOC>;
my $rootPattern;
if(defined $rootLine){
        chomp $rootLine;
        my @tmpRoot = split /\t/,$rootLine;
        $rootPattern = $tmpRoot[2];
}
close LOC;
my %gff;
open GFF,"$gff";
while(<GFF>){
	if(/\tCDS\t/){
		my @tmp = split /\t/;
		$gff{$tmp[3]}{"type"} = 'CDS';
		$gff{$tmp[3]}{"end"} = $tmp[4];
		$gff{$tmp[3]}{"strand"} = $tmp[6];
		$tmp[8] =~ /product=([^\;\n]+)/;
		$gff{$tmp[3]}{"product"} = $1;
		$tmp[8] =~ /Name=([^;]+)/;
		$gff{$tmp[3]}{"id"} = $1;
	}
	if(/\tgene\t/){
		my @tmp = split /\t/;
		$tmp[8]=~/Name=([^\;]+)/;
		$gff{$tmp[3]}{"name"} = $1;
		$tmp[8]=~/locus_tag=([^\;]+)/;
		$gff{$tmp[3]}{"tag"} = $1;
	}
	if(/\trRNA\t/){
		my @tmp = split /\t/;
		$gff{$tmp[3]}{"type"} = 'rRNA';
                $gff{$tmp[3]}{"end"} = $tmp[4];
                $gff{$tmp[3]}{"strand"} = $tmp[6];
                $tmp[8] =~ /product=([^\;\n]+)/;
                $gff{$tmp[3]}{"product"} = $1;
		$gff{$tmp[3]}{"id"} = $1;
		$gff{$tmp[3]}{"tag"} = $1;
	}
	if(/\ttRNA\t/){
		my @tmp = split /\t/;
                $gff{$tmp[3]}{"type"} = 'tRNA';
                $gff{$tmp[3]}{"end"} = $tmp[4];
                $gff{$tmp[3]}{"strand"} = $tmp[6];
                $tmp[8] =~ /product=([^\;\n]+)/;
                $gff{$tmp[3]}{"product"} = $1;
		$gff{$tmp[3]}{"id"} = $1;
		$gff{$tmp[3]}{"tag"} = $1;
	}
}
close GFF;

foreach my $gffEntry (sort {$a<=>$b} keys %gff){
	unless(defined $gff{$gffEntry}{"type"}){
		delete $gff{$gffEntry};
	}   
}

open REC,"$rec";
my @rec;
my %order;
my $i = 0;
my @list;
while(<REC>){
	chomp;
	my @tmp = split;
	my $start = $tmp[1];
	my $end = $tmp[2];
	$order{$i} = $start;
	$rec[$i]{"s"}=$start;
	$rec[$i]{"e"}=$end;
	if($tmp[0] eq $rootPattern){
		$tmp[0] =~tr/AB/BA/;
	}
	$list[$i]{"pattern"} = $tmp[0];
	$list[$i]{"number"} = $tmp[3];
	$i++;
}
close REC;
my @recBackup = @rec;
foreach my $j (sort {$order{$a}<=>$order{$b}} keys %order){
	my $tmp;
	my $found = 0;
	my $foundEnd = 0;
	if($rec[$j]{"s"} < $rec[$j]{"e"}){
		foreach my $start(sort {$a<=>$b} keys %gff){
			unless (defined $gff{$start}{"end"}){
				next;
			}
			if ($rec[$j]{"s"}<=$gff{$start}{"end"} and $rec[$j]{"e"}>=$start){
				$found = 1;
				unless(defined $tmp){
					$tmp = $list[$j]{"start"} = $gff{$start}{"tag"};
					$list[$j]{"product"} =  chASCII($gff{$start}{"product"});
				}	
				else{
					$list[$j]{"end"} = "->$gff{$start}{\"tag\"}";
					my $eproduct = chASCII($gff{$start}{"product"});
					$list[$j]{"eproduct"} = ' -> '. $eproduct;
					$foundEnd = 1;
				}
				unless($foundEnd){
					$list[$j]{"end"} = '';

				}
			}
		}
	}
	else{
		my @first;
		my @last;
		my @result;
		my $key =1;
		foreach my $start(sort {$a<=>$b} keys %gff){
			unless (defined $gff{$start}{"end"}){
				next;
			}
			if ($rec[$j]{"s"}<=$gff{$start}{"end"} or $rec[$j]{"e"}>=$start){
				$found = 1;
				if($key){
					push @first,$start;
					@result = (@last,@first);
				}
				else{
					push @last,$start;
					@result = (@last,@first);
				}
			}
			else{
				$key = 0;
			}
		}
		$list[$j]{"start"} = $gff{$result[0]}{"tag"};
		$list[$j]{"product"} =  chASCII($gff{$result[0]}{"product"});
		if(@result > 1){
			$list[$j]{"end"} = "->$gff{$result[-1]}{\"tag\"}";
			my $eproduct = chASCII($gff{$result[-1]}{"product"});
			$list[$j]{"eproduct"} = ' -> '. $eproduct;
		}
		else{
			$list[$j]{"end"} = '';
		}
	}
	if($found){
		$rec[$j] = undef;
	}
}
foreach my $l (sort {$order{$a}<=>$order{$b}} keys %order){
	if (defined $rec[$l]){
        	$list[$l]{"start"} = 'N/A';
        	$list[$l]{"end"} = '';
		$list[$l]{"product"} = 'N/A';
	}
}
open OUT,">$out";
print OUT "Start\tEnd\tLength\tSNP\tDensity\tLocation\tAnnotation\tPattern\n";
foreach my $l (sort {$order{$a}<=>$order{$b}} keys %order){
	my $start = $recBackup[$l]{"s"};
	my $end = $recBackup[$l]{"e"};
	my $length = $end - $start;
	if($length <0){
		$length = $ref+$length+1;
	}
	else{
		$length++;
	}
	my $snp = $list[$l]{"number"};
	my $density = sprintf "%.6f",$snp/$length;
	my $pattern = $list[$l]{"pattern"};
	my $from = $list[$l]{"start"};
	my $to = $list[$l]{"end"};
	my $product = $list[$l]{"product"};
	chomp $from;
	chomp $to;
	if ($to ne ''){
		my $eproduct = $list[$l]{"eproduct"};
		print OUT "$start\t$end\t$length\t$snp\t$density\t$from$to\t$product$eproduct\t$pattern\n";
	}
	else{
		print OUT "$start\t$end\t$length\t$snp\t$density\t$from\t$product\t$pattern\n";
	}
}
