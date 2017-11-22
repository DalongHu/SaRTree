#!/usr/bin/perl -w
use Bio::SeqIO;
my $indel = shift;
my $gff = shift;
my $ref = shift;
my $loc = shift;
my $out = shift;
unless (defined ($indel) and defined ($gff) and defined $ref and defined $loc and defined ($out)){
	print "Usage: perl $0 <INDEL(list)> <GFF> <REF(length)> <location file(from StrainLocater)> <OUT>\n";
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

open INDEL,"$indel";
my $f = <INDEL>;
chomp $f;
my %indel;
my %list;
while(<INDEL>){
	chomp;
	my @tmp = split;
	my $start = $tmp[0];
	my $end = $tmp[1];
	$indel{$start} = $end;
	$list{$start}{"data"} = $_;
	my $pattern = join "\t", @tmp[2..(@tmp-1)];
	$pattern=~s/\d+/B/g;
	$pattern=~s/-/A/g;
	$pattern=~s/\t//g;
	if($pattern=~/^B/){
		$pattern=~tr/AB/BA/;
	}
	if($pattern eq $rootPattern){
		$pattern=~tr/AB/BA/;
	}
	$list{$start}{"pattern"} = $pattern;
}
close INDEL;

foreach my $location (sort {$a<=>$b} keys %indel){
	my $tmp;
	my $found = 0;
	my $foundEnd = 0;
	foreach my $start(sort {$a<=>$b} keys %gff){
		unless (defined $gff{$start}{"end"}){
			next;
		}
		if ($indel{$location} ne 'N/A'){
			if ($location<=$gff{$start}{"end"} and $indel{$location}>=$start){
				$found = 1;
				unless(defined $tmp){
					$tmp = $list{$location}{"start"} = $gff{$start}{"tag"};
					$list{$location}{"product"} =  chASCII($gff{$start}{"product"});
				}    
				else{
					$list{$location}{"end"} = "->$gff{$start}{\"tag\"}";
					my $eproduct = chASCII($gff{$start}{"product"});
					$list{$location}{"eproduct"} = ' -> '.$eproduct;
					$foundEnd = 1;
				}   
				unless($foundEnd){
					$list{$location}{"end"} = ''; 
				}   
			}   
		}   
		else{
			if ($location<=$gff{$start}{"end"} and $location>=$start){
				$found = 1;
				unless(defined $tmp){
					$tmp = $list{$location}{"start"} = $gff{$start}{"tag"};
					$list{$location}{"product"} =  chASCII($gff{$start}{"product"});
				}   
				else{
					$list{$location}{"end"} = "->$gff{$start}{\"tag\"}";
					 my $eproduct = chASCII($gff{$start}{"product"});
					 $list{$location}{"eproduct"} = ' -> '.$eproduct;
					$foundEnd = 1;
				}   
				unless($foundEnd){
					$list{$location}{"end"} = ''; 
				}   
			}   
		} 
	}
	if($found){
		delete $indel{$location};
	}
}
foreach my $location (sort {$a<=>$b} keys %indel){
        $list{$location}{"start"} = 'N/A';
        $list{$location}{"end"} = '';
	$list{$location}{"product"} = 'N/A';
}
open OUT,">$out";
print OUT "$f\tLocation\tAnnotation\tPattern\n";
foreach my $location (sort {$a<=>$b} keys %list){
	my $data = $list{$location}{"data"};
	my $pattern = $list{$location}{"pattern"};
	my $from = $list{$location}{"start"};
	my $to = $list{$location}{"end"};
	my $product = $list{$location}{"product"};
	my $eproduct = $list{$location}{"eproduct"};
	chomp $from;
	chomp $to;
	if ($to ne ''){
		print OUT "$data\t$from$to\t$product$eproduct\t$pattern\n";
	}
	else{
		print OUT "$data\t$from\t$product\t$pattern\n";
	}
}
