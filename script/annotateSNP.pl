#!/usr/bin/perl -w
use Bio::SeqIO;
my $snp = shift;
my $gff = shift;
my $fna = shift;
my $loc = shift;
my $out = shift;
unless (defined ($snp) and defined ($gff) and defined($fna) and defined $loc and defined ($out)){
	print "Usage: perl $0 <SNP(list)> <GFF> <REF(fasta)> <location file(from StrainLocater)> <OUT>\n";
	print "Copyright: $0 v3.0, wirted by Dalong Hu at 2016/12/22\n~~~~~~~~~~~~~~~~\n";
	exit;
}
sub getPattern{
	my $input = $_[0];
	my @input = split /\t/,$input;
	my $i=1;
	for($j=0;$j<@input;$j++){
		if($input[$j]!~/\d/){
			$input=~s/$input[$j]/$i/g;
			$i++;
			@input = split /\t/,$input;
		}   
	}   
	$input=~s/\t//g;
	$input=~tr/1234/ABCD/;
	return $input;
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
sub codon2trans{
	my $tmp = $_[0];
	my @tmp = split /\-\>/,$tmp;
	my $from = $tmp[0];
	my $to = $tmp[1];
	my @from = split //,$from;
	my @to = split //,$to;
	my $i =0;
	foreach my $base (@from){
		if($base ne $to[$i]){
			if(($base=~/[AG]/ and $to[$i]=~/[AG]/) or ($base =~/[TC]/ and $to[$i]=~/[TC]/)){
				return "transition ";
			}
			else{
				return "transversion ";
			}
		}
		$i++;
	}
}
my %rootSNP;
open LOC,$loc;
my $rootLine = <LOC>;
$rootLine = <LOC>;
my $rootPattern;
my @rootSNP;
if(defined $rootLine){
	chomp $rootLine;
	my @tmpRoot = split /\t/,$rootLine;
	my $side = $tmpRoot[3];
	$rootPattern = $tmpRoot[2]; 
	my $on = $tmpRoot[6];
	my $off = $tmpRoot[7];
	my @on = ();
	if (defined $on){
		@on = split /,/,$on;
	}
	my @off = ();
	if (defined $off){
		@off = split /,/,$off;
	}
	if($side eq 'A'){
		@rootSNP = @on;
	}
	else{
		@rootSNP = @off;
	}
	foreach my $rootSNP(@rootSNP){
		$rootSNP{$rootSNP}++;
	}
}
close LOC;
my %gff;
open GFF,"$gff";
while(<GFF>){
	if(/\tCDS\t/){
		my @tmp = split /\t/,$_;
		$gff{$tmp[3]}{"type"} = 'CDS';
		$gff{$tmp[3]}{"end"} = $tmp[4];
		$gff{$tmp[3]}{"strand"} = $tmp[6];
		$tmp[8] =~ /product=([^\;]+)/;
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
		my @tmp = split /\t/,$_;
		$gff{$tmp[3]}{"type"} = 'rRNA';
                $gff{$tmp[3]}{"end"} = $tmp[4];
                $gff{$tmp[3]}{"strand"} = $tmp[6];
                $tmp[8] =~ /product=([^\;]+)/;
                $gff{$tmp[3]}{"product"} = $1;
		$gff{$tmp[3]}{"id"} = $1;
	}
	if(/\ttRNA\t/){
		my @tmp = split /\t/,$_;
                $gff{$tmp[3]}{"type"} = 'tRNA';
                $gff{$tmp[3]}{"end"} = $tmp[4];
                $gff{$tmp[3]}{"strand"} = $tmp[6];
                $tmp[8] =~ /product=([^\;]+)/;
                $gff{$tmp[3]}{"product"} = $1;
		$gff{$tmp[3]}{"id"} = $1;
	}
}
close GFF;

foreach my $gffEntry (sort {$a<=>$b} keys %gff){
	unless(defined $gff{$gffEntry}{"type"}){
		delete $gff{$gffEntry};
	}
}

open SNP,"$snp";
my $firstline = <SNP>;
chomp $firstline;
my %snp;
my %list;
while(<SNP>){
	chomp;
	my @tmp = split;
	my $location = $tmp[0];
	my $class;
	for (my $i=1; $i<@tmp;$i++){
		unless (defined($class)){
			$class = $tmp[$i];
		}
		else{
			if ($class !~ /$tmp[$i]/){
				$class.=$tmp[$i];
			}	
		}
	}
	$snp{$location}{"class"} = $class;
	$list{$location}{"data"} = $_;
	s/\d+\t//;
	my $pattern = getPattern($_);
	if($pattern eq $rootPattern){
		unless(defined $rootSNP{$location}){
			$pattern=~tr/AB/BA/;
		}
	}
	$list{$location}{"pattern"} = $pattern;
}
close SNP;

my $fasta = Bio::SeqIO->new(-file=>"$fna",-format=>'fasta');
my $ctg = $fasta->next_seq;

foreach my $start (sort {$a<=>$b} keys %gff){
	foreach my $location(sort {$a<=>$b} keys %snp){
		if ($location>=$start and $location<=$gff{$start}{"end"}){
			$list{$location}{"start"} = $start;
			$list{$location}{"end"} = $gff{$start}{"end"};
			$list{$location}{"strand"} = $gff{$start}{"strand"};
			$list{$location}{"id"} = $gff{$start}{"id"};
			$list{$location}{"product"} = $gff{$start}{"product"};
			$list{$location}{"type"} = $gff{$start}{"type"};
			if ($gff{$start}{"type"} eq 'CDS'){
				my $r;
				my $l;
				my $codon;
				my $aa;
				if ($gff{$start}{"strand"} eq '-'){
					$r = $gff{$start}{"end"} - int(($gff{$start}{"end"} - $location)/3)*3;
					$l = $r - 2;
					my $rbase = $ctg->subseq($r,$r);
					$rbase=~tr/[a-z]/[A-Z]/;
					$rbase=~tr/ATGC/TACG/;
					my $lbase = $ctg->subseq($l,$l);
					$lbase=~tr/[a-z]/[A-Z]/;
                                        $lbase=~tr/ATGC/TACG/;
					my $mbase = $ctg->subseq($l+1,$l+1);
                                        $mbase=~tr/[a-z]/[A-Z]/;
                                        $mbase=~tr/ATGC/TACG/;
					$codon = $rbase . $mbase . $lbase;
					if($codon =~/[^ATGC]/){
						$aa = 'X';
					}
					else{
						my $tmp_seq = Bio::Seq->new(-seq=>"$codon");
						$aa = $tmp_seq->translate->seq;
					}
					my @tmp = split //,$snp{$location}{"class"};
					for(my $i=1;$i<@tmp;$i++){
						$tmp[$i]=~tr/[a-z]/[A-Z]/;
						$tmp[$i]=~tr/ATGC/TACG/;
						my $tmp_codon;
						my $tmp_aa;
						if($location == $r){
							$tmp_codon = $tmp[$i] . $mbase .$lbase;
							if($tmp_codon=~/[^ATGC]/){
								$tmp_aa = 'X';
							}
							else{
								my $tmp_obj = Bio::Seq->new(-seq=>"$tmp_codon");
								$tmp_aa = $tmp_obj->translate->seq;
							}
						}
						elsif($location == $l){
							$tmp_codon = $rbase . $mbase .$tmp[$i];
							if($tmp_codon=~/[^ATGC]/){
								$tmp_aa = 'X';
							}
							else{
								my $tmp_obj = Bio::Seq->new(-seq=>"$tmp_codon");
								$tmp_aa = $tmp_obj->translate->seq;
							}
						}
						elsif($location < $r and $location >$l){
							$tmp_codon = $rbase . $tmp[$i] . $lbase;
                                                        if($tmp_codon=~/[^ATGC]/){
								$tmp_aa = 'X';
							}
							else{
								my $tmp_obj = Bio::Seq->new(-seq=>"$tmp_codon");
                                                        	$tmp_aa = $tmp_obj->translate->seq;
							}
						}
						else{
							print "error in $location\n";
							exit;
						}
						my $tmp_trans = "$codon\-\>$tmp_codon";
						$list{$location}{"codon"} .= "$tmp_trans ";
						$list{$location}{"trans"} .= codon2trans($tmp_trans);
						$list{$location}{"aa"} .= "$aa\-\>$tmp_aa ";
						if ($aa eq $tmp_aa){
							$list{$location}{"snp"} .= 's ';
						}
						else{
							$list{$location}{"snp"} .= 'ns ';
						}
					}
						
				}
				else{
				$l = $start+int(($location - $start)/3)*3;
					$r = $l + 2;
					my $rbase = $ctg->subseq($r,$r);
                                        $rbase=~tr/[a-z]/[A-Z]/;
                                        my $lbase = $ctg->subseq($l,$l);
                                        $lbase=~tr/[a-z]/[A-Z]/;
                                        my $mbase = $ctg->subseq($l+1,$l+1);
                                        $mbase=~tr/[a-z]/[A-Z]/;
                                        $codon = $lbase . $mbase . $rbase;
                                        if($codon=~/[^ATGC]/){
						$aa = 'X';
					}
					else{
						my $tmp_seq = Bio::Seq->new(-seq=>"$codon");
                                       		$aa = $tmp_seq->translate->seq;
                                        }
					my @tmp = split //,$snp{$location}{"class"};
                                        for(my $i=1;$i<@tmp;$i++){
                                                $tmp[$i]=~tr/[a-z]/[A-Z]/;
                                                my $tmp_codon;
                                                my $tmp_aa;
                                                if($location == $r){
                                                        $tmp_codon = $lbase . $mbase .$tmp[$i];
                                                        if($tmp_codon=~/[^ATGC]/){
								$tmp_aa = 'X';
							}
							else{
								my $tmp_obj = Bio::Seq->new(-seq=>"$tmp_codon");
								$tmp_aa = $tmp_obj->translate->seq;
							}
                                                }
                                                elsif($location == $l){
                                                        $tmp_codon = $tmp[$i] . $mbase .$rbase;
                                                        if($tmp_codon=~/[^ATGC]/){
								$tmp_aa = 'X';
							}
							else{
								my $tmp_obj = Bio::Seq->new(-seq=>"$tmp_codon");
                                                        	$tmp_aa = $tmp_obj->translate->seq;
							}
                                                }
                                                elsif($location < $r and $location >$l){
                                                        $tmp_codon = $lbase . $tmp[$i] . $rbase;
                                                        if($tmp_codon=~/[^ATGC]/){
								$tmp_aa = 'X';
							}
							else{
								my $tmp_obj = Bio::Seq->new(-seq=>"$tmp_codon");
                                                        	$tmp_aa = $tmp_obj->translate->seq;
							}
                                                }
                                                else{
                                                        print "error in $start\t and \t$location\n";
                                                        exit;
                                                }
												my $tmp_trans ="$codon\-\>$tmp_codon";
                                                $list{$location}{"codon"} .= "$tmp_trans ";
												$list{$location}{"trans"} .= codon2trans($tmp_trans);
                                                $list{$location}{"aa"} .= "$aa\-\>$tmp_aa ";
                                                if ($aa eq $tmp_aa){
                                                        $list{$location}{"snp"} .= 's ';
						}
						else{
                                                        $list{$location}{"snp"} .= 'ns ';
                                                }
                                        }
				}
			}
			else{
				my @tmp = split //,$snp{$location}{"class"};
				for(my $i =1;$i<@tmp;$i++){
					my $tmp_trans = "$tmp[0]\-\>$tmp[$i]";
					$list{$location}{"codon"} .= "$tmp_trans ";
					$list{$location}{"trans"} .= codon2trans($tmp_trans);
				}
				$list{$location}{"aa"} = 'N/A' ;
				$list{$location}{"snp"} = 'nc';
			}		
			delete $snp{$location};
		}
	}
}
foreach my $location (sort {$a<=>$b} keys %snp){
	$list{$location}{"snp"}= 'nc';
	$list{$location}{"aa"} = 'N/A';
        $list{$location}{"start"} = 'N/A';
        $list{$location}{"end"} = 'N/A';
        $list{$location}{"strand"} = 'N/A';
        $list{$location}{"id"} = 'N/A';
        $list{$location}{"product"} = 'N/A';
        $list{$location}{"type"} = 'Unknown';
	my @tmp = split //,$snp{$location}{"class"};
    for(my $i =1;$i<@tmp;$i++){
		my $tmp = "$tmp[0]\-\>$tmp[$i]";
		$list{$location}{"codon"} .= "$tmp ";
		$list{$location}{"trans"} .= codon2trans($tmp);
	}
}
open OUT,">$out";
print OUT "$firstline\tRegion\tGene\tStart\tEnd\tStrand\tMutation Type\tCodon\tAmino Acid\tCluster Type\tProduct\tPattern\n";
foreach my $location (sort {$a<=>$b} keys %list){
	$list{$location}{"product"} = chASCII($list{$location}{"product"});
	chomp $list{$location}{"product"};
	chomp $list{$location}{"id"};
	my $output = "$list{$location}{\"data\"}\t$list{$location}{\"type\"}\t$list{$location}{\"id\"}\t$list{$location}{\"start\"}\t$list{$location}{\"end\"}\t$list{$location}{\"strand\"}\t$list{$location}{\"trans\"}\t$list{$location}{\"codon\"}\t$list{$location}{\"aa\"}\t$list{$location}{\"snp\"}\t$list{$location}{\"product\"}\t$list{$location}{\"pattern\"}\n";
	$output =~s/ \t/\t/g;
	print OUT "$output";
}
