#!/usr/bin/perl -w
my $snp = shift;
my $rec = shift;
my $mode = shift;
my $out =shift;
unless (defined $snp and defined $rec and defined $mode and defined $out){
	print "Usage:perl $0 <SNP> <REC> <MODE(multi|rec|both|either)> <OUT>\n";
	exit;
}
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

open REC, $rec;
my %rec;
while(<REC>){
	my @t = split /\t/;
	unless(defined $rec{$t[1]}){
		$rec{$t[1]} = $t[2];
	}
	else{
		if($t[2]>$rec{$t[1]}){
			$rec{$t[1]} = $t[2];
		}
	}
}
close REC;

open SNP,$snp;
my $f = <SNP>;
open OUT,">$out";
open OUT2,">$out.filter";
print OUT "$f";
print OUT2 "$f";
while(<SNP>){
	my $ent = $_;
	chomp;
	s/(\d+)\t//;
	my $loc = $1;
	my $pat = getPattern($_);
	if($mode eq 'multi'){
		unless($pat=~/C/){
			print OUT "$ent";
		}
		else{
			print OUT2 "$ent";
		}
	}
	elsif($mode eq 'rec'){
		my $check = 1;
		foreach my $key (sort {$a <=> $b} keys %rec){
			if($loc>=$key and $loc <= $rec{$key}){
				$check = 0;
			}
		}
		if($check){
			print OUT "$ent";
		}
		else{
			print OUT2 "$ent";
		}
	}
	elsif($mode eq 'either'){
		unless($pat=~/C/){
			my $check = 1;
			foreach my $key (sort {$a <=> $b} keys %rec){
				if($loc>=$key and $loc <= $rec{$key}){
					$check = 0;
				}
			}
			if($check){
				print OUT "$ent";
			}
			else{
				print OUT2 "$ent";
			}
		}
		else{
			print OUT2 "$ent";
		}
	}
	elsif($mode eq 'both'){
		unless($pat=~/C/){
			print OUT "$ent";
		}
		else{
			my $check = 1;
			foreach my $key (sort {$a <=> $b} keys %rec){
				if($loc>=$key and $loc <= $rec{$key}){
					$check = 0;
				}
			}
			if($check){
				print OUT "$ent";
			}
			else{
				print OUT2 "$ent";
			}
		}
	}
	else{
		print "Wrong term for mode! Only 'rec','multi','either' and 'both' acceptted!\n";
		exit;
	}
}
