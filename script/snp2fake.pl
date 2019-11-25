#!/usr/bin/perl -w
my $in =shift;
my $ref = shift;
my $out = shift;
unless (defined $in and defined $ref and defined $out){
	print "Usage:perl $0 <map(dir)> <reflist> <out>\n";
	exit;
}

my %ref;
open REF,$ref;
my $fref = <REF>;
while(<REF>){
	chomp;
	my @t = split /\t/;
	$ref{$t[0]}= $t[1];
}
close REF;

my @maps = glob"$in/*.snplist";
open OUT,">$out.snp.fasta";
my @seqs = sort {$a <=> $b} keys %ref;
print OUT ">ref\n";
my $refseq = '';
foreach my $loc (@seqs){
	$refseq .= $ref{$loc};
}
print OUT "$refseq\n";

foreach my $map (@maps){
        open IN,"$map";
        my $id = '';
        if($map =~/([^\/]+)\.snplist/){
                $id = $1;
        }
        else{
                print "error on the name of $map\n";
        }
	my %snp;
	my %gap;
	while(<IN>){
		chomp;
		my @t = split /\t/;
		$snp{$t[0]} = $t[2];
	}
	open GAP,"$map.gap";
	my $gaps = <GAP>;
	chomp $gaps;
	my @gaps = split /\t/,$gaps;
	foreach my $gap (@gaps){
		$gap{$gap}++;
	}
	my $snpseq = '';
	my $snplist = '';
	open LST,">$map.extended";
	foreach my $loc (@seqs){
		if(defined $snp{$loc}){
			$snpseq .= $snp{$loc};
			$snplist .= "$snp{$loc}\n";
		}
		elsif(defined $gap{$loc}){
			$snpseq .= 'N';
			$snplist .= "N\n";
		}
		else{
			$snpseq .= $ref{$loc};
			$snplist .= "$ref{$loc}\n";
		}
	}
	print LST "$id\n$snplist";
	print OUT ">$id\n$snpseq\n";

}
