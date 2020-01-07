#!/usr/bin/perl -w
my $in=shift;
my $out = shift;
unless (defined $in and defined $out){
	print "Usage:perl $0 <in> <out>\n";
	exit;
}
open IN,$in;
open OUT,">$out";
my $head = 1;
while(<IN>){
	if(/^##/){
		next;
	}
	if($head){
		my $f = $_;
		chomp $f;
		my @f = split /\t/,$f;
		my @names = @f[9..@f-1];
		my $names = join "\t",@names;
		print OUT "Location/Strain\t$names\n";
		$head = 0;
		next;
	}
	chomp;
	my @t = split /\t/;
	my $data = join "\t",@t[9..@t-1];
	my $ref = $t[3];
	my @pattern = split /,/,$t[4];
	$data =~ s/0/$ref/g;
	for (my $i = 1; $i<=@pattern;$i++){
		my $alt = $pattern[$i-1];
		$data =~s/$i/$alt/g;
	}
	print OUT "$t[1]\t$data\n";
}

