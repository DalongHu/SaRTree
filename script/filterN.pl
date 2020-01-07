#!/usr/bin/perl -w
my $in = shift;
my $cut = shift;
my $out = shift;
unless(defined $in and defined $cut and defined $out){
	print "Usage:perl $0 <in> <cut> <out>\n";
	exit;
}
open IN,$in;
open OUT,">$out";
my $f = <IN>;
print OUT "$f";
$cut = $cut/100.0;
while(<IN>){
	chomp;
	my @t = split /\t/;
	@t = @t[2..@t-1];
	my $count =0;
	foreach my $t (@t){
		if($t eq 'N'){
			$count++;
		}
	}
	my $pro = $count/@t;
	if($pro <= $cut){
		print OUT "$_\n";
	}
	else{
		print "$count\t$pro\t$cut\t$_\n";
	}

}
