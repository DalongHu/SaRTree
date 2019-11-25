#!/usr/bin/perl -w
use Bio::SeqIO;
use FindBin qw($Bin);
use Getopt::Long;

# path of depandencies (Please check and revise them before using!)
my $mauveBin = 'progressiveMauve';
my $fasttreeBin = 'fasttreeMP';


my $input;
my $ref;
my $dir;
my $thread;
my $rec;
my $help;
GetOptions(
	"input|i=s"=>\$input,
	"ref|r=s"=>\$ref,
	"out|o=s"=>\$dir,
	"thread|t:i"=>\$thread,
	"recdetect|d!"=>\$rec,
	"help|h!"=>\$help
);
unless (defined $input and defined $ref and defined $dir){
	print "quicktree pipeline v1.1.\nTo build a tree quickly instead of parsnp before SaRTree V1.3 release.\nDependcies: fasttree, bioperl, mauve\nUsage: perl $0 -i|-input <in(dir with fastas, no reference)> -r|-ref <reference(fasta)> -o|-out <out(empty dir)> -t|-thread [thread (optional, number of cpus to use)] -d|-recdetect (optional, to open recdetect, default off, no parameter required)\nNote: Please check and revise the path of dependencies to fit your platform before using.\nWritten by Dalong Hu 21/Nov/2019.\nUpdated to v1.1 by Dalong Hu 25/Nov/2019\n";
	exit;
}
if(defined $thread){
        unless($thread > 0){ 
                print STDERR "thread should be a positive number\n";
                exit;
        }
}
unless(defined $rec){
	$rec = 0;
}


my $fasta = Bio::SeqIO->new(-file=>$ref,-format=>'fasta');
my $length = $fasta->next_seq->length;
my @samples = glob "$input/*";


#Mauve alignment
print "Alignment by Mauve\n";
system"mkdir $dir/mauve";
my $count_thread = 0;
my $cmd = '';
foreach my $sample(@samples){
	$sample=~/([^\/]+)\.([^\.]+)$/; # make sure using suffix 'fna' or 'fasta', and using strains' names as files' names
	my $name = $1;
	my $type = $2;
	if($type ne 'fna' and $type ne 'fasta' and $type ne 'fa'){
		next;
	}
	print "Starting alignment for $name\n";
	$cmd .= "$mauveBin --output=$dir/mauve/$name.mauve $ref $sample"; # will support HPC in formal version
	if(defined $thread){
		$cmd .= ' & ';
		$count_thread++;
		if($count_thread == $thread){
			`$cmd`;
			$count_thread = 0;
			$cmd = '';
		}
	}
	else{
		`$cmd`;
		$cmd = '';
		print "Finish alignment for $name\n";
	}
}
if($count_thread > 0){
	`$cmd`;
	$count_thread = 0;
	$cmd = '';
}

#Mauve to SNP list
print "Mauve to Map and Ins\n";
system"mkdir $dir/map";
@mauves = glob "$dir/mauve/*.mauve";
$count_thread = 0;
$cmd = '';
foreach my $mauve(@mauves){
	$mauve=~/([^\/]+)\.mauve/;
	my $name = $1;
	print "Converting format of $name\n";

	$cmd .= "perl $Bin/script/mauve2snplist.pl $mauve $dir/map/$name.snplist";
	if(defined $thread){
		$cmd .= ' & ';
		$count_thread++;
		if($count_thread == $thread){
			`$cmd`;
			$count_thread = 0;
			$cmd='';
		}
	}
	else{
		`$cmd`;
		$cmd = '';
		print "finished for $name\n";
	}
}
if($count_thread>0){
	`$cmd`;
	$count_thread = 0;
	$cmd = '';
}

# make whole SNP list
`perl $Bin/script/snp2reflist.pl $dir/map $dir/ref.list`;
`perl $Bin/script/snp2fake.pl $dir/map $dir/ref.list $dir/all`;
my @lists = sort {$a cmp $b} glob "$dir/map/*.extended";
my $list = join " ",@lists;
$list = "$dir/ref.list".' '.$list;
system "paste $list>$dir/all.list";
print "List Finished\n";
#recdetect
if($rec){
	print "RecDetect\n";
	system"perl $Bin/script/RecDetect.pl -s $dir/all.list -l $length -o $dir/rec";
	system"perl $Bin/script/List2Fake.pl $dir/rec.real $dir/final.fasta";
print "RecDetect finshed!\n";
}
else{
	`mv $dir/all.snp.fasta $dir/final.fasta`;
}

#tree construction
print "construction of tree using fasttreeMP\n";
if(defined $thread){
	`bash -c export OMP_NUM_THREADS=$thread`;
}
else{
	`bash -c export OMP_NUM_THREADS=1`;
}
system"$fasttreeBin -nt -gtr -gamma < $dir/final.fasta > $dir/tree.nwk";
print "Done\n";

system"perl $Bin/script/treeview.pl -in $dir/tree.nwk -a min -b";
