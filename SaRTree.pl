#!/usr/bin/perl
use Bio::SeqIO;
use FindBin qw($Bin);
use Getopt::Long;

## parameter definations

$dir; # working dir
$tag; # output tag
$ref; # reference sequence
$gff; # reference gff file for annotation
$date; # date list of isolates
$outgroup; # name of outgroup
$model; # running models
$beast; # run beast automatically
$help; # to show help message
$version; # to show version
$userTree; # using existing tree / for locater
$snp; # SNP list for existing tree for locater
$root; # method to root the tree
$locater; # run strainLocater
$thread;
$coverage;
$config;
$filter;

GetOptions (
	"dir|i=s"=>\$dir,
	"tag|t=s"=>\$tag, # just a tag, not to add path here
	"ref|f=s"=>\$ref,
	"gff|g=s"=>\$gff, # now mandatory, will be optional in next version
	"date|d:s"=>\$date,
	"outgroup|o:s"=>\$outgroup,
	"model|m=s"=>\$model,
	"beast|b:s"=>\$beast,
	"userTree|u:s"=>\$userTree,
	"snp|s:s"=>\$snp,
	"root|r:s"=>\$root,
	"locater|l!"=>\$locater,
	"thread|a:i"=>\$thread,
	"coverage|c:i"=>\$coverage,
	"config|p:s"=>\$config,
	"filter|e:s"=>\$filter,
	"help|h!"=>\$help,
	"version|v!"=>\$version
);# will add multi chromosome support


#------------------------------------------------------------------#
## default paths and parameters should be revised by user
open DEF,"$Bin/Config";
my $default = '';
while(<DEF>){
	chomp;
	my @tmp_default = split /\t/;
	my $default_key = $tmp_default[0];
	my $default_value = $tmp_default[1];
	if($default_value ne '-'){
		$$default_key = $default_value;
	}
	else{
		$$default_key = '';
	}
	$default .= "$default_key\t$$default_key\n";
}
if($thread){
	        $raxMLConfig .= " -T $thread";
		        $beastConfig .= " -threads $thread";
}

#------------------------------------------------------------------#

@config = ('dir','tag','ref','gff','date','outgroup','model','beast','userTree','snp','root','locater','thread','coverage','filter','mauveBin','mauveConfig','raxMLBin','raxMLConfig','raxMLBootstrap','beastBin','beastConfig','treeAnnotatorConfig','recDetectArg');
my %config;
foreach my $config_key(@config){
	if(defined $$config_key){
		$config{$config_key} = $$config_key;
	}
	else{
		$config{$config_key} = '-';
	}
}
## config phrasing
if(defined $config){
	open COF,$config;
	while (<COF>){
		chomp;
		my @t = split /\t/;
		if(defined $config{$t[0]}){
			if($t[1] ne '-'){
				my $tmp_key = $t[0];
				$$tmp_key = $t[1];
			}
		}
		else{
			my $tmp_key = $t[0];
			print STDERR "no parameter $tmp_key acceptable!\n";
			exit;
		}
	}
	close COF;

}
foreach my $config_key(@config){
	if(defined $$config_key){
		$config{$config_key} = $$config_key;
	}
}


$usage = "\n-------------------SaRTree v1.1-------------------\nUsage:\n
==============SaRTree core==============\n
1/ make a working directory (e.g. mkdir test)\n\n
2/ (optional) put reference sequence(fasta file, complete genome, 1 chromosome) and its annotation(gff file) into working directory (use <strainname>.fna and <strainname>.gff as the format of file names, e.g. AB0057.fna and AB0057.gff) (will support multi chromosomes in formal version)\n\n
3/ choose a model which depends on the input files(\"standard\" for mauve result, \"raw\" for raw sequence data, \"formatted\" for SaRTree format mapping or snp list, \"root\" for continuing existing run interupted by rooting the tree)\n\n
for \"standard\" model (the recommended running model)\nmake a directory named \"mauve\"(without the quotes) in the working directory (temporarily only mauve result rupport)\nput all pairwise aligment result by Mauve(btw ref and each strain) into the directory (use <strainname>.mauve as the format of file name, e.g. ZJ06.mauve)\n\n
for \"raw\" model (run mauve in the script, will support multi-threads in formal version)\nmake a directory named \"seq\"(without the quotes) in the working directory\nput all raw samples' sequences (fasta format, use <strainname>.fna as the format of file name, e.g. ZJ06.fna) into the directory\n\n
for \"formatted\" model (not recommended, just for test)\nput SaRTree fotmat snp list file (a tab separated file named as <tag>.list,the \"tag\" should be same as the one using for the -t option, with the format:\nlocation\tref\t-s1-\t-s2-\t...\n1234\tA\tT\tT\t...\n2345\tT\tT\tA\t...\n) and the name list file (tab separated txt file, with format:\nref\tAB0057\n-s1-\tZJ06\n-s2-\t...\n...\n) into working directory\n\n
for \"root\" model\nuse the same commond as the first run and add -u to indicate the rooted tree.\n\n
for \"beast\" model\nuse the same commond as the first run and add -d to indicate the isolation dates and add -u to indicate the nex format output tree from TreeAnnotator.\n\n
4/ run the main script:\nperl $0 [options] -i <work dir> -t <output tag> -f <ref sequence(fna)> -g <ref annotation(gff)> -m <standard/raw/formatted/root/beast> [options]\n\n
==============StrainLocater==============\n
1/ make a working directory (e.g. mkdir test) and make a directory named \"query\" in the working directory\n\n
2/ (optional) put reference sequence(fasta file, complete genome, 1 chromosome) into working directory (use <strainname>.fna as the format of file name, e.g. AB0057.fna) (will support multi chromosomes in formal version)\n\n
3/ choose a model which depends on the input files(\"standard\" for mauve result, \"raw\" for raw sequence data, \"formatted\" for SaRTree format mapping or snp list)\n\n
the details of models are similar to the \"SaRTree core\" section above, but NOTE that all the input files should be put into <work dir>/query/seq (\"raw\" model), <work dir>/query/mauve (\"standard\" model)or <work dir>/query/ (\"formatted\" model)\n
Adding a \"query\" directory separately is to make it possible to run SaRTree core and StrainLocater in the same working directory conveniently\n\n
4/ run the main script: (indicating newick format target tree file and SaRTree format, see the \"formatted\" model section, snp list file for the tree)\n
perl SaRTree.pl [options] -l -i <work dir> -t <output tag> -f <ref sequence(fna)> -u <target tree(nwk)> -s <snp list(SaRTree format)> -m <standard/raw/formatted>\n\n

==============treeview.pl==============

a small script to show a tree on commond line text interface, for details, please run perl script/treeview.pl -h

==============Options==============
	-dir/-i <string>  mandatory, FULL(!) path to the working directory, that is a problem of RAxML, will revise that in next version\n
	-tag/-t <string>  mandatory, output tag, just a tag, not to add path here\n
	-ref/-f <string>  mandatory, full path to reference sequence file (fasta format, named as <strainname>.fna, e.g. AB0057.fna)\n
	-gff/-g <string>  mandatory, full path to annotation file of reference strain (gff format, named as <strainname>.gff, e.g. AB0057.gff)\n
	-model/-m <standard/raw/formatted/root>  mandatory, three input models could be selected depending on the input files\n
	-locater/-l  optional, to run StrainLocater instead of SaRTree main process, when using -l, a target tree and its snp list should also be specified by -u and -s respectively, and a directory named \"query\" should be built in working directory then all input/output files for StrainLocater will be in the directory\n
	-userTree/-u [string]  optional, a newick format tree file, when -l used, -u specifies the target tree to locate. when without -l, -u specifies a user's own phylogenetic tree instead of building a tree in SaRTree process\n
	-snp/-s [string]  optional, when -l used, -s specifies the snp list of the target tree\n
	-root/-r [auto/man/off] optional, the method to root the tree, default \"auto\"; \"auto\": root the tree automatically by StrainLocater based on events (only this method could locate events on the two descendant branches of root correctly); \"man\": program stop after tree building, then user could root the tree manually, but then should restart SaRTree again by indicating \"root\" model as use -u to indicate the rooted tree, and NOTE that the events on two descendant branches of root may not be right, all of the events of them will be counted to \'the branch including reference\'; \"off\": just use the output tree of tree building software, NOT recommended\n 
	-outgroup/-o [string]  optional, outgroup strain's name could be specified by -o; when using \"auto\" model for -r, outgroup will be used by StrainLocater;  when using \"off\" or \"man\" for -r, outgroup will be used by raxml (NOT recommended, should root the tree carefully)\n
	-beast/-b [off/auto/man]  optional, run beast to estimate divergence time, a file with isolattion date should be specified by -d, NOT recommended now, there are some bugs, default off\n
	-date/-d [string]  optional, when -b indicated, -d specofies the isolattion date file with the format(tab separated):\n\tstrain1\t1970\n\tstrain2\t1980\n\t...\n
	-coverage/-c [integer]  optional, an integer btw 0-100, coverage to reference filter, all samples with lower coverage to reference will be put into a \"query\" folder for second run using strainlocater, default 0\n
	-filter/-e [off/rec/multi/both/either] optional, filter out SNPs with low confidence, \"off\" indicates close this function;\"rec\" mode deletes SNPs covered by recombination regions ignroing recombinant SNPs' distribution;\"multi\" mode deletes SNPs presenting more than two types of bases;\"both\" mode deletes only SNPs fitting both feathers detected by \"rec\" and \"multi\" modes;\"either\"  mode deletes SNPs fitting all SNPs fitting either \"rec\" and \"multi\" modes. default \"off\"\n
	-config/-p [string]  optional, a config file with parameters, please use '-' or leave the line out to indicate a default or unavailble parameter. See the config.template.txt in the SaRTree folder for details. When using the config file, all parameters in the file will cover any other input parameters including the commond line input and any default.\n
	-thread/-a [integer]  optional, number of threads to use in multi-threads model.
	-help/-h  optional, show this usage\n
	-version/-v  optional, show version information\n

For details, see README file.\n
Default setting:\n$default\n
\n";

$v_information = "------------------------------------------\nSaRTree is for locating evolutionary events onto phylogenetic trees or building high resolution trees based on mutation events, StrainLocater is for locating new strains onto existing trees.\n\nPlease note that this is a debug test version v1.1, so not all of the functions work well now and many parts of the codes run slowly and redundantly.\n\nThere is no warranties coming with SaRTree, so that users must be responsible for the results generated by SaRTree.\n\nThis is a free software following GNU General Public License <http://www.gnu.org/licenses/gpl.html>\nCopyright: Dalong Hu, University of Sydney (dalong.hu\@sydney.edu.au) Sep-2017\n------------------------------------------\n";





##trouble shooting


if($help){
	print + $usage.$v_information;
	exit;
}
if($version){
	print + $v_information;
	exit;
}
unless($locater){
	unless($model and $dir and $tag and $ref and $gff){
		print STDERR "Missing mandatory options!\n";
		print STDERR $usage.$v_information;
		print "$model\n$dir\n$tag\n$ref\n$gff\n";
		exit;
	}
}
else{
	unless($model and $dir and $tag and $ref){
		print STDERR "Missing mandatory options!\n";
		print STDERR $usage.$v_information;
		exit;
	}
	unless ($snp and $userTree){
		print STDERR "-snp and -userTree are mandatory  when running locater model!\n";
		print STDERR $usage.$v_information;
		exit;
	}
	if($model eq 'root'){
		print STDERR "StrainLocater doesn't support \'root\' model\n";
		exit;
	}
	if(defined $outgroup){
		print STDERR "StrainLocater doesn't support outgroup\n";
		exit;
	}
	$root = 'off';
}
if($model ne 'raw' and $model ne 'formatted' and $model ne 'standard' and $model ne 'root' and $model ne 'beast'){
	print STDERR "models to choose: raw, formatted and standard (or root for second run to rooting,beast for second run to run BEAST)\n";
	print STDERR $usage.$v_information;
	exit;
}

unless(defined $root){
	$root = 'auto';
}
if($root ne 'auto' and $root ne 'man' and $root ne 'off'){
	print STDERR "rooting models to choose: auto, man and off\n";
	exit;
}
if($root eq 'auto'){
	unless(defined $outgroup){
		print STDERR "outgroup needed when using \"auto\" rooting! could indicate it by -o or use \"off\" or \"man\" for -r to ignore outgroup\n";
		exit;
	}
	if ($model eq 'formatted'){
		print STDERR "formatted model only for testing, doesn't support auto rooting\n";
		exit;
	}
}
if($root eq 'man'){
	if(defined $outgroup){
		print STDERR "don't indicate outgroup when using \"man\" for -r\n";
		exit;
	}
}

if($model eq 'root'){
	unless(defined $userTree){
		print STDERR "in root model, should indicate rooted tree by -u\n";
		exit;
	}
}
if($model eq 'beast'){
	unless(defined $date and defined $userTree){
		print STDERR "in beast model, should indicate isolation date file by -date and TreeAnnoator output nex file by -u\n";
	}
}

unless(defined $beast){
	$beast = 'off';
}
else{
	if($beast ne 'off' and $beast ne 'auto' and $beast ne 'man'){
		print STDERR "beast models to choose: auto, man and off\n";
		exit;
	}
}

if($beast eq 'auto' xor $date){
	print STDERR "-beast auto and -date should be indicated at the same time!\n";
	exit;
}

if(defined $coverage){
	unless($coverage =>0 and $coverage <=100){
		print STDERR "-coverage shoule be btw 0 and 100!\n";
		exit;
	}
}
if(defined $thread){
	unless($thread > 0){
		print STDERR "thread should be a positive number";
		exit;
	}
}

if(defined $filter){
	if($filter ne 'multi' and $filter ne 'rec' and $filter ne 'either' and $filter ne 'both' and $filter ne 'off'){
		print STDERR "filter modes are one of rec,multi,either,both and off";
		exit;
	}
}
else{
	$filter = 'off'; # default off, 'both' is recommended to filter SNPs both with multi-type bases and included in recombination regions
}


sub name2Index{
	my $name = $_[0];
	my $nameFile = $_[1];
	my %nameFile;
	open NAMEFILE,$nameFile;
	while(<NAMEFILE>){
		chomp;
		my @tmp = split;
		$nameFile{$tmp[1]} = $tmp[0];
	}
	my $index = $nameFile{$name};
	return $index;
}


## main script strats:

my $fasta = Bio::SeqIO->new(-file=>$ref,-format=>'fasta');
my $length = $fasta->next_seq->length;
my $tmp_tree;
my $out_tree;

if($locater){
	$dir = $dir."\/query";	
}

if($model eq 'root'){
	system "mv $dir/result/$tag.name.nwk $dir/tmp/$tag.name.unrooted.nwk";
	goto ROOT;
}
if($model eq 'beast'){
	system "rm $dir/result/$tag.fake.fasta";
	system "cp $userTree $dir/tmp/$tag.date.nex";
	system "perl $Bin/script/Name2Index.pl $dir/tmp/$tag.date.nex $dir/tmp/$tag.name.txt $dir/tmp/$tag.tmp.nex";
	$tmp_tree = "$dir/tmp/$tag.tmptree.nwk";
	$out_tree = "$dir/result/$tag.name.nwk";
	goto BEAST;
}


system"mkdir $dir/tmp";
system"mkdir $dir/result";
open OUTC,">$dir/result/$tag.config.txt";
foreach my $config_key(@config){
	if(defined $$config_key and $$config_key ne ''){
		print OUTC "$config_key\t$$config_key\n";
	}
	else{
		print OUTC "$config_key\t-\n";
	}
}
#open LOG, ">$dir/log.txt";



#MauveScript - 'raw' model, will include HPC support
if ($model eq 'raw'){
	my @samples = glob "$dir/seq/*.fna"; # make a seq dir and put all input fna into
	print "Alignment by Mauve\n";
	system"mkdir $dir/mauve";
	my $count_thread = 0;
	my $cmd = '';
	foreach my $sample(@samples){
		$sample=~/([^\/]+)\.fna/; # make sure using suffix 'fna', and using strains' names as files' names
		$name = $1;
		print "Starting alignment for $name\n";
		$cmd .= "$mauveBin $mauveConfig --output=$dir/mauve/$name.mauve $ref $sample"; # will support HPC in formal version
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
}

if ($model ne 'formatted'){
#Mauve2Map
	print "Mauve to Map and Ins\n";
	system"mkdir $dir/map";
	@mauves = glob "$dir/mauve/*.mauve";
	my $count_thread = 0;
	my $cmd = '';
	foreach my $mauve(@mauves){
		$mauve=~/([^\/]+)\.mauve/;
		my $name = $1;
		print "Converting format of $name\n";
		
		$cmd .= "perl $Bin/script/Mauve2Map.pl $mauve $dir/map/$name.map"; # will support HPC in formal version
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

#Map2List
	print "Map to List\n";
	system"perl $Bin/script/MakeRefList.pl $ref $dir/tmp $dir/tmp/$tag.name.txt";
	system"mkdir $dir/list";
	if($root eq 'auto'){ # for auto rooting
		system "mv $dir/map/$outgroup.map* $dir/tmp/";
		system"perl $Bin/script/Map2List.pl $dir/tmp/$outgroup.map $length $dir/tmp/$outgroup.list outgroup $dir/tmp/outgroup.name.txt";
	}
	@maps = glob "$dir/map/*.map";
	my $index = 1;
	$count_thread = 0;
	$cmd = '';
	foreach my $map(sort {$a cmp $b} @maps){
		$map=~/([^\/]+)\.map/;
		my $name = $1;
		print "Coverting format of $name\n";
		$cmd .="perl $Bin/script/Map2List.pl $map $length $dir/list/$name.list -s$index- $dir/tmp/$tag.name.txt";
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
			print "finished for $name\n";
		}
		$index++;
	}
	if($count_thread>0){
		`$cmd`;
		$count_thread = 0;
		$cmd = '';
	}

#coverage filter (haven't implement multi threads!)
	if(defined $coverage){
		system "mkdir $dir/query";
		system "mkdir $dir/query/list";
		system "mkdir $dir/query/map";
		system "perl $Bin/script/coverage.pl $dir/list $dir/tmp/coverage.sum.txt";
		open COV,"$dir/tmp/coverage.sum.txt";
		open COT,">$dir/tmp/coverage.out.txt";
		while(<COV>){
			chomp;
			my @t = split /\t/;
			if($t[2]<$coverage){
				system "mv $dir/list/$t[1].list $dir/query/list/";
				system "mv $dir/map/$t[1].map* $dir/query/map/";
				print COT "$t[1]\t$t[2]\n";
			}
		}

	}


#Make Whole List
	print "Make Whole List\n";
	if($root eq 'auto'){ #for auto rooting
		system "paste $dir/tmp/ref.list $dir/tmp/$outgroup.list > $dir/tmp/outgroup.list";
		system"perl $Bin/script/List2SNP.pl $dir/tmp/outgroup.list $dir/tmp/outgroup.snp";
		system"perl $Bin/script/Index2Name-title.pl $dir/tmp/outgroup.snp $dir/tmp/outgroup.name.txt $dir/tmp/outgroup.name.snp";
	}
	my @lists = sort {$a cmp $b} glob "$dir/list/*.list";
	my $list = join " ",@lists;
	$list = "$dir/tmp/ref.list".' '.$list;
	system "paste $list>$dir/tmp/$tag.list";
	print "List Finished\n";
}
else{# need to provide mapping list or snp list with 'ref' '-s1-'... names and a name list.
	system"cp $dir/$tag.name.txt $dir/tmp/$tag.name.txt";
	system"cp $dir/$tag.list $dir/tmp/$tag.list"; # will supprot varient formats to save run time
}

system"head -n 1 $dir/tmp/$tag.list > $dir/tmp/order.list";
system"mv $dir/tmp/$tag.name.txt $dir/tmp/$tag.name.old.txt";
system"perl $Bin/script/orderNameList.pl $dir/tmp/$tag.name.old.txt $dir/tmp/order.list $dir/tmp/$tag.name.txt";

if ($locater){
	print "Starting StrainLocater!\n";
	system"perl $Bin/script/List2SNP.pl $dir/tmp/$tag.list $dir/tmp/$tag.snp";
	system"perl $Bin/script/Index2Name-title.pl $dir/tmp/$tag.snp $dir/tmp/$tag.name.txt $dir/tmp/$tag.name.snp";
	system"perl $Bin/script/strainLocater.pl -l $dir/tmp/$tag.name.snp -s $snp -t $userTree -o $dir/result/$tag";
	print "StrainLocater finished!\n";
#	system"rm -r $dir/tmp $dir/mauve $dir/map $dir/list"; # for debugging, no temp file deleted
	exit;
}

#RecDetect will support other softwares
print "RecDetect\n";
system"perl $Bin/script/List2Match.pl $dir/tmp/$tag.list $dir/tmp/$tag.match";
system"perl $Bin/script/List2SNP.pl $dir/tmp/$tag.match $dir/tmp/$tag.snp";
system"perl $Bin/script/RecDetect.pl -s $dir/tmp/$tag.snp -l $length -o $dir/tmp/$tag $recDetectArg";
print "RecDetect finshed!\n";

#Filter out low confidence SNPs
print "Filtering low condidence SNPs...\n";
unless($filter eq 'off'){
	system"mv $dir/tmp/$tag.real $dir/tmp/$tag.real.tmp";
	system"perl $Bin/script/FiltComplex.pl $dir/tmp/$tag.real.tmp $dir/tmp/$tag.rec $filter $dir/tmp/$tag.real";
}



ROOT: # the label to restart SaRTree for manually rooting model

unless($userTree){
#tree building
	print "Building phylogenetic tree!\n";
	system"perl $Bin/script/List2Fake.pl $dir/tmp/$tag.real $dir/tmp/$tag.fake.fasta";
	system"perl $Bin/script/fasta2phylip.pl $dir/tmp/$tag.fake.fasta $dir/tmp/$tag.fake.phy";
	if($root ne 'auto'){
		if(defined $outgroup){
			$outgroup = name2Index($outgroup,"$dir/tmp/$tag.name.txt");
			$raxMLConfig.=" -o $outgroup";
		}
	}
	system"$raxMLBin $raxMLConfig -s $dir/tmp/$tag.fake.phy -n target.tree.nwk -w $dir/tmp";
	system"$raxMLBin $raxMLBootstrap $raxMLConfig -s $dir/tmp/$tag.fake.phy -n bootstrap.tree.nwk -w $dir/tmp";
	system"$raxMLBin -f b $raxMLConfig -t $dir/tmp/RAxML_bestTree.target.tree.nwk -z $dir/tmp/RAxML_bootstrap.bootstrap.tree.nwk -n final.tree.nwk -w $dir/tmp";
	system"perl $Bin/script/Index2Name.pl $dir/tmp/RAxML_bipartitions.final.tree.nwk $dir/tmp/$tag.name.txt $dir/result/$tag.name.nwk";
	$tmp_tree = "$dir/tmp/RAxML_bipartitions.final.tree.nwk";
	$out_tree = "$dir/result/$tag.name.nwk";
	print "Tree finished\n";
	if($root eq 'man'){
		print "SaRTree stops for tree rooting manually\nOutput tree is saved in $dir/result/$tag.name.nwk, please root it and rerun SaRTree with \"root\" for model and indicate rooted tree file by -u\n";
		exit;
	}
	if($root eq 'auto'){
		print "tree rooting\n";
		system"perl $Bin/script/Index2Name-title.pl $dir/tmp/$tag.real $dir/tmp/$tag.name.txt $dir/tmp/$tag.name.real";
		system"mv $dir/result/$tag.name.nwk $dir/tmp/$tag.name.unrooted.nwk";
		system"perl $Bin/script/strainLocater.pl -l $dir/tmp/outgroup.name.snp -s $dir/tmp/$tag.name.real -t $dir/tmp/$tag.name.unrooted.nwk -o $dir/tmp/outgroup";
		my $rootingresult = `perl $Bin/script/rooting.pl $dir/tmp/$tag.name.unrooted.nwk $dir/tmp/outgroup.onTree.txt $dir/result/$tag.name`; #need to implement!! and also the information of the two branches to root
		if ($rootingresult =~ /fail/){
			print STDERR "fail to locate outgroup onto the tree, please root the tree manually and use root model to rerun SaRTree\n";
			exit;
		}
		system"mv $dir/result/$tag.name.tmp.* $dir/tmp/";
		system"perl $Bin/script/Name2Index.pl $out_tree $dir/tmp/$tag.name.txt $dir/tmp/$tag.root.id.nwk";
		$tmp_tree = "$dir/tmp/$tag.root.id.nwk";
		print "\nUnrooted tree:\n";
		system "perl $Bin/script/treeview.pl -in $dir/tmp/$tag.name.unrooted.nwk -a min -b";
		print "\nRooted tree:\n";
		system "perl $Bin/script/treeview.pl -in $dir/result/$tag.name.nwk -a min -b";
		print "\nTree rooting finished!\n";
	}
}
else{
	system"perl $Bin/script/Name2Index.pl $userTree $dir/tmp/$tag.name.txt $dir/tmp/usertree.id.nwk";
	$tmp_tree = "$dir/tmp/usertree.id.nwk";
	system "cp $userTree $dir/result/$tag.name.nwk";
	$out_tree = "$dir/result/$tag.name.nwk";
}

system "cp $tmp_tree $dir/tmp/$tag.tmptree.nwk";
$tmp_tree ="$dir/tmp/$tag.tmptree.nwk";


if ($beast eq 'auto'){
#Beast
	print "Building tree with divergence date\n";
	system"perl $Bin/script/Name2Index-list.pl $date $dir/tmp/$tag.name.txt $dir/tmp/$tag.date.txt"; ###havn't done.!!!!!!!!!!!!!!!1
	system"perl $Bin/script/MakeSimpleXml.pl $dir/tmp/$tag.fake.fasta $dir/tmp/$tag.date.txt $dir/tmp/$tag.xml";
	system"$beastBin $beastConfig $dir/tmp/$tag.xml";
	$beastBin=~/(.*)\/[^\/]+/;
	my $beastBinPath =$1;
	system"$beastBinPath/treeannotator $treeAnnotatorConfig -target $tmp_tree $dir/tmp/$tag.fake.fasta.trees $dir/tmp/$tag.tmp.nex";
	system"perl $Bin/script/Index2Name.pl $dir/tmp/$tag.tmp.nex $dir/tmp/$tag.name.txt $dir/result/$tag.date.nex";
	print "Beast finished!\n";
}
elsif($beast eq 'man'){
	system"perl $Bin/script/Index2Name.pl $dir/tmp/$tag.fake.fasta $dir/tmp/$tag.name.txt $dir/result/$tag.fake.fasta";
	print "SaRTree stops for running BEAST manually\nFake multialignment is saved in $dir/result/$tag.fake.fasta, target tree is saved in $dir/result/$tag.name.nwk please use them to run BEAST and rerun SaRTree with \"beast\" for model and indicate isolation date file by -date and output tree with divergence dates generated by TreeAnnotator in nex format by -u\n";
	exit;
}

BEAST: # the label to restart SaRTree for manually BEAST model



#Summary and Annotation of Result
print "Annotation of Result!\n";
if($root ne 'auto'){
	system"perl $Bin/script/Index2Name-title.pl $dir/tmp/$tag.real $dir/tmp/$tag.name.txt $dir/tmp/$tag.name.real";
	system"perl $Bin/script/getRoot.pl $out_tree $dir/tmp/$tag.name.real $dir/tmp/outgroup.onTree.txt";
}
else{
	open TEST,"$dir/tmp/outgroup.onTree.txt";
	my $testLocationFile = <TEST>;
	$testLocationFile = <TEST>;
	unless(defined $testLocationFile){
		system"perl $Bin/script/Index2Name-title.pl $dir/tmp/$tag.real $dir/tmp/$tag.name.txt $dir/tmp/$tag.name.real";
		system"perl $Bin/script/getRoot.pl $out_tree $dir/tmp/$tag.name.real $dir/tmp/outgroup.onTree.txt";
	}
}
system"perl $Bin/script/annotateSNP.pl $dir/tmp/$tag.real $gff $ref $dir/tmp/outgroup.onTree.txt $dir/tmp/$tag.snp.annotated.xls";
system"perl $Bin/script/annotateRec.pl $dir/tmp/$tag.rec $gff $length $dir/tmp/outgroup.onTree.txt $dir/tmp/$tag.rec.annotated.xls";
system"perl $Bin/script/sumInfo.pl $dir/tmp/$tag.snp.annotated.xls $dir/tmp/$tag.snp.sum";
system"perl $Bin/script/sumInfo.pl $dir/tmp/$tag.rec.annotated.xls $dir/tmp/$tag.rec.sum";
if($model ne 'formatted'){
	system"perl $Bin/script/Ins2List.pl $dir/map $ref $dir/tmp/$tag.name.txt $dir/tmp/$tag.ins.sum";
	system"perl $Bin/script/List2N.pl $dir/tmp/$tag.list $dir/tmp/$tag.N.list";
	system"perl $Bin/script/List2NOMAP.pl $dir/tmp/$tag.list $dir/tmp/$tag.NOMAP.list";
	system"perl $Bin/script/List2Del.pl $dir/tmp/$tag.list $dir/tmp/$tag.del.tmp.list";
	system"perl $Bin/script/sumDel.pl $dir/tmp/$tag.del.tmp.list $dir/tmp/$tag.del.snp.list $dir/tmp/$tag.del.sum";
	system"cat $dir/tmp/$tag.ins.sum $dir/tmp/$tag.del.sum > $dir/tmp/$tag.indel.seq.xls";
	system"perl $Bin/script/Index2Name-title.pl $dir/tmp/$tag.indel.seq.xls $dir/tmp/$tag.name.txt $dir/result/$tag.indel.seq.name.xls";
	system"perl $Bin/script/countIndel.pl $dir/tmp/$tag.indel.seq.xls $dir/tmp/$tag.indel.xls";
	system"perl $Bin/script/annotateIndel.pl $dir/tmp/$tag.indel.xls $gff $length $dir/tmp/outgroup.onTree.txt $dir/tmp/$tag.indel.annotated.xls";
	system"perl $Bin/script/sumInfo.pl $dir/tmp/$tag.indel.annotated.xls $dir/tmp/$tag.indel.sum";
	system"perl $Bin/script/tree_annotate_indel.pl $tmp_tree $dir/tmp/$tag.event.nex $dir/tmp/$tag.snp.sum $dir/tmp/$tag.rec.sum $dir/tmp/$tag.indel.sum $dir/tmp/outgroup.onTree.txt $dir/tmp/$tag.name.txt";
	if($beast ne 'off'){
		system"perl $Bin/script/sumBranchWithDate.pl $out_tree $dir/tmp/$tag.name.nex $dir/result/branch.xls $dir/tmp/$tag.snp.annotated.xls $dir/tmp/$tag.rec.annotated.xls $dir/tmp/$tag.indel.annotated.xls $dir/tmp/outgroup.onTree.txt $dir/tmp/$tag.name.txt $date";
	}
	else{
		system"perl $Bin/script/sumBranch.pl $out_tree $dir/result/branch.xls $dir/tmp/$tag.snp.annotated.xls $dir/tmp/$tag.rec.annotated.xls $dir/tmp/$tag.indel.annotated.xls $dir/tmp/outgroup.onTree.txt $dir/tmp/$tag.name.txt";
	}
	system"perl $Bin/script/finalXls.pl $dir/tmp/$tag.indel.annotated.xls $dir/result/branch.xls.pattern2branch.txt $dir/tmp/$tag.name.txt $dir/result/$tag.indel.final.xls";
}
else{
	system"perl $Bin/script/tree_annotate.pl $tmp_tree $dir/tmp/$tag.event.nex $dir/tmp/$tag.snp.sum $dir/tmp/$tag.rec.sum $dir/tmp/outgroup.onTree.txt $dir/tmp/$tag.name.txt";
	system"perl $Bin/script/sumBranchNoIndel.pl $out_tree $dir/result/branch.xls $dir/tmp/$tag.snp.annotated.xls $dir/tmp/$tag.rec.annotated.xls $dir/tmp/outgroup.onTree.txt $dir/tmp/$tag.name.txt";
}
system"perl $Bin/script/Index2Name.pl $dir/tmp/$tag.event.nex $dir/tmp/$tag.name.txt $dir/result/$tag.event.name.nex";
system"perl $Bin/script/finalXls.pl $dir/tmp/$tag.snp.annotated.xls $dir/result/branch.xls.pattern2branch.txt $dir/tmp/$tag.name.txt $dir/result/$tag.snp.final.xls";
system"perl $Bin/script/finalXls.pl $dir/tmp/$tag.rec.annotated.xls $dir/result/branch.xls.pattern2branch.txt $dir/tmp/$tag.name.txt $dir/result/$tag.rec.final.xls";

#system"rm -r $dir/tmp $dir/mauve $dir/map $dir/list"; # for debugging, no temp file deleted

print "All finished!\n";
