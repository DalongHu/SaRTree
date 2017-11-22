#!/usr/bin/perl
print "Installing SaRTree v1.0..\n";

#perl modules
print "Checking perl Packages..\n";
use File::Find;
my $re = shift || ".";
my %modules;
find sub {
	    return unless my ($x) = $File::Find::name =~ m{\./(.*\.pm)$};
	        $x =~ s,/,::,g;
		    $modules{$x}++ if /\Q$re\E/i;
}, map "$_/.", grep -d && /^[^.]/, @INC;

my @required = ('FindBin','Getopt::Long','Bio::SeqIO','Bio::TreeIO','Bio::AlignIO');
foreach $key (@required){
	if(defined $modules{"$key.pm"}){
		print "Module $key checked!\n";
	}
	else{
		while(1){
			print "Do you want to install $key now?(yes/no)\n";
			my $yes = <STDIN>;
			chomp $yes;
			if($yes eq 'y' or $yes eq 'yes'){
				print `cpan -i $key`;
				last;
			}
			elsif($yes eq 'n' or $yes eq 'no'){
				print "exit..\n";
				exit;
			}
		}
	}
}
print "perl Modules are checked!\n";

#mauve
print "Checking progressiveMauve and raxml..\n";
my $path2mauve = 'progressiveMauve';
my $mauve = `$path2mauve --version`;
while(!(defined $mauve)){
	print "Please install progressiveMauve and rerun this script to indicate the path.\n";
	print "Download it from http://darlinglab.org/mauve/download.html and follow the installation steps, or on Ubuntu system, please try 'sudo apt-get install progressivemauve' to install it.\n";
	print "Didn't find progressiveMauve, please indicate the path-to-progressiveMauve or input 'no' to exit:\n";
	print "example: /usr/bin/progressiveMauve\n";
	$path2mauve = <STDIN>;
	chomp $path2mauve;
	if($path2mauve eq 'no'){
		print "Please install progressiveMauve and rerun this script to indicate the path.\n";
		print "Download it from http://darlinglab.org/mauve/download.html and follow the installation steps, or on Ubuntu system, please try 'sudo apt-get install progressivemauve' to install it.\n";
		exit;
	}
	$mauve = `$path2mauve --version`;
}
print "progreassiveMauve has been found at\n $path2mauve\n";

#raxml
my $path2raxml = 'raxmlHPC';
my $raxml = `$path2raxml -v`;
while(!(defined $raxml)){
	print "Please install raxml and rerun this script to indicate the path.\n";
	print "Download it from https://github.com/stamatak/standard-RAxML and follow the installation steps, or on Ubuntu system, please try 'sudo apt-get install raxml' to install it.\n";
	print "Didn't find raxml, please indicate the path-to-raxml or input 'no' to exit:\n";
	print "example: /usr/bin/raxmlHPC or /usr/bin/raxml\n";
	$path2raxml = <STDIN>;
	chomp $path2raxml;
	if($path2raxml eq 'no'){
		print "Please install raxml and rerun this script to indicate the path.\n";
		print "Download it from https://github.com/stamatak/standard-RAxML and follow the installation steps, or on Ubuntu system, please try 'sudo apt-get install raxml' to install it.\n";
		exit;
	}
	$raxml = `$path2raxml -v`;
}
print "raxml has been found at\n $path2raxml\n";


#beast
while(1){
	print "Do you want to set path to BEAST?(optional)\n";
	print "input 'yes' or 'no' to indicate\n";
	my $yes = <STDIN>;
	chomp $yes;
	if($yes eq 'yes' or $yes eq 'y'){
		print "Setting BEAST!\n";
		last;
	}
	elsif($yes eq 'no' or $yes eq 'n'){
		goto LAST;
	}
}
my $path2beast = 'beast-mcmc';
my $beast = `$path2beast -version`;
while(!(defined $beast)){
	print "Please install BEAST v1.8.4 (Now SaRTree doesn't suppot BEAST v2.0 and later versions) and rerun this script to indicate the path.\n";
	print "Download it from http://beast.community/install_on_unix and follow the installation steps, or on Ubuntu system, please try 'sudo apt-get install beast-mcmc' to install it.\n";
	print "Didn't find BEAST, please indicate the path-to-BEAST or input 'no' to exit:\n";
	print "example: /usr/bin/beast-mcmc\n";
	$path2beast = <STDIN>;
	chomp $path2beast;
	if($path2beast eq 'no'){
		print "Please install BEAST v1.8.4 and rerun this script to indicate the path.\n";
		print "Download it from http://beast.community/install_on_unix and follow the installation steps, or on Ubuntu system, please try 'sudo apt-get install raxml' to install it.\n";
		exit;
	}
	$beast = `$path2beast -v`;
}
print "BEAST has been found at\n $path2beast\n";

LAST:
`mv Config Config.backup`;
open COF,"template/config.template";
open OUT,">Config";
while(<COF>){
	if(/^mauveBin/){
		print OUT "mauveBin\t$path2mauve\n";
		next;
	}
	if(/^raxMLBin/){
		print OUT "raxMLBin\t$path2raxml\n";
		next;
	}
	if(/^beastBin/){
		if(defined $beast){
			print OUT "beastBin\t$path2beast";
			next;
		}
		else{
			print OUT $_;
			next;
		}
	}
	print OUT $_;
}
print "Installation finished! All the pathes and settings are stored in Config file in the SaRTree main folder. By revising that Config file could change the settings later. Please run example.sh to test\n";
