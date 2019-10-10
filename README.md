------------
##SaRTree v1.2
------------

##Description:
------------------------------------------

###==============Main problems to solve==============

1 At present, a huge amount of phylogenetic analysis are published in variant areas. However, very few phylogenetic trees could be fellow because, normally, it is impossible to provide evolutionary events on branches of trees.

2 A very normal question is that when getting one or some new strains' sequences of a species with published phylogenetic analysis, how to add it to existing or published phylogenetic tree to get the phylogenetic information of new strains?

Due to different sequencing platforms, different sequencing qualities, different sequencing libraries, different assembly softwares with different parameters, different mapping methods with different parameters, different phylogenetic analysis methods/models/softwares/parameters and so many other differences on situations or research methods when doing phylogenetic analysis, it is always impossible and of course very non-convenient to rebuild published trees. Sometimes, it is hard to keep the original structure of a published tree after adding new samples, especially for new strains which have only draft genome or low quality sequences. And sometimes, it is even very difficult to rebuild a published tree using the original samples' sequences. 

It is back to the question 1, normally published trees are hard or maybe impossible to use or compare.

3 Is it possible to build a "tracable" phylogenetic tree database to make it easy to get phylogenetic locations of new samples?

4 How to know what happened on the branches of a tree? Sometimes, information of branches鈥?length are not enough for phylogenetic analysis, detailed evolutionary information may reveal more functional issues.

###==============SaRTree: A pipeline to solve the problems above==============

1 First, the main idea of SaRTree is to get detailed evolutionary events on each branch of a tree, so that comparing with the events, new samples could be located on the existing tree.

2 By a statistics method, recombination events could be detected by RecDetect so that it is possible to get really mutational SNPs for further analysis and, as a result, high resolution phylogenetic reconstruction and divergence time estimation could be possibly implemented.

3 Any third-part software should be use for assembly, alignment, mapping, tree building, divergence date estimating etc. It gives user more freedom to fit different complex situations. SaRTree only provide a platform or pipeline with necessary programs to finish the same thing with different method. Of course, we have recommended software list, but it is hard to say they are good for all situations.

4 Using a simple algorithm in StrainLocater, new samples with variant backgrounds could be located to existing tree with SaRTree output files or same format files. It is very fast and convenient to use an existing tree and then the "tracable tree database" project could start and be possible to use.

------------------------------------------

##Version & License:
------------------------------------------
SaRTree is for locating evolutionary events onto phylogenetic trees or building high resolution trees based on mutation events, StrainLocater is for locating new strains onto existing trees.

Please note that this is a debug test version v1.2, so not all of the functions work well now and many parts of the codes run slowly and redundantly.

There is no warranties coming with SaRTree, so that users must be responsible for the results generated by SaRTree.

This is a free software following GNU General Public License <http://www.gnu.org/licenses/gpl.html>. The license file is also included in this package.

Copyright and Contact: Dalong Hu, University of Sydney (dalong.hu@sydney.edu.au) Jan-2018

------------------------------------------

##Citation:
------------------------------------------
The paper hasn't published, please wait...

------------------------------------------

##Install:
------------------------------------------
SaRTree is wrote by perl, it is no need to install it. Just unzip the package and put the directory to anywhere you want and install all dependencies before using. **After installed all dependencies, please revise the parameters in "Config" file to setup default settings for third-party software.**

###Dependencies Install:

####Conda users:
Please run **'sh install-conda.sh'** to install all the dependencies automatically. Then a environment called **'env_sartree'** will be created and please use **'source activate env_sartree'** to access and run SaRTree.

####Ubuntu users:
Please try **'sudo perl install-ubuntu.pl'** and follow the guide.

####Other users:
Please run the script **"install.pl"** to install perl modules automatically and install raxml, mauve and beast manually and then revise the paths and parameters of them in the configure file **"Config"** manually.

####Manually install:
Ubuntu users:
Please try **'sudo cpan -i FindBin Getopt::Long Bio::SeqIO Bio::TreeIO Bio::AlignIO'** to install perl packages and try **'sudo apt-get install mauve-aligner raxml beast-mcmc'** to install third-party software.

Conda users:
Please add channels "conda-forge", "defaults" ,"r", "bioconda" and "esteinig" by 
	**'conda config --add channels conda-forge r bioconda esteinig'**
	
And then create a running environment with dependencies installing by 
	**'conda create -y --name env_sartree perl perl-yaml-libyaml perl-bioperl perl-findbin perl-getopt-long raxml mauve beast'**
	
Then run 
	**'source activate env_sartree'** 
	
to access SaRTree environment and run SaRTree.

------------------------------------------

##Dependencies:
------------------------------------------

###==============Mandatory==============
1. Operating System: any GNU-Linux system with normal commands (cp, mkdir, cat, cut and paste)
2. perl v5.10.1 or later versions
3. perl Modules:
    Bio
    Bio::SeqIO
    Bio::TreeIO
    Bio::AlignIO
    FindBin
    Getopt::Long
4. Mauve v2.4.0 or later versions
5. RAxML v8.1.17 or later versions

**Please run script "install.pl" to get the perl Modules by cpan while your may need to provide root permission to install some of them (try "sudo perl install.pl" if there is a "permission denied" error occured) and set up the path to dependencies software.**

###==============Optional============== (please revise the paths and parameters of them in the configure file "Config")

1.bwa/bowtie (not supported now, will add them in the next version)
2.BEAST v1.8.4 (not support BEAST 2 now, will fix in the next version)
3.FigTree (any version, for tree viewing)

(All the earlier versions of the software above might be OK. We haven't tested all the early versions.)

------------------------------------------


##Usage:
------------------------------------------

###==============SaRTree core==============

1. make a working directory (e.g. mkdir test)


2. (optional) put reference sequence(fasta file, complete genome, 1 chromosome) and its annotation(gff file) into working directory (use "strainname".fna and "strainname".gff as the format of file names, e.g. AB0057.fna and AB0057.gff) (will support multi chromosomes in formal version)


3. choose a model which depends on the input files("standard" for mauve result, "raw" for raw sequence data, "formatted" for SaRTree format mapping or snp listi, for continuing existing run interupted by rooting the tree)

    for "standard" model, (the recommended running model)
make a directory named "mauve"(without the quotes) in the working directory (temporarily only mauve result rupport) put all pairwise aligment result by Mauve(btw ref and each strain) into the directory (use "strainname".mauve as the format of file name, e.g. ZJ06.mauve)

    for "raw" model, (run mauve in the script, will support multi-threads in formal version) make a directory named "seq"(without the quotes) in the working directory
put all raw samples' sequences (fasta format, use "strainname".fna as the format of file name, e.g. ZJ06.fna) into the directory

    for "formatted" model, (not recommended, just for test) put SaRTree format snp list file (a tab separated file named as "tag".list,the "tag" should be same as the one using for the -t option, with the format:

    location	ref	-s1-	-s2-	...

    1234	A	T	T	...

    2345	T	T	A	...

    ) and the name list file (tab separated txt file, with format:

    ref	AB0057

    -s1-	ZJ06

    -s2-	...

    ...

    ) into working directory

    for "root" model,
use the same commond as the first run and add -u to indicate the rooted tree.

    for "beast" model,
use the same commond as the first run and add -d to indicate the isolation dates and add -u to indicate the nex format output tree from TreeAnnotator.

4. run the main script:

    perl SaRTree.pl -i &lt;work dir> -t &lt;output tag> -f &lt;ref sequence(fna)> -g &lt;ref annotation(gff)> -m &lt;standard/raw/formatted/root/beast> [options]


###==============StrainLocater==============

1. make a working directory (e.g. mkdir test) and make a directory named "query" in the working directory


2. (optional) put reference sequence(fasta file, complete genome, 1 chromosome) into working directory (use "strainname".fna as the format of file name, e.g. AB0057.fna) (will support multi chromosomes in formal version)


3. choose a model which depends on the input files("standard" for mauve result, "raw" for raw sequence data, "formatted" for SaRTree format mapping or snp list)


    the details of models are similar to the "SaRTree core" section above, but NOTE that all the input files should be put into "work dir"query/seq ("raw" model), "work dir"/query/mauve ("standard" model)or "work dir"/query/ ("formatted" model)
Adding a "query" directory separately is to make it possible to run SaRTree core and StrainLocater in the same working directory conveniently


4. run the main script: (indicating newick format target tree file and SaRTree format, see the "formatted" model section, snp list file for the tree)

    perl SaRTree.pl [options] -l -i &lt;work dir> -t &lt;output tag> -f &lt;ref sequence(fna)> -u &lt;target tree(nwk)> -s &lt;snp list(SaRTree format)> -m &lt;standard/raw/formatted>

###==============treeview.pl==============

a small script to show a tree on commond line text interface, for details, please run perl script/treeview.pl -h

###==============Options==============

	-dir/-i <string>  mandatory, **FULL(!)** path to the working directory, that is a problem of RAxML, will revise that in next version

	-tag/-t <string>  mandatory, output tag, just a tag, not to add path here

	-ref/-f <string>  mandatory, full path to reference sequence file (fasta format, named as "strainname".fna, e.g. AB0057.fna)

	-gff/-g <string>  mandatory, full path to annotation file of reference strain (gff format, named as "strainname".gff, e.g. AB0057.gff)

	-model/-m <standard/raw/formatted/root/beast>  mandatory, three input models could be selected depending on the input files, for details see description above

	-locater/-l  optional, to run StrainLocater instead of SaRTree main process, when using -l, a target tree and its snp list should also be specified by -u and -s respectively, and a directory named "query" should be built in working directory then all input/output files for StrainLocater will be in the directory

	-userTree/-u [string]  optional, a newick format tree file, when -l used, -u specifies the target tree to locate. when without -l, -u specifies a user's own phylogenetic tree instead of building a tree in SaRTree process

	-snp/-s [string]  optional, when -l used, -s specifies the snp list of the target tree

	-root/-r [auto/man/off] optional, the method to root the tree, default "auto"; "auto": root the tree automatically by StrainLocater based on events (only this method could locate events on the two descendant branches of root correctly); "man": program stop after tree building, then user could root the tree manually, but then should restart SaRTree again by indicating "root" model as use -u to indicate the rooted tree, and NOTE that the events on two descendant branches of root may not be right, all of the events of them will be counted to 'the branch including reference'; "off": just use the output tree of tree building software, NOT recommended

	-outgroup/-o [string]  optional, outgroup strain's name could be specified by -o; when using "auto" model for -r, outgroup will be used by StrainLocater;  when using "off" or "man" for -r, outgroup will be used by raxml (NOT recommended, should root the tree carefully)

	-beast/-b [off/auto/man] optional, run beast to estimate divergence time, a file with isolattion date should be specified by -d, not recommended now, there are some bugs, default off

	-date/-d [string]  optional, when -b used, -d specifies the isolattion date file with the format(tab separated):
	strain1	1970
	strain2	1980
	...

	-coverage/-c [integer]  optional, an integer btw 0-100, coverage to reference filter, all samples with lower coverage to reference will be put into a "query" folder for second run using strainlocater, default 0

	-filter/-e [off/rec/multi/both/either]	optional, filter out SNPs with low confidence, "off" indicates close this function;"rec" mode deletes SNPs covered by recombination regions ignroing recombinant SNPs' distribution;"multi" mode deletes SNPs presenting more than two types of bases;"both" mode deletes only SNPs fitting both detected by "rec" and "multi" modes;"either" mode deletes SNPs fitting all SNPs fitting either "rec" and "multi" modes. default "off"

	-config/-p [string]  optional, a config file with parameters, please use '-' or leave the line out to indicate a default or unavailble parameter. See the config.template.txt in the SaRTree folder for details. When using the config file, all parameters in the file will cover any other input parameters including the commond line input and any default.

	-thread/-a [integer]  optional, number of threads to use in multi-threads model.

	-help/-h  optional, show this usage

	-version/-v  optional, show version information

------------------------------------------

##Output Interpretation:
------------------------------------------

After finishing all the process of the main script, all the meaningful output files will be put in a directory named "result" which is built automatically in the working directory. The files in "tmp" directory are temporary files which are used to trace back to find bugs, recommend to delete the directory after checking all the result.

The information within the files are as below:
(will rename all the output files in formal version to make them more reasonable)

###==============SaRTree==============

    example/result/<tag>.name.nwk: when using RAxML to build the phylogenetic tree, this is a newick format tree file with the original branch length and bootstrap re sult based on the mutation events
    example/result/<tag>.event.name.nex: the tree with all evolutionary events located on branches in nexus format, "events" section of each branch shows the located events, using the format "__(number of mutations)__(number of recombinations)__(number of indels)__", which are easy to show by a third party software FigTree using its "node bar" or "branch bar" function.
    example/result/<tag>.date.nex: created when -b used, nexus format tree file with divergence dates as branches lengths when BEAST is using, for more detail of it, please see documents of BEAST, and to show the data, a third party software FigTree is recommended.
    example/result/branch.xls: a summary of all evolutionary events and other information of each branch
    example/result/branch.xls.pattern2branch.xls: correspondence between pattern type and branch index
    example/result/<tag>.snp.final.xls: a summary of all mutation events' information including branches located, annotation and brief profiles
    example/result/<tag>.rec.final.xls: a summary of all recombination events' information including branches located, annotation and brief profiles
    example/result/<tag>.indel.final.xls: created when "raw" or "standard" model used, a summary of all indel events' information including lengths of sequences of each sample at the indel regions, branches located, annotation and brief profiles
    example/result/<tag>.indel.seq.name.xls: created when "raw" or "standard" model used, a summary of all indel events' sequences of each sample at the indel region.

###==============StrainLocater==============

    example/query/result/<tag>.nodes.txt: description of each indexed node
    example/query/result/<tag>.ourOfTree.txt: the strains without good enough location on the tree are excluded out of the result
    example/query/result/<tag>.onTree.txt: a summary of the strains could be located on the target tree and their locations
    example/query/result/<tag>.tree.nex: a nexus tree file with information of new strains' locations, "located" section of each branch shows the located strains, which could be showed by third party software FigTree using "node bar" or "branch bar" function.
------------------------------------------

##Example:
------------------------------------------
A set of Acinetobacter baumannii genomes is included as an example in the "example" directory.

**The script "example.sh" includes all the basic commands to run this example, could just try "sh example.sh" to test.**

###==============input files details==============

	example/ACICU.fna	:reference genome sequence file
	example/ACICU.gff	:reference genome annotation file
	example/date.txt	:(for BEAST, USELESS in this version)2-column tab-separated text file with isolation dates information
	example/seq/<name>.fna	:all the other strains' complete genome sequence files for SaRTree core
	example/query/seq/<name>.fna	:all the query strains' draft genome sequence files for StrainLocater

###==============example commands==============

    cd /home/user/SaRTree-v1.2
    perl SaRTree.pl -i /home/user/SaRTree-v1.2/example/ -t test -f example/ACICU.fna -g example/ACICU.gff -m raw -o 307-0294
    perl SaRTree.pl -l -i /home/user/SaRTree-v1.2/example/ -t test -f example/ACICU.fna -u example/result/test.name.nwk -s example/result/test.snp.final.xls -m raw

After all processes finished, all output files for SaRTree core should be found in "example/result" directory and output files for StrainLocater in "example/query/result".

------------------------------------------

###==============example for manually rooting==============

    cd /home/user/SaRTree-v1.2
    perl SaRTree.pl -i /home/user/SaRTree-v1.2/example/ -t test -f example/ACICU.fna -g example/ACICU.gff -m raw -r man

(Then root the output tree "/home/user/SaRTree-v1.2/example/result/test.name.nwk" by hand and save the rooted tree as nwk format for example "/home/user/SaRTree-v1.2/example/test.rooted.nwk" )

    perl SaRTree.pl -i /home/user/SaRTree-v1.2/example/ -t test -f example/ACICU.fna -g example/ACICU.gff -m root -r man -u /home/user/SaRTree-v1.2/example/test.rooted.nwk

Then will go on the following steps. The "beast" model is similar.

------------------------------------------

###==============example config==============

    cd /home/user/SaRTree-v1.2
    perl SaRTree.pl -p example/test.config.txt

With only a config file indicating all usefull parameters, SaRTree could be repeated easily

------------------------------------------

##Other Notices:
------------------------------------------

1. For tree painting, SaRTree is not focusing on this part, only provide a small script named treeview.pl to show phylogenetic tree on commond line text interface. All the useful information within a tree is included in nexus format tree files. User could choose any tree-viewing software or user's own scripts which could open nexus format files to visualize the output data. It is recommended to use a third party software "FigTree" <http://tree.bio.ed.ac.uk/software/figtree/> to show the information.
2. Due to the complexity of the parameters of BEAST, it is recommended to run BEAST independently.
3. It is also recommended to run MAUVE independently by user.
4. All the scripts included in "script" directory could be ran independently for debugging.
5. No temporary file will be deleted automatically for debugging.
------------------------------------------

##Things Coming in Next Version:
------------------------------------------

1. supporting for more recombination-event-detecting software
2. supporting for mapping software and result input
3. supporting for multi chromosomes
4. supporting for parameters of RecDetect
5. supporting for complex algorithms for StrainLocater
6. optimizing all codes and performance
7. supporting for more tree-building software
8. supporting for step-wise running, so that user could run BEAST, RAxML and other third party softwares independently
9. trouble shooting and input/output checking
------------------------------------------
