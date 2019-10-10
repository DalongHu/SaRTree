#!/bin/bash
pwd_path=`pwd`
perl SaRTree.pl -c Config-example -i $pwd_path/example/ -t test -f example/ACICU.fna -g example/ACICU.gff -m raw -o 307-0294
perl SaRTree.pl -c Config-example -l -i $pwd_path/example/ -t test -f example/ACICU.fna -u example/result/test.name.nwk -s example/result/test.snp.final.xls -m raw
