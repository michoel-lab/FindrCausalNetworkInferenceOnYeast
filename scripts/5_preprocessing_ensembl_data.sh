#!/bin/sh
#  Processing of Ensembl file listing for yeast:
#
echo '# Start processing yeast ensembl data'

inpath=../data/input/
outpath=../data/output/

file1=Saccharomyces_cerevisiae.R64-1-1.83.gff3
out1=cleaned_genes_ensembl83_yeast_step1.txt
out2=cleaned_genes_ensembl83_yeast_step2.txt

sed -n -e '/gene/{p;n;}' $inpath$file1 | sed -e 's/;/ /g' -e 's/description.*/ /' > $outpath$out1
#sed  -e 's/:/ /g' -e "s/[[:space:]]\+/ /g" $outpath$out1 > $outpath$out2
sed -n  '/ID=gene/p' $outpath$out1 | sed -n  '/Name/p'   > $outpath$out2

echo '# Done'
# EOF
