bioawk -c fastx '{OFS=""} /chrM/ {print ">"$name,"\n",$seq}' ~/ref/genomes/sacCer3/sacCer3.fa > sacCer3.chrM.fa
tRNAscan-SE -a sacCer3.chrM.trnas.fa sacCer3.chrM.fa
python convert-mt.py sacCer3.chrM.trnas.fa > sacCer3.chrM.trnas.adapted.fa
