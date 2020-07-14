#20200302

mkdir db && cd db

perl integrall-utils  download  220  integrall | bash
perl integrall-utils  parse     220  integrall >>integrall.txt
tar czvf integrall.tar.gz  integrall && rm  integrall

cut -f1 integrall.txt  | perl -ane  ' print qq{efetch -db nuccore -format gb  -id $F[0]     \>\> integrall.gb\n}'    | parallel -j 3
cut -f1 integrall.txt  | perl -ane  ' print qq{efetch -db nuccore -format fasta  -id $F[0]  \>\> integrall.fasta\n}' | parallel -j 3

cat integrall.fasta | perl -ne 's/>(\S+)\.\d+/>$1/; print' |seqtk seq -C > ../integrall.fasta

pigz -p 40  integrall.fasta

zcat integrall.gb.gz | gbff-utils taxonomy - | taxon-utils  translate  -c2  /biostack/database/taxonomy/taxon.map - > db.taxonomy.txt
fastx-utils view -cl integrall.fasta >db.annotation.txt

bp_genbank2gff3.pl  --quiet integrall.gb   --outdir  ../convert
pigz -p 40  integrall.gb

cd ..
grep "GenBank" convert/integrall.gb.gff  >integrall.gff
rm -rf  convert

#mobile_element
grep "mobile_element" integrall.gff | perl -ane 'print qq{$F[0]\t} . ($F[3] - 1) .qq{\t$F[4]\t$F[0]:$F[3]-$F[4]\t0\t$F[6]\n}' > mobile_element.bed
grep "mobile_element" integrall.gff | perl -ane '$t="-"; $t=$1 if(/mobile_element_type=(.+)/); print qq{$F[0]:$F[3]-$F[4]\t$t\n}' > mobile_element.type
bedtools  getfasta  -fi  integrall.fasta  -bed  mobile_element.bed  -name  -s  -fo  mobile_element.fasta
perl -F: -ane 'print qq{$F[0]\t$_}' mobile_element.type | tsv-utils annotation  taxonomy.txt - | cut -f2,3,4 |  perl -ane 'chomp; @F=split /\t/; print qq{$F[0]\t$F[1] \[$F[2]\]\n}' >mobile_element.annotation.txt

#CDS
grep -P "\tCDS\t" integrall.gff | grep  "protein_id=" | perl -ane  'my ($id)=$_=~/protein_id=(.+?)\./; print qq{$F[0]\t} . ($F[3] - 1) .qq{\t$F[4]\t$id\t0\t$F[6]\n}' >integrall.CDS.bed
bedtools  getfasta  -fi  integrall.fasta  -bed  integrall.CDS.bed  -name  -s  -fo  integrall.CDS.fasta
zgrep -P "\tCDS\t" integrall.gff | grep  "protein_id="  | cut -f1,9  | perl -ane ' my ($id,$name,$product)=("","",""); $id=$1 if($_=~/protein_id=(.+?)\./); $name=$1 if($_=~/Name=(.+?);/); $product=$1 if($_=~/product=(.+?);/); print qq{$id\t$name $product\t$F[0]\n}' | sed 's/"" //'  | tsv-utils annotation  -c 3  taxonomy.txt - | perl -ane 'chomp; @F=split /\t/; print qq{$F[0]\t$F[1],$F[2],\[$F[3]\]\n}' >integrall.CDS.annotation.txt

#mobile_element && CDS
bedtools  intersect  -wao  -a  mobile_element.bed  -b  integrall.CDS.bed   > mobile_element.CDS.bed
cut -f10 mobile_element.CDS.bed | fastx-utils subseq  integrall.CDS.fasta  -  >mobile_element.CDS.fasta

#keywords
grep -P "\tCDS\t" integrall.gff | grep  "protein_id=" > CDS.gff
grep "integrase" CDS.gff > db.gff
grep "integron"  CDS.gff  >> db.gff
grep "cassette"  CDS.gff  >> db.gff
cat  db.gff | cut -f1,9 | perl -ane  'my ($id)=$_=~/protein_id=(.+?)\./;print qq{$F[0]\t$id\n}' | sort | uniq >integron.CDS.txt
rm CDS.gff
cut -f2  integron.CDS.txt | fastx-utils subseq  integrall.CDS.fasta.gz - >integron.CDS.fasta

#retain
cat  mobile_element.CDS.bed integron.CDS.txt  | cut -f1|sort | uniq | tsv-utils subset -r  db.annotation.txt - | grep -P  "IntI|integron|cassette"  > retain.txt
cut -f1 retain.txt | tsv-utils subset integrall.CDS.bed - > retain.bed
cut  -f4   retain.bed | fastx-utils subseq  integrall.CDS.fasta.gz  - >retain.fasta

#db
cat retain.fasta.gz  mobile_element.CDS.fasta.gz   integron.CDS.fasta.gz | fastx-utils  dedup  - >integron.ffn
fastx-utils view  integrall.ffn | tsv-utils subset  curation/integrall.CDS.annotation.txt.gz - >misc/annotation.txt

makeblastdb  -dbtype nucl -title integrall -hash_index  -out  integrall  -in  integrall.ffn
