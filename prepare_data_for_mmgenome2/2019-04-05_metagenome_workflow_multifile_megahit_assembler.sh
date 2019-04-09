#########################################################################################
#											#
# Example solution for generating a metagenome 					        #
# assembly and coverage files for binning in the mmgenome2 package                      #
# (https://github.com/KasperSkytte/mmgenome2) 						#
#											#
# What you need to run this                                                             #
# - install the dependencies on your server and update the paths in this script         #
# - supply a file named "samples" with the sequencing IDs for your paired end data      #
#    (no file extension)                                                                #
# - update the path to the place where the script can find your sequencing data         #
# - update the number of threads/CPUs that the software is allowed to use               #
# 											#
#											#
# Best wishes the mmgenome2 team 							#
#											#
#########################################################################################



#############################
####      Settings       ####
#############################

THREADS=60; # Number of threads to use
DATAPATH=/space/users/rkirke08/Desktop/rkirkegaard/MetaGenomes/ZymoMock/data/;

# Adaptors
NEX_ADP1=CTGTCTCTTATACACATCT # Illumina Nextera adaptor sequences
NEX_ADP2=CTGTCTCTTATACACATCT # Illumina Nextera adaptor sequences
TRU_ADP1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA # Illumina TruSeq and NEB Nebnext adaptor sequences
TRU_ADP2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT # Illumina TruSeq and NEB Nebnext adaptor sequences
TRIMLENGTH=80; # Remove reads shorter than X bp (adjust based on your sequencing length)
QUALITY=20; # Remove reads with a phred score below


#############################
####    Dependencies     ####
#############################
# samples (file with sequencing IDs)
# Software | Tested and working with version:
# R (v. 3.4.4)
MEGAHIT=megahit; # Path to megahit (v 1.1.3) 
PROKKA=prokka; # Path to prokka (v. 1.12)
BARRNAP=barrnap; # Path to barrnap
PRODIGAL=prodigal; # Path to prodigal (v. 2.6.2)
HMMSEARCH=hmmsearch; # HMMER (v. 3.0)
MINIMAP2=/space/users/smk/Software/minimap2-2.5/minimap2-2.5_x64-linux/minimap2; # Path to minimap2 (v. 2-2.5)
SAMTOOLS=samtools; # Path to samtools (v. 1.3.1)
USEARCH=usearch10;  # Path to usearch (v. 10)
CUTADAPT=cutadapt; # Path to  cutadapt (v. 1.16)
KAIJU=/space/users/smk/Software/kaiju/bin; # Path to kaiju (v. 1.6.0)
KAIJU_DB=/space/users/smk/Software/kaiju/database; # Path to kaiju database
MMPATH=MMscripts/; # Path to some scripts

#############################
####      Workflow       ####
#############################

date | tee -a log.txt;
echo "Metagenome workflow started" | tee -a log.txt;

# Create directories
mkdir -p temp
mkdir -p results
mkdir -p data

#######################
#### Preprocessing ####
#######################
echo "Finding your samples and copying them to the data folder" | tee -a log.txt;
while read samples
do
NAME=$samples;
find $DATAPATH -name $NAME*R*fastq.gz 2>/dev/null -type f -exec cp {} . \;
(ls $NAME*R1*fastq.gz >> /dev/null 2>&1) || echo "$NAME R1 not found" | tee -a log.txt;
(ls $NAME*R2*fastq.gz >> /dev/null 2>&1) || echo "$NAME R2 not found" | tee -a log.txt;
for file in $NAME*R1*fastq.gz; do cat $file >> data/$NAME.R1.fastq.gz; done
for file in $NAME*R2*fastq.gz; do cat $file >> data/$NAME.R2.fastq.gz; done
done < samples
echo "Number of SEQIDs" | tee -a log.txt;
wc -l samples | tee -a log.txt;
echo "Number of fastq files" | tee -a log.txt;
ls data/ | grep "fastq.gz" | wc -l | tee -a log.txt;


gunzip data/*.fastq.gz;
rm *fastq.gz;

# Remove phiX
date | tee -a log.txt;
echo "PhiX filtering started" | tee -a log.txt;
# # https://www.drive5.com/usearch/manual/cmd_filter_phix.html
# usearch10 -filter_phix reads_fwd.fq -reverse reads_rev.fq -output fil_fwd.fq -output2 fil_rev.fq
while read samples
do
NAME=$samples;
$USEARCH -filter_phix data/$NAME.R1.fastq -reverse data/$NAME.R2.fastq -output data/$NAME.R1.filtered.fastq -output2 data/$NAME.R2.filtered.fastq -threads $THREADS
done < samples

date | tee -a log.txt;
echo "Adaptor trimming started" | tee -a log.txt;
# Trim adaptors with cutadapt in parallel (http://seqanswers.com/forums/showthread.php?t=65982)
find data/*.R1.filtered.fastq | sed 's/.R1.filtered.fastq$//' | parallel "$CUTADAPT -a $NEX_ADP1 -a $TRU_ADP1 -A $NEX_ADP2 -A $TRU_ADP2 -m $TRIMLENGTH -q $QUALITY -o {}.R1.trimmed.fastq --paired-output {}.R2.trimmed.fastq {}.R1.filtered.fastq {}.R2.filtered.fastq"

echo "combine trimmed reads" | tee -a log.txt;
cat data/*.R1.trimmed.fastq > temp/combined.trimmed.R1.fq
cat data/*.R2.trimmed.fastq > temp/combined.trimmed.R2.fq

##################
#### Assembly ####
##################
date | tee -a log.txt;
echo "Metagenome assembly started" | tee -a log.txt;
$MEGAHIT -1 temp/combined.trimmed.R1.fq -2 temp/combined.trimmed.R2.fq --k-list 43,71,99,127 \
-t $THREADS -m 0.7 -o temp/megahit 
# Remove contigs shorter than 1000 bp and clean contig names
awk '!/^>/ { next } { getline seq } length(seq) >= 1000 { print ">" ++i "\n" seq }' \
temp/megahit/final.contigs.fa > results/assembly.fasta

###########################
#### Generate coverage ####
###########################
date | tee -a log.txt;
echo "Coverage generation started" | tee -a log.txt;
while read samples
do
NAME=$samples;
echo $NAME | tee -a log.txt;
echo "Map reads to assembly using minimap2 (https://github.com/lh3/minimap2)" | tee -a log.txt;
$MINIMAP2 -ax sr -t $THREADS results/assembly.fasta data/$NAME.R1.trimmed.fastq data/$NAME.R2.trimmed.fastq | $SAMTOOLS view --threads $THREADS -Sb - | $SAMTOOLS sort --threads $THREADS - -o temp/$NAME.aln.sorted.bam
$SAMTOOLS depth temp/$NAME.aln.sorted.bam > temp/$NAME.depth.txt
FILEend="_cov.csv";
Rscript $MMPATH/calc_coverage_script_datatable.r temp/$NAME.depth.txt results/$NAME$FILEend $THREADS
done < samples

$SAMTOOLS merge --threads $THREADS - temp/*.aln.sorted.bam | $SAMTOOLS view --threads $THREADS -h -o temp/combinedmapping.bam -

#########################
#### Additional info ####
#########################
date | tee -a log.txt;
echo "Metagenome annotation started" | tee -a log.txt;
$PROKKA --outdir temp/prokka --cpus $THREADS --metagenome --prefix prokka results/assembly.fasta

#################
#### get 16S ####
#################
date | tee -a log.txt;
echo "Extracting rRNA sequences" | tee -a log.txt;
$BARRNAP --outseq results/rRNA.fa < results/assembly.fasta

#############################
#### predicting proteins ####
#############################
date | tee -a log.txt;
echo "Finding essential genes - Predicting proteins (Prodigal)" | tee -a log.txt;
$PRODIGAL -a temp/temp.orfs.faa -i results/assembly.fasta -m -o temp/temp.txt -p meta -q
cut -f1 -d " " temp/temp.orfs.faa > temp/assembly.orfs.faa

date | tee -a log.txt;
echo "Finding essential genes - running HMM search" | tee -a log.txt;
$HMMSEARCH --tblout temp/orfs_hmm.txt --cut_tc --notextw \
  --cpu $THREADS $MMPATH/essential.hmm temp/assembly.orfs.faa > /dev/null

echo "scaffold,hmm.id" > results/essential.csv
tail -n+4  temp/orfs_hmm.txt |\
  sed 's/ * / /g' | cut -f1,4 -d " " |\
  sed -e 's/_/ /' -e '/^#/ d' | tr " " "," |\
  cut -f1,3 -d"," >> results/essential.csv

date | tee -a log.txt;
echo "Kaiju protein taxonomic classification" | tee -a log.txt;
printf "\nTaxonomic classification of contigs using Kaiju\n\n"
$KAIJU/kaiju -p -z $THREADS -t $KAIJU_DB/nodes.dmp -f $KAIJU_DB/kaiju_db.fmi \
  -i temp/assembly.orfs.faa -o temp/kaiju.out
  $KAIJU/addTaxonNames -u -r phylum -t $KAIJU_DB/nodes.dmp -n $KAIJU_DB/names.dmp \
  -i temp/kaiju.out -o temp/kaiju_names.out

date | tee -a log.txt;
echo "Majority vote contig classification" | tee -a log.txt;
echo "scaffold,phylum" > results/tax.csv
cat temp/kaiju_names.out | \
  sed -e 's/_/\t/' -e '/NA;/d' -e 's/; //'  | \
  cut -f2,5 | \
  awk -F "\t" '{a[$1","$2]++} END{OFS = ","; for (i in a) print i, a[i]}' - | \
  awk -F "," '{if (c[$1]<$3){a[$1]=$1; b[$1]=$2; c[$1]=$3}; d[$1]+=$3} \
  END{OFS = ","; for (i in a){if (c[i] >= 2 && c[i]/d[i] > 0.51) print a[i], b[i] }}' - |\
  sort -n -t, -k1,1 >> results/tax.csv

date | tee -a log.txt;
echo "Generating connection network" | tee -a log.txt;
$SAMTOOLS view -h --threads $THREADS -o temp/combinedmapping.sam temp/combinedmapping.bam
perl $MMPATH/network.pl -infile temp/combinedmapping.sam -minconnections 2 -outcon results/network.txt


date | tee -a log.txt;
echo "Metagenome workflow finished" | tee -a log.txt;
echo "Program finished (Check that no errors or warnings were produced)" | tee -a log.txt;
