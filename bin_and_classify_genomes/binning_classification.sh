#!/bin/bash
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>log.out 2>&1

###############
# Description #
###############
# The purpose of this script is to provide:
# automated binning
# automated classification
# automated annotation
# data ready for mmgenome2
# 
# Given:
# an input assembly
# a directory with trimmed illumina reads for binning

# Settings
INPUTASSEMBLY="/shared-nfs/RHK/Projects/2020_mmlong/hicanu/results/assembly.fasta"
ILMDIR=/shared-nfs/RHK/Projects/2020_mmlong/data/trimmed_data/;
PBDIR=/shared-nfs/RHK/Projects/2020_mmlong/data/
NPDIR=/shared-nfs/RHK/Projects/2020_mmlong/data/
THREADS=100;
TMPDIR=/scratch/tmp_users/RHK

# Modules
MODULE_SAMTOOLS=SAMtools/1.10-foss-2018a;
MODULE_MINIMAP2=Minimap2/2.15-foss-2018a;
MODULE_RACON=Racon/1.3.3-claaudia-amd-foss-2018a;
MODULE_METABAT=MetaBAT/2.12.1-foss-2018a;
MODULE_PROKKA=prokka/1.14.0-foss-2018a-BioPerl-1.7.2;
MODULE_CHECKM=CheckM/1.1.2-foss-2018a-Python-3.6.4;
MODULE_GTDB=GTDBTk/1.0.2-foss-2018a-Python-3.6.4;
MODULE_KAIJU=Kaiju/1.7.0-foss-2018a
MODULE_HMM=HMMER/3.2.1-foss-2018a
MODULE_JAVA=Java/13.0.1
ESSENTIAL=/space/users/smk/Desktop/2017_long_read/pipeline/mmlong/mmlong_tools/databases/essential.hmm;
KAIJU_DB=/shared-nfs/RHK/databases/kaiju/

# Environment
mkdir -p temp/mapping/
mkdir -p results/
module purge

################
# Run workflow #
################

REF=$INPUTASSEMBLY
cp $REF results/

####################
# Map pacbio reads #
####################
date >> log.txt
echo "Map pacbio reads" >> log.txt
module load $MODULE_MINIMAP2
module load $MODULE_SAMTOOLS
while read PBsamples
do
PBSAMPLE=$PBsamples
OUTPUTFILE=temp/mapping/$PBSAMPLE.cov.bam
if [ -s $OUTPUTFILE ]; then echo "$OUTPUTFILE has already been generated"; 
else
minimap2 -ax asm20 --secondary=no -t $THREADS $REF $PBDIR/$PBSAMPLE.fastq |\
  samtools view --threads $THREADS -Sb -F 0x104 - |\
  samtools sort --threads $THREADS - > $OUTPUTFILE
module purge
fi
done < PBsamples
if [ -s $OUTPUTFILE ]; then echo "Successfully generated $OUTPUTFILE" >> log.txt; else echo "Failed generating $OUTPUTFILE" >> log.txt; exit; fi

####################
# Map nanopore reads #
####################
date >> log.txt
echo "Map nanopore reads" >> log.txt
module load $MODULE_MINIMAP2
module load $MODULE_SAMTOOLS
while read NPsamples
do
NPSAMPLE=$NPsamples
OUTPUTFILE=temp/mapping/$NPSAMPLE.cov.bam
if [ -s $OUTPUTFILE ]; then echo "$OUTPUTFILE has already been generated"; 
else
minimap2 -ax map-ont --secondary=no -t $THREADS $REF $NPDIR/$NPSAMPLE.fastq |\
  samtools view --threads $THREADS -Sb -F 0x104 - |\
  samtools sort --threads $THREADS - > $OUTPUTFILE
module purge
fi
done < NPsamples
if [ -s $OUTPUTFILE ]; then echo "Successfully generated $OUTPUTFILE" >> log.txt; else echo "Failed generating $OUTPUTFILE" >> log.txt; exit; fi

####################
# Map illumina reads #
####################
date >> log.txt
echo "Map illumina reads" >> log.txt
module load $MODULE_MINIMAP2
module load $MODULE_SAMTOOLS
while read ILMsamples
do
ILMSAMPLE=$ILMsamples
OUTPUTFILE=temp/mapping/$ILMSAMPLE.cov.bam
if [ -s $OUTPUTFILE ]; then echo "$OUTPUTFILE has already been generated";  
else
minimap2 -ax sr -t $THREADS $REF $ILMDIR/$ILMSAMPLE"_ilmtrim.fq" |\
  samtools view --threads $THREADS -Sb -F 0x104 - |\
  samtools sort --threads $THREADS - > $OUTPUTFILE
fi
done < ILMsamples
module purge
if [ -s $OUTPUTFILE ]; then echo "Successfully generated $OUTPUTFILE" >> log.txt; else echo "Failed generating $OUTPUTFILE" >> log.txt; exit; fi

##############################
# Bin genomes using metabat2 #
##############################
OUTPUTFILE=results/readcov.txt;
if [ -s $OUTPUTFILE ]; then echo "$OUTPUTFILE has already been generated"; 
else
mkdir -p temp/metabat2/bins/
module load $MODULE_METABAT
jgi_summarize_bam_contig_depths --percentIdentity 90 --outputDepth results/readcov.txt temp/mapping/*cov.bam
metabat2 -i $REF -a $OUTPUTFILE -t $THREADS -o temp/metabat2/bins/
module purge
fi
if [ -s $OUTPUTFILE ]; then echo "Successfully generated $OUTPUTFILE" >> log.txt; else echo "Failed generating $OUTPUTFILE" >> log.txt; exit; fi

###################################
# Classify genome bins using GTDB #
###################################
date >> log.txt
echo "Classify dereplicated bins using GTDB" >> log.txt
if [ "$(ls -A temp/gtdb/)" ]; then echo "gtdb has already been generated" >> log.txt;  
else
module load $MODULE_GTDB
mkdir -p /scratch/tmp_users/RHK
export TMPDIR=/scratch/tmp_users/RHK
gtdbtk classify_wf --cpus 20 --genome_dir temp/metabat2/bins/ --out_dir temp/gtdb --extension .fa
cp temp/gtdb/classify/*.summary.tsv results/
module purge
fi

#########################################
# Check quality of genomes using checkm #
#########################################
date >> log.txt
echo "Check genome completeness using checkm" >> log.txt
OUTPUTFILE=results/checkm.tsv
if [ -s $OUTPUTFILE ]; then echo "$OUTPUTFILE has already been generated";  
else
module load $MODULE_CHECKM
mkdir -p TMPDIR
checkm lineage_wf -x fa -t $THREADS --tmpdir TMPDIR --pplacer_threads 10 --reduced_tree --tab_table temp/metabat2/bins/ temp/checkm_results -f results/checkm.tsv
module purge
fi
if [ -s $OUTPUTFILE ]; then echo "Successfully generated $OUTPUTFILE" >> log.txt; else echo "Failed generating $OUTPUTFILE" >> log.txt; exit; fi

####################################
# Annotate metagenome using prokka #
####################################
module load $MODULE_PROKKA
module load $MODULE_JAVA
prokka --cpus $THREADS --metagenome --outdir temp/prokka $REF
module purge

################################
# Annotate contigs using kaiju #
################################
OUTPUTFILE=results/tax.csv
if [ -s $OUTPUTFILE ]; then echo "$OUTPUTFILE has already been generated";  
else
module load $MODULE_HMM
module load $MODULE_KAIJU
# Essential genes
date >> log.txt
echo "Detection of essential genes using HMM" >> log.txt
cat temp/prokka/PROKKA*.faa > temp/prokka/orfs.faa
hmmsearch --tblout temp/hmm.orfs.txt --cut_tc --notextw \
--cpu $THREADS $ESSENTIAL temp/prokka/orfs.faa > /dev/null
echo "scaffold,hmm.id" > temp/essential.csv
tail -n+4  temp/hmm.orfs.txt | sed 's/ * / /g' | \
cut -f1,4 -d " " | sed -e 's/_/ /' -e '/^#/ d' | tr " " "," | cut -f1,3 -d"," \
>> temp/essential.csv

# Kaiju protein taxonomic classification
date >> log.txt
echo "Taxonomic classification of contigs using Kaiju" >> log.txt
kaiju -p -z $THREADS -t $KAIJU_DB/nodes.dmp -f $KAIJU_DB/kaiju_db_refseq.fmi \
-i temp/orfs.faa -o temp/kaiju.out
addTaxonNames -u -r phylum -t $KAIJU_DB/nodes.dmp -n $KAIJU_DB/names.dmp \
-i temp/kaiju.out -o temp/kaiju.names.out

# Majority vote contig classification
echo "scaffold,phylum" > temp/tax.csv
cat temp/kaiju.names.out | \
sed -e 's/_/\t/' -e '/NA;/d' -e 's/; //'  | \
cut -f2,5 | \
awk -F "\t" '{a[$1","$2]++} END{OFS = ","; for (i in a) print i, a[i]}' - | \
awk -F "," '{if (c[$1]<$3){a[$1]=$1; b[$1]=$2; c[$1]=$3}; d[$1]+=$3} \
END{OFS = ","; for (i in a){if (c[i] >= 2 && c[i]/d[i] > 0.51) print a[i], b[i] }}' - |\
sort -n -t, -k1,1 >> temp/tax.csv
cp temp/tax.csv results/
cp temp/essential.csv results/
module purge
fi
if [ -s $OUTPUTFILE ]; then echo "Successfully generated $OUTPUTFILE" >> log.txt; else echo "Failed generating $OUTPUTFILE" >> log.txt; exit; fi

##########################
# Extract 16S rRNA genes #
##########################
OUTPUTFILE=results/16S.fa
if [ -s $OUTPUTFILE ]; then echo "$OUTPUTFILE has already been generated";  
else
module load $MODULE_BARRNAP
date >> log.txt
echo "Detect 16S rRNA genes with Barrnap" >> log.txt
barrnap $REF --threads $THREADS --outseq temp/rRNA.fa
grep -A1 ">16S" temp/rRNA.fa > $OUTPUTFILE
module purge
fi
if [ -s $OUTPUTFILE ]; then echo "Successfully generated $OUTPUTFILE" >> log.txt; else echo "Failed generating $OUTPUTFILE" >> log.txt; exit; fi


date >> log.txt
echo "End of workflow" >> log.txt
