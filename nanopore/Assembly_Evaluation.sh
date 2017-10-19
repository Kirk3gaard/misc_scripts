#!/bin/sh

# Genome assembly evaluation with E. Coli K12 MG1655 (Ref U00096.2, 4639675 bp)
# We have illumina data from SRA (SRR2627175, PE 270 Mbp)
# and 1D R9.4 nanopore data (183 Mbp)
# Miniasm assembly 
# + racon correction 
# + pilon polishing
# CANU assembly
# + nanopolish
# + pilon polishing
# Unicycler assembly

# data
NANOPORE=/space/users/rkirke08/Desktop/pikachutest/basecalling/EColiK12MG1655/basecalled1D_EcoliK12MG1655_alba2_0_1/workspace/pass/NANO_subset.fastq;
NANOPOREFA=NANO_subset.fa;
ILLUMINA_R1=/space/users/rkirke08/Desktop/pikachutest/basecalling/EColiK12MG1655/PE_data/SUBSET_R1.fastq;
ILLUMINA_R2=/space/users/rkirke08/Desktop/pikachutest/basecalling/EColiK12MG1655/PE_data/SUBSET_R2.fastq;
ILLUMINA=/space/users/rkirke08/Desktop/pikachutest/basecalling/EColiK12MG1655/PE_data/SRR2627175_subsampled.fastq;

# Settings
THREADS=24;

# PATHS
MINIMAPPATH=/space/users/rkirke08/Desktop/Miniasm/minimap/;
MINIMAP2PATH=/space/users/rkirke08/Desktop/software/minimap2/;
MINIASMPATH=/space/users/rkirke08/Desktop/Miniasm/miniasm/;
RACONPATH=/space/users/rkirke08/Desktop/software/racon/racon/racon/bin/;
UNICYCLERPATH=/space/users/rkirke08/Desktop/software/Unicycler/;
CANUPATH=/opt/canu-1.6/Linux-amd64/bin/;
NANOPOLISHPATH=/opt/nanopolish/;

###########################################
# Assembly + polishing with various tools #
###########################################

# Miniasm assembly (v. 0.2-r137-dirty)
$MINIMAPPATH/minimap -Sw5 -L100 -m0 -t $THREADS $NANOPORE $NANOPORE | gzip -1 > reads.paf.gz
$MINIASMPATH/miniasm -f $NANOPORE reads.paf.gz > miniasm_assembly.gfa
awk '/^S/{print ">"$2"\n"$3}' miniasm_assembly.gfa > miniasm_assembly.fa

# Racon correction (v. )
### Assembly consensus correction
$MINIMAP2PATH/minimap2 -t $THREADS -x map-ont miniasm_assembly.fa $NANOPORE > mappings1.paf  
$RACONPATH/racon -t $THREADS $NANOPORE mappings1.paf miniasm_assembly.fa consensus1.fasta  

## Multiple rounds of correction gives better accuracy ##
$MINIMAP2PATH/minimap2 -t $THREADS -x map-ont consensus1.fasta $NANOPORE > mappings2.paf   
RACONPATH/racon -t $THREADS $NANOPORE mappings2.paf consensus1.fasta consensus2.fasta 
awk '/^>/{print ">" ++i; next}{print}'  consensus2.fasta > racon2x_miniasm_assembly.fasta

# Pilon polish Racon2x corrected assembly (https://github.com/broadinstitute/pilon/releases/)
# Create index files
bwa index racon2x_miniasm_assembly.fasta
# Map reads using BWA
bwa mem -p racon2x_miniasm_assembly.fasta $ILLUMINA_R1 $ILLUMINA_R2 > aln.sam
# Convert SAM file to BAM
samtools view -Sb -o out.bam aln.sam
samtools sort -m 1000000000 out.bam sorted
samtools index sorted.bam
# Pilon polishing with illumina
java -jar pilon.jar --genome racon2x_miniasm_assembly.fasta --unpaired sorted.bam --output pilon_polished --outdir pilon_polished --threads $THREADS


# Unicycler assembly (v. 0.4.1)
$UNICYCLERPATH/unicycler-runner.py -s $ILLUMINA -l $NANOPORE --threads $THREADS --spades_path /space/users/rkirke08/Desktop/spades/SPAdes-3.10.1-Linux/bin/spades.py --no_correct --no_pilon --racon_path /space/users/rkirke08/Desktop/software/racon/racon/racon/bin/racon -o unicycler_dir


# CANU assembly (v. 1.6)
$CANUPATH/canu -p ecoliK12MG1655 -d CANU-K12MG1655 genomeSize=4.8m -nanopore-raw $NANOPORE

# CANU+Nanopolish (v. 0.8.3)
$NANOPOLISHPATH/nanopolish index -d /space/users/rkirke08/Desktop/rkirkegaard/MinION/datadump/20170405_RHK_MinIONrun27_EcoliK12MG1655_rad002/ $NANOPOREFA
# Index the draft genome
bwa index $ASSEMBLY
# Align the basecalled reads to the draft sequence
bwa mem -x ont2d -t $THREADS $ASSEMBLY $NANOPOREFA | samtools sort -o reads.sorted.bam -T reads.tmp -
samtools index reads.sorted.bam
# -P and -t needs to be adjusted so that P*t does not exceed the number of available threads
python $NANOPOLISHPATH/scripts/nanopolish_makerange.py $ASSEMBLY | parallel --results nanopolish.results -P 6 \
    $NANOPOLISHPATH/nanopolish variants --consensus polished.{1}.fa -w {1} -r $NANOPOREFA -b reads.sorted.bam -g $ASSEMBLY --threads 4 --min-candidate-frequency 0.1
python $NANOPOLISHPATH/scripts/nanopolish_merge.py polished.*.fa > polished_genome.fa



# Pilon polish nanopolished CANU assembly
# Create index files
bwa index CANU_nanopolished_genome.fa

# Map reads using BWA
bwa mem -p CANU_nanopolished_genome.fa $ILLUMINA_R1 $ILLUMINA_R2 > CANU_Nanopolished_aln.sam

# Convert SAM file to BAM
samtools view -Sb -o CANU_Nanopolished_out.bam CANU_Nanopolished_aln.sam
samtools sort -m 1000000000 CANU_Nanopolished_out.bam -o CANU_Nanopolished_sorted.bam
samtools index CANU_Nanopolished_sorted.bam


# Pilon polishing with illumina 
java -jar pilon.jar --genome CANU_nanopolished_genome.fa --unpaired CANU_Nanopolished_sorted.bam --output pilon_polished_CANU_nanopolished --outdir pilon_polished_CANU_nanopolished --threads $THREADS

##############################
# Comparison                #
##############################
#######
# ANI #(Perl script from: https://github.com/chjp/ANI)
#######

# Miniasm
perl ANI.pl -bl /space/sharedbin/blast-2.2.25/bin/blastall -fd /space/sharedbin/blast-2.2.25/bin/formatdb -qr miniasm_assembly.fa -sb EcoliK12MG1655_reference.fasta -od ANI/result_raw
# Racon1x corrected
perl ANI.pl -bl /space/sharedbin/blast-2.2.25/bin/blastall -fd /space/sharedbin/blast-2.2.25/bin/formatdb -qr consensus1.fasta -sb EcoliK12MG1655_reference.fasta -od ANI/result_racon1x
# Racon2x corrected
perl ANI.pl -bl /space/sharedbin/blast-2.2.25/bin/blastall -fd /space/sharedbin/blast-2.2.25/bin/formatdb -qr racon2x_miniasm_assembly.fasta -sb EcoliK12MG1655_reference.fasta -od ANI/result_racon2x
# Pilon polished
perl ANI.pl -bl /space/sharedbin/blast-2.2.25/bin/blastall -fd /space/sharedbin/blast-2.2.25/bin/formatdb -qr pilon_polished.fasta -sb EcoliK12MG1655_reference.fasta -od ANI/result_pilon
# Unicycler
perl ANI.pl -bl /space/sharedbin/blast-2.2.25/bin/blastall -fd /space/sharedbin/blast-2.2.25/bin/formatdb -qr unicycler_assembly.fasta -sb EcoliK12MG1655_reference.fasta -od ANI/result_unicycler
# CANU
perl ANI.pl -bl /space/sharedbin/blast-2.2.25/bin/blastall -fd /space/sharedbin/blast-2.2.25/bin/formatdb -qr ecoliK12MG1655.contigs.fasta -sb EcoliK12MG1655_reference.fasta -od ANI/result_canu
# CANU+Nanopolish
perl ANI.pl -bl /space/sharedbin/blast-2.2.25/bin/blastall -fd /space/sharedbin/blast-2.2.25/bin/formatdb -qr CANU_nanopolished_genome.fa -sb EcoliK12MG1655_reference.fasta -od ANI/result_canu_nanopolish
# CANU+Nanopolish+Pilon
perl ANI.pl -bl /space/sharedbin/blast-2.2.25/bin/blastall -fd /space/sharedbin/blast-2.2.25/bin/formatdb -qr pilon_polished_CANU_nanopolished.fasta -sb EcoliK12MG1655_reference.fasta -od ANI/result_canu_nanopolish_pilon
# spades
perl ANI.pl -bl /space/sharedbin/blast-2.2.25/bin/blastall -fd /space/sharedbin/blast-2.2.25/bin/formatdb -qr best_spades_graph.fasta -sb EcoliK12MG1655_reference.fasta -od ANI/result_spades
# Reference
perl ANI.pl -bl /space/sharedbin/blast-2.2.25/bin/blastall -fd /space/sharedbin/blast-2.2.25/bin/formatdb -qr EcoliK12MG1655_reference.fasta -sb EcoliK12MG1655_reference.fasta -od ANI/result_reference

##############
# Annotation #
##############
# Miniasm assembly
prokka miniasm_assembly.fa --outdir miniasm_anno
# Racon1x corrected
prokka consensus1.fasta --outdir racon1x_anno
# Racon2x corrected
prokka racon2x_miniasm_assembly.fasta --outdir racon2x_anno
# Racon+Pilon corrected
prokka pilon_polished.fasta --outdir pilon_anno
# CANU
prokka CANU_ecoliK12MG1655.contigs.fasta --outdir CANU_anno
# CANU+nanopolish
prokka CANU_nanopolished_genome.fa --outdir CANU_Nanopolish_anno
# CANU+nanopolish+pilon
prokka pilon_polished_CANU_nanopolished.fasta --outdir CANU_Nanopolish_pilon_anno
# Unicycler
prokka unicycler_assembly.fasta --outdir unicycler_anno
# Unicycler
prokka best_spades_graph.fasta --outdir spades_anno
# Reference
prokka EcoliK12MG1655_reference.fasta --outdir reference_anno

### Convert gff to table
# Spades
perl gff.prokka.to.table.pl -i spades_anno/PROKKA_10172017.gff -o tables/spades_annotation.txt
# Miniasm
perl gff.prokka.to.table.pl -i miniasm_anno/PROKKA_10132017.gff -o tables/miniasm_annotation.txt
# Racon1x
perl gff.prokka.to.table.pl -i racon1x_anno/PROKKA_10172017.gff -o tables/racon1x_annotation.txt
# Racon2x
perl gff.prokka.to.table.pl -i racon2x_anno/PROKKA_10132017.gff -o tables/racon2x_annotation.txt
# Racon2x+pilon
perl gff.prokka.to.table.pl -i pilon_anno/PROKKA_10132017.gff -o tables/racon2xpilon_annotation.txt
# CANU
perl gff.prokka.to.table.pl -i CANU_anno/PROKKA_10152017.gff -o tables/canu_annotation.txt
# CANU+nanopolish
perl gff.prokka.to.table.pl -i CANU_Nanopolish_anno/PROKKA_10162017.gff -o tables/canunanopolish_annotation.txt
# CANU+nanopolish+pilon
perl gff.prokka.to.table.pl -i CANU_Nanopolish_pilon_anno/PROKKA_10172017.gff -o tables/canunanopolishpilon_annotation.txt
# UNIcycler
perl gff.prokka.to.table.pl -i unicycler_anno/PROKKA_10132017.gff -o tables/unicycler_annotation.txt
# Reference
perl gff.prokka.to.table.pl -i reference_anno/PROKKA_10132017.gff -o tables/reference_annotation.txt

##############
# Run Checkm # (CheckM v1.0.7)
##############
checkm lineage_wf ./bins ./output --threads $THREADS --extension fasta --file checkm_report.txt
