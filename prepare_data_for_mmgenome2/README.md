# How to generate input files for mmgenome2?
Following the release of our [mmgenome2 package](https://github.com/KasperSkytte/mmgenome2) it has been clear that it can be tricky to generate some of the files needed for the binning process. 
This repository contains code that will:
  
* generate an assembly from illumina data using the [megahit](https://github.com/voutcn/megahit) assembler, 
* map reads to the assembly using [minimap2](https://github.com/lh3/minimap2) to generate coverage, 
* predict open reading frames using [prodigal](https://github.com/hyattpd/Prodigal), 
* find essential genes using [hmmer](http://hmmer.org/), 
* classify contigs taxonomically using [kaiju](https://github.com/bioinformatics-centre/kaiju).

It is only meant as an example of how the files can be generated and it is important to highlight that it can be modified in many ways. 
The assembler can be exchanged with your preferred solution ([spades](http://cab.spbu.ru/software/spades/), [wtdbg2](https://github.com/ruanjue/wtdbg2), [miniasm](https://github.com/lh3/miniasm/), [canu](https://canu.readthedocs.io/en/stable/index.html)), the kind of reads mapped to generate coverage can be any type you have (illumina, pacbio, nanopore etc.). 
Furthermore, any additional information that might make your binning easier can be added e.g. automatic binning ([metabat2](https://bitbucket.org/berkeleylab/metabat), [autometa](https://bitbucket.org/jason_c_kwan/autometa), etc), whether contigs are circular or linear, whether the contigs are plasmids or not ([PlasFlow](https://github.com/smaegol/PlasFlow)) etc.
