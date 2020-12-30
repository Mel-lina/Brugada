# Fine-mapping clustered GWAS hits enhances the identification of disease risk and protective genetic variants 
This repository contains all the scripts (Perl, Shell and R) used for the current project on Brugada syndrome cases.


## Pipeline
Here, we describe all the steps followed to obtain the results from this project.
It is important to acknowledge that each script starts with a description of the binaries and other files needed. Also, please do not move the location of the scripts.


### Short read sequencing
* First, execute ```Short_read_alignment.pl```. This script is intended to preprocess FASTQ files from all 86 Brugada syndrome samples and align them to the human reference genome GRCh37/hg19. 
* Second, execute ```GATK_variant_Calling.pl```. This script performs a variant calling using GATK v3.8-0 and also applies the GATK variant quality score recallibration (VQSR).  

### Ancestry admixture analysis
NOTE: For ancestry admixture analysis of Brugada syndrome and Wellderly subjects, we used a different 1000 Genomes reference panel than for GTEx subjects.
* First, execute ```prepare_1KG.sh```. This script pre-processes 1000 Genomes VCF files needed for ancestry admixture analysis of Brugada syndrome and Wellderly subjects.
* Second, execute ```ancestry_BrS.sh``` and ```ancestry_Wdy.sh``` for ancestry admixture analysis of Brugada syndrome and Wellderly subjects, respectively, using 1000Genomes Phase3 as reference panel. 
* Third, execute ```tSNE_BrS.R``` and ```tSNE_Wdy.R``` to run t-SNE on Brugada syndrome and Wellderly subjects and plot ancestry admixture results.
* Finally, execute ```ancestry_GTEx.sh``` for ancestry admixture analysis of GTEx subjects. 

### Long read sequencing
* First, execute ```MinION_alignment_and_filtering.pl``` for alignment and filtering of long-read sequencing data (generated on the MinION device).
* Second, execute ```whatshap_phase.pl``` for haplotype phasing of long reads.
* Finally, execute ```whatshap_haplotag.pl``` to tag phased haplotypes for visualization on IGV.
