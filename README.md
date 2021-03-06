# Vahine

Scripts for Vahine metagenomics and metatranscriptomics:

* `generate_ACT_comparison_crunch.sh`: given two genome in fasta format, generate the ACT comparison file.

* `reorder_gbk_according_to_fna.py`: given a ordered fna file (usually generated by Mauve), order the genbank file based on the contig order in the fna file, which can be loaded into Mauve to do genome alignment. 

* `download_extract_virus_fna.sh`: automaticlly pull all virus genomes from NCBI, then extract the fna files.  

* `extract_pe_fastq_by_contigs.py`: given a SAM file, a subset of reference contigs in fasta format, fwd and rev fastq reads, extract the read pairs aligned to the subset contigs.

* `slice_bam_by_contigs.sh`: given a BAM file and a subset of reference contigs in fasta format, slice bam to keep records aligned to the reference contigs.

* `calc_cov.sh`: given a subset of reference contigs in fasta format, the forward and reverse reads in fastq format, calculate the read coverage for each contig.

* `calc_cov_in_bam_depth.py`: script to calculate coverage of each contig in a samtools depth file. It was embedded into `calc_cov.sh`, but can also be used standalone.  
