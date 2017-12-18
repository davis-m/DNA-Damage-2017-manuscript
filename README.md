Coverage plots and restriction enzyme sites

coverage_plots/
This directory contains scripts used to coordinate the analysis of aligned
reads to produce comparative coverage plots.

align.Organise_Bowtie2_alignments.pl 
This script will take a set of FASTQ files and run Bowtie2 to align each to a
genome. Alignments are performed on a server farm.

align.sam_folder_to_bam.pl
This script coordinates the conversion of sam files to bam files.

filter.bam_files.sh
This script will apply a series of filters to a directory of bam files,
including the removal of mapped duplicates and low quality.  NOTE: all read
manipulation is performed in tmp/.

sub.rename_reads.pl
This script reformats read IDs.

sub.sam_to_bam_streamlined.pl
Takes a sam file and uses samtools to convert this to a sorted and indexed bam
file.

sub.submit_Bowtie2_jobs.sh
This script organises alignments with Bowtie2 in tmp/.

coverage.create_relative_heatmap.pl
This script will accept a description file containing BAM files for comparison
and it will convert these BAMs into coverage files and subsequently counts
matrices using deeptools, in regions specified in an annotation description
file. Finally the matrices will be combined into regional profiles.

dtCompare.R
This script was kindly donated by Tommaso Leonardi. This script reads matrices
produced by deeptools from bigWig files and generates multiple heatmaps and
profiles in a single plot. Please see tleonardi on GitHub for the full
repository.

rename_standard_tally_output.pl
This script will take a FASTQ file from the standard SequenceImp output. It
will then replace the read IDs.  This will ensure IDs in paired files match
(removes trinucleotide score). Read repeat numbers will be retained in the
format (x\d+$).

Useful_tools.pm
A perl module containing useful subroutines.



RE_sites/
This directory contains the scripts used to explore restriction enzyme sites.

compare_AsiSI_sites_to_general_set.R
This script takes the genome wide AsiSI site set and compares it to the
active (targeted) subset.

bed.remove_redundancy.R
This script will take a set of comma separated BED files ordered according to
hierarchy.  Redundancy will be removed from each in order. A set of labels of
the same length as the input files will be used to write the non-redundant
spans to file. 

split.divide_bed_by_features.R
The script will take a BED file, find overlaps with BED file annotations in a
strand independent manner and then split to multiple output files depending
upon underlying annotations. Output files will be redundant.

split.divide_bed_by_single_feature_overlap.pl
This script will take 2 bed files of features and separate the first BED file
into elements that overlap the features in the second and those that don't. A
base name for output is also required. Features in each output will also be
merged to remove overlaps Merging is non-strand specific and collapses to 3
columns.

extract.annotation_sets.pl
This script will take a GTF file from Ensembl with five_prime_utr,
three_prime_utr and CDS features.  It will extract these features, introns,
promoters and intergenic regions for @fav_biotypes features.  It requires
access to BEDtools and tab delineated chromosome class (Ensembl) and chromosome
length (seqimp) files.

motifs.find_RE_motif_coordinates.pl
This script will search for specified motifs in a fasta file.

sample.random_lines.sh
This script will take a file, a number of lines to randomly sample and a
number of samples to take. It will then produce random subsamples.

x.general_R_functions.R
Contains a set of useful R functions.



