#!/usr/local/bin/perl

###########################################################################
#
# This script is part of a set of tools which aims to facilitate the
# classification of small RNAs in a streamlined and efficient manner.
# Copyright (C) 2016 EMBL - European Bioinformatics Institute
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program (Please see the COPYING file for details).
# If not, see <http://www.gnu.org/licenses/>.
#
############################################################################


# This script will accept a description file containing BAM files for comparison and it
# will convert these BAMs into coverage files and subsequently counts matrices using deeptools,
# in regions specified in an annotation description file. Finally the matrices
# will be combined into regional profiles.

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename "dirname";
use Cwd "abs_path";

# SUBROUTINES
use lib dirname(dirname(abs_path($0)))."/lib";
use Small_rna::Useful_tools "read_desc";
use Small_rna::Useful_tools "read_option_file";
use Small_rna::Useful_tools "pipesystem";
use Small_rna::Useful_tools "exclude_options";

# TO DO: Sort criteria and sam sort order requirements to be defined.
# Check no conflict in options
# Add scaled regions option
# Add option to skip bigwig creation


my $help = 0;

# Tom and Gio's Heatmap plotter
my $dtCompare = "$ENV{SRNA_PATH}/coverage/dtCompare.R";

# Inputs will be a directory.
my $inDir = "";

# Output will also be a directory.
my $outDir = "";

# A description file listing names and the BAM file to find in the input directory
my $description = "";

# BAM coverage options file.
my $bigwig_coverage_options = "";

# BAM compare options file.
my $bigwig_compare_options = "";

# Compute matrix options file.
my $matrix_options = "";

# Heatmap optionf file.
my $heatmap_options = "";

# Regions intended for comparison.
# The directory
my $regionDir = "";
# The files
my $regionDesc = "";

# A regional BED file for plotting.

if (! GetOptions(
   "help"                           => \$help,
   "--in=s"                         => \$inDir,
   "--out=s"                        => \$outDir,
   "--desc=s"                       => \$description,
   "--bigwig-coverage-options=s"    => \$bigwig_coverage_options,
   "--bigwig-compare-options=s"     => \$bigwig_compare_options,
   "--matrix-options=s"             => \$matrix_options,
   "--heatmap-options=s"            => \$heatmap_options,
   "--regions-dir=s"                => \$regionDir,
   "--regions-desc=s"               => \$regionDesc
)){
   die "Could not parse options\n";   
}

&help_message() if $help;

sub help_message{
   <<EOH;
   --help            <FLAG>   This message
   --in              <DIR>    Input directory containing BAM files
   --out             <DIR>    Directory for outputs
   --desc            <FILE>   Description file containing a sample name and BAM file for each sample
   --bigwig-compare-options   <FILE>   A file containing all of the options for bamCompare for sample vs. control analysis
   --bigwig-coverage-options  <FILE>   A file containing all of the options for bamCoverage for single file analysis
   --matrix-options  <FILE>   A file containing all of the options for computeMatrix
   --heatmap-options <FILE>   A file containing all of the options for dtCompare
   --regions-dir     <DIR>    A directory containing the regions for the relative computation of coverage
   --regions-desc    <FILE>   A list of annotation BED files in the regions directory to be used of the comparison with a 'name' for each
EOH
}

die "Require input directory (--in)\n"                            if length($inDir)==0;
die "Require output directory (--out)\n"                          if length($outDir)==0;
die "Require description file (--desc)\n"                         if length($description)==0;
#die "Require file of bamCoverage options (--bigwig-options)\n"    if length($bigwig_options)==0;
die "Require file of computeMatrix options (--matrix-options)\n"  if length($matrix_options)==0;
die "Require file of dtCompare.R options (--heatmap-options)\n"   if length($heatmap_options)==0;
die "Require regions directory containing annotation BED files (--regions-dir)\n" if length($regionDir)==0;
die "Require a file listing all interesting annotation BED files in the regions-dir with two columns: unique 'name' and 'bed (--regions-desc)\n'" if length($regionDesc)==0;

die "Input directory does not exist\n" if (! -d $inDir);
die "Output directory does not exist\n" if (! -d $outDir);
die "Description file empty or does not exist\n" if (! -s $description);
die "computeMatrix options file does not exist\n" if (! -e $matrix_options);
die "dtCompare.R options file does not exist\n" if (! -e $heatmap_options);
die "Regions directory does not exist\n" if (! -d $regionDir);
die "Regions description file does not exist or is empty\n" if (! -s $regionDesc);


# Description file read in

my $description_ref = &read_desc($description,["name","bam_1","bam_2"],"name");

# Read in annotation description file

my $regions_ref = &read_desc($regionDesc,["name","bed"],"name");


# bamCoverage options file read in
my $bam_cov_options = "";
if(length($bigwig_coverage_options)!=0){
   die "bamCoverage options file does not exist (--bigwig-coverage-options)\n" if (! -e $bigwig_coverage_options);
   print STDERR "Reading bamCoverage options: $bigwig_coverage_options\n";
   $bam_cov_options = &read_option_file($bigwig_coverage_options);
   &exclude_options($bam_cov_options,["-b", "--outFileName", "--outFileFormat"]);
}

# bamCompare options file read in
my $bam_comp_options = "";
if(length($bigwig_compare_options)!=0){
   print STDERR "Reading bamCompare options: $bigwig_compare_options\n";
   die "bamCompare options file does not exist (--bigwig-compare-options)\n" if (! -e $bigwig_compare_options);
   $bam_comp_options = &read_option_file($bigwig_compare_options);
   &exclude_options($bam_comp_options,["-b1","-b2","--outFileName","--outFileFormat"]);
}

# compute matrix options file read in

print STDERR "Reading computeMatrix options: $matrix_options\n";
my $count_matrix_options = &read_option_file($matrix_options);
&exclude_options($count_matrix_options,["--reference-point","--scoreFileName","--outFileName","--outFileNameMatrix","--regionsFileName","--sortRegions","--missingDataAsZero"]);

# heatmap options file read in

print STDERR "Reading heatmap options: $heatmap_options\n";
my $plot_heatmap_options = &read_option_file($heatmap_options);
&exclude_options($plot_heatmap_options,["vanilla","--files","--labels","--outFile","--outTable"]);

#exit;

# Loop through BAM files in description file

# Compute single Bigwig but compare to multiple BED files to create multiple heatmaps.
   
# Create Bigwig file for each sample in the description file
# Rotate through the samples and for each generate a BIGWIG file
# For each sample ID record the Bigwig to a hash

my %bigwig_record;
for my $sampleID (keys %{$description_ref}){
   
   print STDERR "\n######### Preparing Bigwig: $sampleID\n";
   #$loopcount++;


   # Design sample specific Bigwig and using names in description file
   my $bigwigOut = "$outDir/".$sampleID.".bw";
   die "$bigwigOut already exists\n" if (-e $bigwigOut);


   # Sorted at heatmap stage
   # Also file for sample order from first countmatrix
   # my $count_ordered_bed = "$outDir/".$sampleID.".matrix_reordered_bed.bed";
   # die "$count_ordered_bed already exists\n" if (-e $count_ordered_bed);

   # Run bigwig conversion
   
   # Check they are in the input directory
   my $inBAM = "$inDir/".${$description_ref}{$sampleID}{"bam_1"};
   die "$inBAM can't be found\n" if (! -s $inBAM);


   if (${$description_ref}{$sampleID}{"bam_2"} eq "-"){
      die "Require a bamCoverage options file to be specified for single bam analysis (--bigwig-coverage-options)\n" if length($bigwig_coverage_options)==0;
      my $bigwig_call = "bamCoverage -b $inBAM --outFileName $bigwigOut $bam_cov_options --outFileFormat bigwig";
      print STDERR "\nBigWig conversion:\n$bigwig_call\n";
      die "BigWig conversion failed\n" if &pipesystem($bigwig_call);
      die "Can't find conversion call output: $bigwigOut\n" if (! -s $bigwigOut);
   }else{
      die "Require a bamCompare options file to be specified for input vs. control analysis (--bigwig-compare-options)\n" if length($bigwig_compare_options)==0;
      my $inBAM2 = "$inDir/".${$description_ref}{$sampleID}{"bam_2"};
      die "$inBAM2 can't be found\n" if (! -s $inBAM2);
      my $bigwig_call = "bamCompare -b1 $inBAM -b2 $inBAM2 --outFileName $bigwigOut $bam_comp_options --outFileFormat bigwig";
      print STDERR "\nBigWig conversion:\n$bigwig_call\n";
      die "BigWig conversion failed\n" if &pipesystem($bigwig_call);
      die "Can't find conversion call output: $bigwigOut\n" if (! -s $bigwigOut);
   }

   $bigwig_record{$sampleID} = $bigwigOut;
}

#print Dumper \%bigwig_record;

#exit;

# CHECKED

# Push each bigwig to a hash with a sample name key
# For each annotation BED file in the annot directory compute the matrix for all samples in the hash
# Matrix output base will be file stem concatenated to annotation file stem
# Create a heatmap for all sets
# Each heatmap base will be annotation file stem




# ADD NEW LOOP HERE FOR LOOPING THROUGH SEPERATE ANNOTATIONS
for my $annot (keys %{$regions_ref}){ 
   print STDERR "\n######### Generating heatmap for new annotation: $annot\n\n";

   # Check each annotation file exists
   my $regions_file = "$regionDir/".$regions_ref->{$annot}{bed};
   die "$regions_file does not exist or is empty\n" if (! -s $regions_file);
   
   # Create output and tracking files for each annotation pass
   my $label_record = "";
   my $matrix_record = "";
   my $heatmap_plot_base = "$outDir/".join(".", keys %bigwig_record).".$annot" ;
   
   die "$heatmap_plot_base"."_profile.pdf already exists\n" if (-e ("$heatmap_plot_base"."_profile.pdf"));   
   die "$heatmap_plot_base"."_heatmap.pdf already exists\n" if (-e ("$heatmap_plot_base"."_heatmap.pdf"));   
   
   my $heatmap_table= "$heatmap_plot_base.heatmap.txt";
   die "$heatmap_table already exists\n" if (-e $heatmap_table);   

   print STDERR "Output base: $heatmap_plot_base\n"; 

   for my $sample_choice(keys %bigwig_record){
      my $countmatrix = "$outDir/".$sample_choice.".$annot.count_matrix.txt.gz";
      die "$countmatrix already exists\n" if (-e $countmatrix);   
      my $R_parsable_counts = "$outDir/".$sample_choice.".$annot.count_matrix.Rparse.tab";
      die "$R_parsable_counts already exists\n" if (-e $R_parsable_counts);   
      my $output_region_order = "$outDir/".$sample_choice.".$annot.region_order.bed";
      die "$output_region_order already exists\n" if (-e $output_region_order);   
   
      # Compute matrix
      my $computeMatrix_call = "computeMatrix  reference-point  --scoreFileName $bigwig_record{$sample_choice} --outFileName $countmatrix --outFileNameMatrix $R_parsable_counts --outFileSortedRegions $output_region_order --regionsFileName $regions_file --sortRegions keep --missingDataAsZero $count_matrix_options";
      
      # NOTE: Now sorted by Tom's script 080716
      
      print STDERR "\nCreating count matrix for $sample_choice:\n$computeMatrix_call\n";
      die "Could not create the countMatrix:\n" if &pipesystem($computeMatrix_call);
      die "Could not find the count matrix file for R: $R_parsable_counts\n" if (! -s $R_parsable_counts);
      die "Could not find the file containing the region order: $output_region_order\n" if (! -s $output_region_order);
   
      $label_record = "$label_record\t$sample_choice";
      $matrix_record = "$matrix_record\t$R_parsable_counts";
   
      # Now plotted and sorted by Tom's script
   }
   
   # Now plotted and sorted by Tom's script
   print STDERR "\nCreating the coverage heatmap and combined profile:\n";
   my $plotHeatMap_call = "Rscript --vanilla $dtCompare --files $matrix_record --labels $label_record --outFile $heatmap_plot_base --outTable $heatmap_table $plot_heatmap_options";
   print STDERR "$plotHeatMap_call\n";
   die "Could not create the coverage heatmap\n" if &pipesystem($plotHeatMap_call);
   
   #last;
}

print STDERR "\nDone\n";
# CHECKED

