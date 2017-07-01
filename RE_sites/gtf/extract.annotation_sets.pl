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

# This script will take a GTF file from Ensembl with five_prime_utr, three_prime_utr and CDS features.
# It will extract these features, introns, promoters and intergenic regions for @fav_biotypes features.
# It requires access to BEDtools and tab delineated chromosome class (Ensembl) and chromosome length 
# (seqimp) files.

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $gtf = "";
my $base = "";
my $chrlenfile="";
my $chrclassfile = "";
my $test = 0;
my $howmany = -1;


# GTF
# Sequence class file
# Output base

if(! GetOptions(
   "gtf=s"  => \$gtf,  
   "base=s" => \$base,
   "len=s"  => \$chrlenfile,
   "class=s"=> \$chrclassfile,
   "test"   => \$test,
   "do=i"   => \$howmany
)){
   die "Could not parse the command line options\n";
}

die "Need an annotation GTF to sort regions\n" if length($gtf)==0;
die "GTF file does not exist\n" if (! -e $gtf);
die "Need a file base name\n" if length($base)==0;
die "Need a chromosome class file\n" if length($chrclassfile)==0;
die "Chromosome class file is missing or empty\n" if (! -s $chrclassfile);
die "Need a chromosome length file\n" if length($chrlenfile)==0;
die "Chromosome length file is missing or empty\n" if (! -s $chrlenfile);

# Keep only biotypes within a set list - noncodnig classes taken from mouse annotation in Ensembl and edited to reflect human gtf (v84).
my @fav_biotypes = ("protein_coding", "3prime_overlapping_ncrna", "antisense", "lincRNA", "non_coding", "processed_transcript", "sense_intronic", "sense_overlapping","macro_lncRNA","bidirectional_promoter_lncrna");
my %fav_biotypes = map {$_ => 0} @fav_biotypes;

# Load chromosome class file

sub read_classes{

   print STDERR "Reading chromosome class file\n";
   my $classfile = shift;
   my %chrClasses;
   open(CLASS, "< $classfile") || die "Could not open the chromosome class file\n";
   my $lineNo =0 ;
   while(<CLASS>){
      chomp;
      $lineNo++;
      if($lineNo==1){
         die "Not expected header\n" if (! /Chr\tClass/);
         next;
      }else{
         my @lineSplit = split "\t";
         die "Chromosome class file line has the incorrect number of fields: $_\n" if scalar(@lineSplit)!=2;
         $chrClasses{$lineSplit[0]} = $lineSplit[1];
      }
   }
   return(\%chrClasses);
}

# Calculate introns and promoters from an array conprised of transcript start and stop coordinates 

sub sort_promoters_and_introns{
   my $feature_chr = shift;
   my $feature_strand = shift;
   my $feature_id = shift;
   my $exon_array_ref = shift;
   my $gtftracker = shift;

   my @sorted_block_exons = sort { $a->[0] <=> $b->[0] } @{$exon_array_ref};

   # Each time a new transcript is seen and $transcriptID changes:
      # Determine first exon according to strand - write upstream promoter region to feature BED file

   if($feature_strand eq "+"){
      print PROMOTER "$feature_chr\t".($sorted_block_exons[0]->[0]-1001)."\t".($sorted_block_exons[0]->[0]-1)."\t$feature_id\t1000\t$feature_strand\n";
   }elsif($feature_strand eq "-"){
      print PROMOTER "$feature_chr\t".($sorted_block_exons[$#sorted_block_exons]->[1])."\t".($sorted_block_exons[$#sorted_block_exons]->[1]+1000)."\t$feature_id\t1000\t$feature_strand\n";
   }else{
      die "Line $gtftracker: transcript strand is unexpected\n";
   }

   # calculate intronic regions.
      # Record last base of each exon, move note first base of next and push  the difference to annotation array as intron
   
   if(scalar(@sorted_block_exons) > 1){
      my $exonNo = 0;
      my $splice_start;
      foreach(@sorted_block_exons){
         $exonNo++;
         if($exonNo == 1){
            $splice_start = $_->[1]
         }else{
            print INTRON "$feature_chr\t$splice_start\t".($_->[0]-1)."\t$feature_id\t1000\t$feature_strand\n"; 
            $splice_start = $_->[1];
         }
      }
   }
}

# Read in the chromosome class file
my $class_hash_ref = &read_classes($chrclassfile);

#print Dumper \$class_hash_ref;

# Only retain chromosome and scaffold classes in the chromosome class hash.
my %class_hash_preserve;
for my $theChr (keys(%{$class_hash_ref})){

   if((${$class_hash_ref}{$theChr} eq "chromosome") || (${$class_hash_ref}{$theChr} eq "scaffold")){
      $class_hash_preserve{$theChr} = ${$class_hash_ref}{$theChr};
   }

}

#print Dumper \%class_hash;

# CHECKED

# Parse Ensembl GTF file

open(GTF, "< $gtf") || die "Could not open the Ensembl GTF for a first pass\n";

my $gtfline=0;
my %allTranIDs;
my $currenttranscriptID = "";
my %chromosomes_classes = map{$class_hash_preserve{$_} => 0} keys %class_hash_preserve;   # Record which chromosomes have annotation
my %geneID_tracker;
my $current_block_line=0;
my $running_gene="";

# print Dumper \%chromosomes_classes;
# CHECKED

print STDERR "First pass: checking transcript annotations are grouped\n";
while (<GTF>){
   $gtfline++;
   # Skip comments
   next if (/^#/);

   # FIRST PASS:
   
   chomp;
   my ($chr,$source,$feature,$start,$end,$score,$strand,$frame,$attribute) = split "\t";
   die "Line $gtfline: expecting 9 fields in each GTF line\n" if (! defined($attribute));

   my $newtranscriptID = "";
   
   # Ensure GTF includes scaffold annotation - skip if chromosome not in preservation type
   if(exists($class_hash_preserve{$chr})){
      $chromosomes_classes{$class_hash_preserve{$chr}}++;
   }else{
      next;
   }

   #print "attribute: $attribute\n";
   
   # Ensure that all transcript exons are sequentially annotated - or die
   if($feature eq "gene"){
      if($attribute=~/gene_id\s"(ENS\w+)";\s/){
         $geneID_tracker{$1}=$gtfline;
      }else{
         die "Line $gtfline: can't identify gene ID\n";
      }

   }elsif($feature eq "exon"){
      if($attribute=~/transcript_id\s"(ENS\w+)";\s/){
         $newtranscriptID = $1;
      }else{
         die "Line $gtfline: can't identify transcript ID\n";
      }


      if ($newtranscriptID eq $currenttranscriptID){
         $current_block_line = $gtfline;
         next;
      }else{

         # Exons from a new transcript
         if(exists($allTranIDs{$newtranscriptID})){
            die "Line $gtfline: transcript IDs split - $newtranscriptID\n";
         }else{
            
            # Each block of transcripts should be accompanied by a 'gene' line
            # Check gene line precedes transcript blocks
            if($attribute=~/gene_id\s"(ENS\w+)";\s/){
               die "Line $gtfline: gene line not yet identified for $newtranscriptID\n" if(! exists($geneID_tracker{$1}));
              
               # Can be multiple transcripts per gene. Only check gene annotation directly precedes when gene ID associated with transcript sets changes. 
               if($running_gene ne $1){
                  die "Line $gtfline: gene line for new transcript block does not directly precede it\n" if ($current_block_line > $geneID_tracker{$1});
               }
               $running_gene=$1;            
   
               $current_block_line = $gtfline;
            }else{
               die "Line $gtfline: can't identify gene ID\n";
            }

            $currenttranscriptID = $newtranscriptID;
            $allTranIDs{$newtranscriptID}++;
         }
      }
   }else{
      next;
   }

}

close(GTF) || die "Could not close the GTF";

# CHECKED

# Check at least one scaffold is recorded
if(! $test){
   foreach(keys %chromosomes_classes){
      die "$_: no annotation found\n" if ($chromosomes_classes{$_} < 1);
   }
}

# Second pass: Record features of interest in an array to print to BED file
print STDERR "Second pass: recording interesting annotations\n";

$gtfline=0;
my %allbiotypes;

my $gene_span_file = "$base.gene_spans.bed";
my $five_file = "$base.five_prime_utr.bed";
my $three_file = "$base.three_prime_utr.bed";
my $cds_file = "$base.cds.bed";
my $promoter_file = "$base.promoter.bed";
my $intron_file = "$base.intron.bed";
my $exon_file = "$base.exon.bed";

print STDERR "Opening files for saving annotation sets\n";
open(GTF, "< $gtf") || die "Could not open the Ensembl GTF for a second pass\n";
open(GENESPAN,">$gene_span_file") || die "Could not open the gene span file for annotation\n";
open(FIVE, "> $five_file") || die "Could not open the 5' UTR file\n";
open(THREE, "> $three_file") || die "Could not open the 3' UTR file\n";
open(CDS, "> $cds_file") || die "Could not open the CDS file\n";
open(PROMOTER, "> $promoter_file") || die "Could not open the promoter file\n";
open(INTRON, "> $intron_file") || die "Could not open the intron file\n";
open(EXON, "> $exon_file") || die "Could not open the exon file\n";

# CHECKED

my $three_utr_count=0;
my $five_utr_count=0;
my $cds_count=0;
my $exon_count=0;

# Variables for tracking exons for each transcript
my $block_transcript ="";
my @block_exons = ();
my $block_strand ="";
my $block_chromosome="";
my $block_gene="";


# SECOND PASS:

while(<GTF>){
   $gtfline++;
   
   last if (($howmany > 0)&&($gtfline > $howmany));
   
   # Skip comments
   next if (/^#/);

   chomp;
   my ($chr,$source,$feature,$start,$end,$score,$strand,$frame,$attribute) = split "\t";
   
   # Throw sequences not on a 'chromosome' or 'supercontig'.
   next if(! exists($class_hash_preserve{$chr}));
   
   # Throw features with an incorrect biotype
   my $biotype = "";
   if($attribute =~ /gene_biotype\s+"(\w+)";/){
      $biotype=$1;
   }else{
      die "Line $gtfline: could not find the gene biotype\n";
   }
   $allbiotypes{$biotype}++;  
   next if (! exists($fav_biotypes{$biotype}));
   $fav_biotypes{$biotype}++;   
 
   # For each annotation save the gene ID and transcript ID
   my $geneID;
   my $transcriptID;
   if($attribute=~/gene_id\s"(ENS\w+)";\s/){
      $geneID=$1;
   }else{
      die "Line $gtfline: can't salvage the gene ID\n";
   }
   
   # Write gene spans to BED file for intergenic calculation 
   if($feature eq "gene"){
      print GENESPAN "$chr\t".($start-1)."\t$end\t$geneID\t1000\t$strand\n";
      next;
   }

   # CHECKED

   # Find transcript ID
   if($attribute=~/transcript_id\s"(ENS\w+)";\s/){
      $transcriptID=$1
   }else{
      die "Line $gtfline: can't salvage the transcript ID\n";
   }
   

   # Each feature type will be written to a seperate BED file
   # Preserve - coding 5utr, 3utr and CDS in BED format.
   # Elements will be named by transcript ID

   (print FIVE "$chr\t".($start-1)."\t$end\t$transcriptID\t1000\t$strand\n")&&($five_utr_count++)&&(next) if ($feature eq "five_prime_utr");
   (print THREE "$chr\t".($start-1)."\t$end\t$transcriptID\t1000\t$strand\n")&&($three_utr_count++)&&(next) if ($feature eq "three_prime_utr");
   (print CDS "$chr\t".($start-1)."\t$end\t$transcriptID\t1000\t$strand\n")&&($cds_count++)&&(next) if ($feature eq "CDS");

   # CHECKED

   #print $feature."\t".$gtfline."\n";
   if($feature eq "exon"){
      $exon_count++;
     
      print EXON "$chr\t".($start-1)."\t$end\t$transcriptID\t1000\t$strand\n";

      #print "Line $gtfline: ($exon_count) $block_transcript - $transcriptID\n";

      if ($block_transcript eq $transcriptID){

         # Check strands and chromosomes are the same for each transcript
         die "Line $gtfline: same transcript must have same gene ID\n" if ($block_gene ne $geneID);
         die "Line $gtfline: same transcript must have same strand\n" if ($block_strand ne $strand);
         die "Line $gtfline: same transcript must have same chromosome\n" if ($block_chromosome ne $chr);
         
         # Push 'exons' to an array for each new transcript
         push @block_exons,[$start,$end];
      
      }else{
         
         # If first exon skip to new transcript block set up
         if($exon_count!=1){
         
            # SORT OLD TRANSCRIPT EXONS:
            # Order exons in the transcript based on start coordinates
           
            &sort_promoters_and_introns($block_chromosome, $block_strand, $block_transcript, \@block_exons, $gtfline);

         }

         # ORGANISE NEW TRANSCRIPT EXON AND SET UP
         
         # Record strand etc. - to ensure all are the same
         $block_strand = $strand;
         $block_chromosome = $chr;
         @block_exons = ();
         
         push @block_exons,[$start,$end];

         $block_transcript = $transcriptID;
         $block_gene = $geneID;
      }
   }
   # Transcript to gene to biotype table to be created
}

# The last transcript
&sort_promoters_and_introns($block_chromosome, $block_strand, $block_transcript, \@block_exons, $gtfline);



#print "\n".$gtfline."\n";

print STDERR "Closing files\n";

close(GTF) || die "Could not close the GTF\n";
close(GENESPAN) || die "Could not close the gene span BED file\n";
close(FIVE) || die "Could not close the 5' UTR file\n";
close(THREE) || die "Could not close the 3' UTR file\n";
close(CDS) || die "Could not close the CDS file\n";
close(PROMOTER) || die "Could not close the promoter file\n";
close(INTRON) || die "Could not close the intron file\n";
close(EXON) || die "Could not close the exon file\n";

# If 5utr, 3utr or CDS annotation missing - die
die("UTR or CDS annotation missing") if (($five_utr_count==0)||($three_utr_count==0)||($cds_count==0));


# Use the genespans file to calculate intergenic regions
print STDERR "Generating intergenic regions file\n";

my $intergenic_file = "$base.intergenic.bed";

#my $intergenic_call = "complementBed -i $gene_span_file -g $chrlenfile > $intergenic_file";
#die "Couldn't find the intergenic regions:\n$intergenic_call\n" if system($intergenic_call);

# Read in the intergenic regions file as they are calculated from standard input and add additional fields before printing to file (BED 6 format).

open (INTER, "complementBed -i $gene_span_file -g $chrlenfile |") || die "Could not open the complementBed connection for intergenic regions\n";
open (INTEROUT, "> $intergenic_file") || die "Could not open the file for the intergenic regions output: $intergenic_file\n";

while(<INTER>){
   chomp;
   print INTEROUT "$_\tintergenic\t.\t*\n";
}

close(INTER) || die "Could not close the intergenic region connection\n";
close(INTEROUT) || die "Could not close the intergenic regions output\n";

#print Dumper \%fav_biotypes;
#print Dumper \%allbiotypes;

##!genome-build GRCh38.p5
###!genome-version GRCh38
###!genome-date 2013-12
###!genome-build-accession NCBI:GCA_000001405.20
###!genebuild-last-updated 2015-10
#1       ensembl_havana  gene    69091   70008   .       +       .       gene_id "ENSG00000186092"; gene_version "4"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTHUMG00000001094"; havana_gene_version "2";
#1       ensembl_havana  transcript      69091   70008   .       +       .       gene_id "ENSG00000186092"; gene_version "4"; transcript_id "ENST00000335137"; transcript_version "3"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTHUMG00000001094"; havana_gene_version "2"; transcript_name "OR4F5-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30547"; havana_transcript "OTTHUMT00000003223"; havana_transcript_version "2"; tag "basic"; transcript_support_level "NA";
#1       ensembl_havana  exon    69091   70008   .       +       .       gene_id "ENSG00000186092"; gene_version "4"; transcript_id "ENST00000335137"; transcript_version "3"; exon_number "1"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTHUMG00000001094"; havana_gene_version "2"; transcript_name "OR4F5-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30547"; havana_transcript "OTTHUMT00000003223"; havana_transcript_version "2"; exon_id "ENSE00002319515"; exon_version "1"; tag "basic"; transcript_support_level "NA";
#1       ensembl_havana  CDS     69091   70005   .       +       0       gene_id "ENSG00000186092"; gene_version "4"; transcript_id "ENST00000335137"; transcript_version "3"; exon_number "1"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTHUMG00000001094"; havana_gene_version "2"; transcript_name "OR4F5-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30547"; havana_transcript "OTTHUMT00000003223"; havana_transcript_version "2"; protein_id "ENSP00000334393"; protein_version "3"; tag "basic"; transcript_support_level "NA";
#1       ensembl_havana  start_codon     69091   69093   .       +       0       gene_id "ENSG00000186092"; gene_version "4"; transcript_id "ENST00000335137"; transcript_version "3"; exon_number "1"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTHUMG00000001094"; havana_gene_version "2"; transcript_name "OR4F5-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30547"; havana_transcript "OTTHUMT00000003223"; havana_transcript_version "2"; tag "basic"; transcript_support_level "NA";
#1       ensembl_havana  stop_codon      70006   70008   .       +       0       gene_id "ENSG00000186092"; gene_version "4"; transcript_id "ENST00000335137"; transcript_version "3"; exon_number "1"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTHUMG00000001094"; havana_gene_version "2"; transcript_name "OR4F5-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30547"; havana_transcript "OTTHUMT00000003223"; havana_transcript_version "2"; tag "basic"; transcript_support_level "NA";
