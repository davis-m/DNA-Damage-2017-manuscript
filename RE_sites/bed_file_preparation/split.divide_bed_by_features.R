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

# The script will take a BED file, find overlaps with BED file annotations in a
# strand independent manner and then split to multiple output files depending
# upon underlying annotations. Output files will be redundant.

# TO DO:
# Make strand specificity optional

library("GenomicRanges")
library("tools")

env_path <- Sys.getenv(x="SRNA_PATH")
source(paste(env_path,"/tools/x.general_R_functions.R",sep=""))

strand_specific <- FALSE

args <- R.utils::commandArgs(asValue=TRUE)

inputBED <- check_argument(args,"inbed","character")           # Input BED file to be divided
annotDir <- check_argument(args,"annotdir","character")        # Annotation directory
annotBEDbase <- check_argument(args,"annotstem","character")   # Common stem for the names of all annotation files to be used for the comparison
outDir   <- check_argument(args,"out","character")             # Ouput directory
lenFile <- check_argument(args,"len","character")              # Matching chromosome length file, 2 columns, header - chromosome name tab length.

# Interpret strand specificity flag

if(! is.null(args$strand_specific)){
   stopifnot(is.logical(args$strand_specific))
   strand_specific <- args$strand_specific                     # Flag to determine if strand information should be used: Default - No
   print("Using strand specific overlaps")
}else{
   print("Working without strand information")
}

# CHECKED

# Read in the chromosome length file
check_file_exists(lenFile)
chrLen <- read_Ensembl_ChrLen(lenFile)
chrLen <- organiseChrLen(chrLen)

# Read in the input BED file and convert to a GR object

check_file_exists(inputBED)
inBED <- read_bed(inputBED)

# BED file does not accept * in strand column so convert here.
if(any(inBED$strand==".")){
   inBED[inBED$strand==".","strand"] <- "*"
}

if(!all(inBED$strand %in% c("+","-","*"))){
   stop("Unrecognised strand info in input BED file\n")
}

inBEDGR <- BEDtoGR_V3(inBED,chrLen,stranded=TRUE,element_names=inBED$name)

# Strand specificity: Convert all strands to "*"? - strand independent - sort with the findOverlaps command - need to retain info for output
# strand(inBEDGR) <- "*"

# input file base to conatenate for output

inBase <- basename(inputBED)
inStem <- file_path_sans_ext(inBase)


# Collect all annotation files with the same base name for comparisons

allAnnotBEDs <- list.files(path=annotDir,pattern=paste(annotBEDbase,".*\\.bed",sep=""))
if(length(allAnnotBEDs)<1){
   stop(paste("No annotation BED files found in ",annotDir," with pattern '",annotBEDbase,".*\\.bed'",sep=""))
}

# CHECKED

# Cycle through annotation files
for (annotType in allAnnotBEDs){
   print(paste("Finding overlaps: ",annotType))

   # Organise base for the annot files (remove .bed)
   annotTypeStem <- sub("(^[\\w\\-\\.]+).bed","\\1",annotType,perl=TRUE)
   
   # Organise output name - concatenate file basename to the output base
   if(strand_specific){
      outputFile <- paste(outDir,"/",inStem,".",annotTypeStem,".strand_spec.bed",sep="")
   }else{
      outputFile <- paste(outDir,"/",inStem,".",annotTypeStem,".strand_non_spec.bed",sep="")
   }

   if(file.exists(outputFile)){
      stop(paste("Output file already exists:",outputFile))
   }

   #print(annotTypeBase)

   # Read in each annot file
   annotBED <- read_bed(paste(annotDir,annotType,sep="/"))

   # BED file does not accept * in strand column so convert here.
   if(any(annotBED$strand==".")){
      annotBED[annotBED$strand==".","strand"] <- "*"
   }

   if(!all(annotBED$strand %in% c("+","-","*"))){
      stop("Unrecognised strand info in annotation BED file\n")
   }

   annotBEDGR <- BEDtoGR_V3(annotBED,chrLen,stranded=TRUE,element_names=annotBED$name)
   rm("annotBED")

   # Convert to BED format
   if(strand_specific){
      # Find overlapping inputBED sites (strand specific)
      featureOverlaps <- findOverlaps(inBEDGR,annotBEDGR,ignore.strand=FALSE)
   }else{
      # Find overlapping inputBED sites (strand non-specific)
      featureOverlaps <- findOverlaps(inBEDGR,annotBEDGR,ignore.strand=TRUE)
   }
   
   # Define new element names by concatenating element names with transcript IDs
   
   feature_overlap_sets <- as.data.frame(featureOverlaps)
   if(nrow(feature_overlap_sets) > 0){
      feature_overlap_sets$subjectHits <- mcols(annotBEDGR)[feature_overlap_sets$subjectHits,"names"]
      feature_overlap_sets <- split(feature_overlap_sets$subjectHits,feature_overlap_sets$queryHits)
      feature_overlap_sets <- sapply(feature_overlap_sets,function(x){paste(unique(x),collapse="//")})

      overlappedGRanges <- inBEDGR[as.integer(names(feature_overlap_sets))]
      mcols(overlappedGRanges)$overlapped <- feature_overlap_sets
      mcols(overlappedGRanges)$names <- paste(mcols(overlappedGRanges)$names,mcols(overlappedGRanges)$overlapped,sep=".")
   
      # Write overlapping sites to BED file
      writeGRtoBED(outputFile,overlappedGRanges)

   }else{
      print(paste("No oerlaps for",annotType))
   }
   rm("annotBEDGR")
}

sessionInfo()

# CHECKED

