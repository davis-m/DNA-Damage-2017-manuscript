#########################################################################

# This script is part of a set of tools which aims to facilitate the
# classification of small RNAs in a streamlined and efficient manner.
# Copyright (C) 2016 EMBL - European Bioinformatics Institute

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program (Please see the COPYING file for details).
# If not, see <http://www.gnu.org/licenses/>.

###########################################################################

library("GenomicRanges")

# This script will take a set of comma separated BED files ordered according to hierarchy.
# Redundancy will be removed from each in order. A set of labels of the same length as the input files
# will be used to write the non-redundant spans to file. 

env_path <- Sys.getenv(x="SRNA_PATH")
source(paste(env_path,"/tools/x.general_R_functions.R",sep=""))

args <- R.utils::commandArgs(asValue=TRUE)

if(!is.null(args$out)){
   outDir <- check_argument(args,"out","character")
}else{
   stop("Require and output directory")
}

if(!is.null(args$beds)){
   beds <- check_argument(args,"beds","character")
}else{
   stop("Require a set of comma seperated BED files")
}

if(!is.null(args$labels)){
    labels <- check_argument(args,"labels","character")
}else{
   stop("Require a set of labels the same length as the initial BED files")
}

if(!is.null(args$chrlens)){
    chrlens <- check_argument(args,"chrlens","character")
}else{
   stop("Require a file of chromosome lengths")
}

# Check files and  directies exist
check_file_exists(chrlens)
check_file_exists(outDir)

# Read in the chromosome length file
chrlens <- read_Ensembl_ChrLen(chrlens)
chrlens <- organiseChrLen(chrlens)

# Split beds and labels
all_beds <- unlist(strsplit(beds,","))
all_labels <- unlist(strsplit(labels,","))

# Check there are at least 2 bed files

if(length(all_beds) < 2){
   stop("Require at least 2 BED files to remove redundancy")
}

# Check they are the same length
if(length(all_labels)!=length(all_beds)){
   stop("There must be one label per BED file")
}

# Check labels and files are unique
if(any(duplicated(all_labels))){
   stop("All labels must be unique")
}
if(any(duplicated(all_beds))){
   stop("All bed files must be unique")
}


# Create named vector of BED files
names(all_beds) <- all_labels

rm("all_labels","beds","labels")

# Subroutine for writing
write_intervals <- function(intervals_to_write,all_intervals,output_file){   
   
   # Check GRanges has some ranges in it or warn and next
   if(length(intervals_to_write)==0){
      write("No intervals remain to write to file",file=stderr())
      return()
   }else{
      
      # Check no overlap with all_ranges
      overlap_count <- countOverlaps(intervals_to_write,all_intervals)
      if(any(overlap_count > 0)){
         stop("Redundant intervals remain in file to write")
      }
      
      # Print output file
      write(paste("\nOutput file: ",output_file,"\n",sep=""),file=stderr())

      # Output file metrics to include: max and min widths, number of ranges
      max_width <- max(width(intervals_to_write))
      min_width <- min(width(intervals_to_write))
      number_of_intervals <- length(intervals_to_write)

      write(paste("Remaining interval metrics:\nmax width - ",max_width,"\nmin width - ",min_width,"\nnumber of intervals remaining - ",number_of_intervals,"\n\n",sep=""),file=stderr())

      # Write remaining ranges to file.
      writeGRtoBED(output_file,intervals_to_write)      
   }
}


# For each:
bed_number <- 0
all_ranges <- GRanges()

for(bed_set in names(all_beds)){
   bed_number <- bed_number+1   
   
   write(paste("### Removing redundancy:",bed_set),file=stderr())

   # Check BED file exists
   check_file_exists(all_beds[bed_set])   

   # Read in BED file
   intervals <- read_bed(all_beds[bed_set])
   write(paste("Number of intervals:",nrow(intervals)),file=stderr())

   # Convert BED to GR
   intervals[intervals$strand==".", "strand"] <- "*"
   intervals_GR <- BEDtoGR_V3(intervals,chrlens,stranded=TRUE,element_names=intervals$name)   
   rm("intervals")

   # Create output file name from label and outDir
   outfile <- paste(outDir,"/",bed_set,".no_redundancy.bed",sep="")

   # Check that the output file does not exist
   if(file.exists(outfile)){
      stop(paste(outfile,"already exists"))
   }
   
   # CHECKED   

   # If first create all_ranges object
   if(bed_number == 1){
      
      mcols(intervals_GR)$names <- paste(bed_set,".",seq(1,length(intervals_GR)),sep="")
         
      # Write all intervals to non-redundant BED file
      write_intervals(intervals_GR,all_ranges,outfile)
      
      # Add ranges to all_ranges
      all_ranges <- c(all_ranges,intervals_GR)

   }else{
      # If subsequent BED file
      # Find diff between all_ranges and file -> diff_ranges object
      subset_GR <- setdiff(intervals_GR,all_ranges,ignore.strand=FALSE)
      
      # Add back names
      mcols(subset_GR)$names <- paste(bed_set,".",seq(1,length(subset_GR)),sep="")

      # Write diff_ranges to non-redundant BED file
      write_intervals(subset_GR,all_ranges,outfile)
      
      # Add ranges to all_ranges
      all_ranges <- c(all_ranges,intervals_GR)
   }
}
write("Done",file=stderr())


sessionInfo()
