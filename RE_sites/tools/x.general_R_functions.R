###########################################################################

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

# Function to write a GRanges object with a 'names' metadata column to file
writeGRtoBED <- function(outFile=c(),theGR=GRanges()){
   if(length(outFile)==0){
      stop("writeGRtoBED: require an 'outFile'")
   }
   if(length(theGR)==0){
      stop("writeGRtoBED: require 'theGR'")
   }
   if (!is(theGR,"GRanges")){
      stop("writeGRtoBED: 'theGR' must be 'GRanges' class")
   }
   if(! "names" %in% colnames(mcols(theGR))){
      stop("theGR requires a names metadata column")
   }

   if(file.exists(outFile)){
      stop(paste(outFile,"already exists"))
   }

   outDF <- as.data.frame(theGR)
   outDF$score <- 1000

   # 1 based coordinates to 0-based end-open
   outDF$start <- outDF$start-1

   # * not defined in BED format
   outDF$strand <- as.character(outDF$strand)
   outDF$seqnames <- as.character(outDF$seqnames)

   if(any(outDF$strand=="*")){
      outDF[outDF$strand=="*","strand"] <- "."
   }

   # Write table to file
   write.table(outDF[c("seqnames","start","end","names","score","strand")],file=outFile,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}


# This function will take a seqimp GRanges converted alignment set and split it into a set of unique and multimapped alignments
split_by_unique_alignments <- function(read_GR){

   if(class(read_GR) != "GRanges"){
      stop("read_GR must be a GRanges object")
   }
   if(! ("locations" %in% names(values(read_GR)))){
      stop("Need a locations column in the GR metadata")
   }
   if(!is.integer(values(read_GR)[,"locations"])){
      stop("The locations column must be a set of integers")
   }

   read_list <- GRangesList()
   read_list[["single_align"]] <- read_GR[values(read_GR)$locations == 1]
   read_list[["multi_align"]] <- read_GR[values(read_GR)$locations > 1]
   read_list[["all_align"]] <- read_GR

   return(read_list)
}


# Functions for checking command line options and the existance of files
check_argument <- function(argument_list,argument_name,arg_class){
   if(!is.null(argument_list[[argument_name]])){
      arg_return <- argument_list[[argument_name]]
   }else{
      stop(paste("Require argument: --",argument_name,sep=""))
   }
   arg_return <- as(arg_return,arg_class)
   return(arg_return)
}

check_file_exists <- function(file_to_check){
   if(!file.exists(file_to_check)){
      stop(paste("Missing file or directory:",file_to_check,sep="\n"))
   }
}


# Load a sample description table
load_sample_description <- function(fileName,expected_headers,nameCol=""){
   # fileName - path to sample table file
   # expected_headers - required headers to check
   # nameCol - column to use as the rownames

   if(!(file.exists(fileName))){
      stop(paste("File does not exist: ",fileName,sep=""))
   }
   temp_table <- read.table(file=fileName, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   if(! all(expected_headers %in% colnames(temp_table))){
      stop(paste("Some expected header names are missing: ",fileName,sep=""))
   }

   if(nchar(nameCol)>0){
      if(!nameCol %in% expected_headers){
         stop(paste("Requested name column not amongst those expected: ",fileName," - ",nameCol, sep=""))
      }else if(any(duplicated(temp_table[,nameCol]))){
         stop("Duplicated names in nameCol")
      }else{
         rownames(temp_table) <- temp_table[,nameCol]
      }
   }

   return(temp_table)
}

# Convert a GR object to its antisense
antisense_me <-  function(sense_version){
   if(any(!strand(sense_version) %in% c("+","-"))){
      stop("Unexpected lack of a strand")
   }
   antisense_version <- sense_version
   pos_strands <- strand(antisense_version) == "+"
   strand(antisense_version[pos_strands]) <- "-"
   strand(antisense_version[!pos_strands]) <- "+"
   return(antisense_version)
}

## usage: foo<-magicLoad("toload1.RData")
## usage: bar<-magicLoad("toload2.RData","aaa")
# Function provided by Oliver Davis
magicLoad<-function(file,pattern=""){
        env<-new.env()
        load(file,envir=env)
        length(ls(env,pattern=pattern)) > 0 || stop("cannot find an object matching ",pattern)
        return(eval(as.name(ls(env,pattern=pattern)[1]),envir=env))
}


# Read in a 6 column BED file
read_bed <-function(bed_file){
   bed_read <- read.table(bed_file,sep="\t",header=FALSE,as.is=TRUE)
   colnames(bed_read) <- c("seq","start","end","name","score","strand")
   return(bed_read)
}


# In version 2 I have added a name option to add a name to each element in the BED file
# In version 3 I have made seqnameConversion optional
BEDtoGR_V3 <- function(BEDRData,lengthFun,stranded=FALSE,classif="raddish",element_names=c(),convert=FALSE){
   if(missing(BEDRData)){stop("BEDtoGR: requires BED data")}
   if(missing(lengthFun)){stop("BEDtoGR: requires chromosome length data")}
   if((!"seq"%in%colnames(BEDRData))|(!"start"%in%colnames(BEDRData))|(!"end"%in%colnames(BEDRData))){
      stop("BEDtoGR: BED data must have 'seq', 'start' and 'end' headers")
   }
   if(convert){
      BEDRData$seq <- seqnameConverter(BEDRData$seq)
   }
   if(!stranded){
     if(classif!="raddish"){
         propGR  <-   GRanges(
            seqnames       = Rle(factor(BEDRData$seq,levels=names(lengthFun))),
            ranges         = IRanges(start=(BEDRData$start+1), end=BEDRData$end),
            strand         = "*",
            seqlengths     = lengthFun,
            classification = rep(classif,nrow(BEDRData))
         )
      }else{
         propGR  <-   GRanges(
            seqnames       = Rle(factor(BEDRData$seq,levels=names(lengthFun))),
            ranges         = IRanges(start=(BEDRData$start+1), end=BEDRData$end),
            strand         = "*",
            seqlengths     = lengthFun
         )
      }
   }else{
     if(!"strand"%in%colnames(BEDRData)){stop("BEDtoGR: Bed data must have 'strand' header")}
     if(classif!="raddish"){
         propGR  <-   GRanges(
            seqnames       = Rle(factor(BEDRData$seq,levels=names(lengthFun))),
            ranges         = IRanges(start=(BEDRData$start+1), end=BEDRData$end),
            strand         = BEDRData$strand,
            seqlengths     = lengthFun,
            classification = rep(classif,nrow(BEDRData))
         )
      }else{
         propGR  <-   GRanges(
            seqnames       = Rle(factor(BEDRData$seq,levels=names(lengthFun))),
            ranges         = IRanges(start=(BEDRData$start+1), end=BEDRData$end),
            strand         = BEDRData$strand,
            seqlengths     = lengthFun
         )
      }
   }
   if((length(element_names) != 0)){
      if(length(element_names) != length(propGR)){
         stop("Names passed to BEDtoGR_V2 must be the same length as number of elements in the BED file")
      }else{
         values(propGR)$names <- element_names
      }
   }
   return(propGR)
}

# Read a UCSC chr length file - created from 2bit file (No header)
read_UCSC_ChrLen <- function(chrLenFile=""){
   if (nchar(chrLenFile)==0){stop("No chromosome length file passesd")}
   chromLengthTab <- read.table(chrLenFile, sep="\t", as.is=TRUE, header = FALSE,colClasses = c("character", "integer"))
   colnames(chromLengthTab) <- c("Name","Length")
   return(chromLengthTab)
}

# Ensembl Chr length file - created by SequenceImp annotation
read_Ensembl_ChrLen <- function(chrLenFile=""){
   if (nchar(chrLenFile)==0){stop("No chromosome length file passesd")}
   chromLengthTab <- read.table(chrLenFile, sep="\t", as.is=TRUE, header = TRUE,colClasses = c("character", "integer"))
   return(chromLengthTab)
}  


### Organise SequenceImp annotation chromosome length table into a named list (v12-028)

organiseChrLen <- function(chrTable){
   if(!is.data.frame(chrTable)){stop("No chromosome length table specified")}
   chrTable <- chrTable[order(chrTable[,1]),]
   tempLens <- chrTable[,2]
   names(tempLens) <- chrTable[,1]
   return(tempLens)
}

# Convert UCSC identifiers into Ensembl IDs
# chrNames - a vector of chromosome identifiers
seqnameConverter <- function(chrNames){
   newSequences <- sub("^chr","",chrNames,ignore.case=TRUE,perl=TRUE)
   newSequences <- sub("^M$","MT",newSequences,ignore.case=TRUE,perl=TRUE)
   return(newSequences)
}

# Read in a typical GTF file
# GTFfile - The path to a GTF file
read_gtf <- function(GTFfile){
   check_file_exists(GTFfile)
   aGTFtable <- read.table(file=GTFfile, header=FALSE, sep="\t",colClasses=c("character","character","character","integer","integer","character","character","character","character"))
   colnames(aGTFtable) <- c("seq","source","feature","start","end","score","strand","frame","attribute")
   return(aGTFtable)
}


