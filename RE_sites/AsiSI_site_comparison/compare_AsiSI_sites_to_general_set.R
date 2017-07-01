###########################################################################

# This script was written to compare 2 sets of AsiSI sites to each other.
# Copyright (C) 2017 EMBL - European Bioinformatics Institute

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

# This script takes the genome wide AsiSI site set and compares it to the
# active (targeted) subset.


library("GenomicRanges")
library("RColorBrewer")

setwd("") # Needs to be set


# APPEND FULL PATHS:

all_AsiSI_sites <- "AsiSI_sites.GRCh38.prim_assembly.bed"
repaired_AsiSI_sites <- "AsiSI_sites_repair.bed"
annotation_set_dir <- "formatted_annotation/"
chrLenFile = "human_1.chrom_len.tab"
outDir = "Active_AsiSI_site_sets/"

# Subroutines
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

# Read in a 6 column BED file
read_bed <-function(bed_file){
   bed_read <- read.table(bed_file,sep="\t",header=FALSE,as.is=TRUE)
   colnames(bed_read) <- c("seq","start","end","name","score","strand")
   return(bed_read)
}

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

read_Ensembl_ChrLen <- function(chrLenFile=""){
   if (nchar(chrLenFile)==0){stop("No chromosome length file passesd")}
   chromLengthTab <- read.table(chrLenFile, sep="\t", as.is=TRUE, header = TRUE,colClasses = c("character", "integer"))
   return(chromLengthTab)
}

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

# Organise the chromosome length information
chrlens <- read_Ensembl_ChrLen(chrLenFile)
chrlens <- organiseChrLen(chrlens)

# Read in the annotation information
annotation_files <- list.files(annotation_set_dir,pattern="homo_sapiens.GRCh38.ens-v84*",full.names=TRUE)

annotList <- GRangesList()
for (annotSet in annotation_files){
      annotBED <- read_bed(annotSet)
      dataset_name <- gsub(".*homo_sapiens.GRCh38.ens-v84\\.(\\w+)\\.bed","\\1",annotSet,perl=TRUE)
      if(any(names(annotList)==dataset_name)){
         stop("Repeated dataset")
      }
      annotList[[dataset_name]] <- BEDtoGR_V3(annotBED,chrlens,stranded=TRUE) 
}

# CHECKED


# Read in all the AsiSI sites
all_AsiSI <- read_bed(all_AsiSI_sites)
allAsiSIGR <- BEDtoGR_V3(all_AsiSI,chrlens,element_names=paste("all_sites",seq(1,nrow(all_AsiSI)),sep="_"))
names(allAsiSIGR) <- allAsiSIGR$names

# Read in sites that are targetted
new_targetted_sites <- read_bed(repaired_AsiSI_sites)
targettedAsiSIGR <- BEDtoGR_V3(new_targetted_sites,chrlens,element_names=paste("targetted_sites",seq(1,nrow(new_targetted_sites)),sep="_"))
names(targettedAsiSIGR) <- targettedAsiSIGR$names

# Pair targetted sites with a site within the 'whole set'
paired_sites <- as.data.frame(findOverlaps(targettedAsiSIGR,allAsiSIGR))
paired_sites$queryHits <- names(targettedAsiSIGR[as.integer(paired_sites$queryHits)])
paired_sites$subjectHits <- names(allAsiSIGR[as.integer(paired_sites$subjectHits)])
rownames(paired_sites) <- paired_sites$subjectHits

# CHECKED


# Check that all the sites are 8 nt long
if(!all(width(intersect(allAsiSIGR,targettedAsiSIGR)) == 8)){
   stop("Some overlaps not the expected length")
}

annotate_sites <- as.data.frame(matrix(data=0,ncol=length(annotList),nrow=length(allAsiSIGR),dimnames=list("rows"=allAsiSIGR$names,"columns"=names(annotList))))

site_features <- as.data.frame(findOverlaps(allAsiSIGR,annotList,ignore.strand=TRUE))
site_features$queryHits <- allAsiSIGR[as.integer(site_features$queryHits)]$names
site_features$subjectHits <- names(annotList)[as.integer(site_features$subjectHits)]
features_per_site <- split(site_features$subjectHits,site_features$queryHits)

for (i in names(features_per_site)){
   annotate_sites[i,features_per_site[[i]]] <- 1
}

prefered_hierarchy <- c("promoter","five_prime_utr","cds","three_prime_utr","exon","intron","intergenic")
# NOTE gene_spans excluded as redundant given other features

already_seen <- c()
sites_per_feature <- list()
targetted_sites_per_feature <- list()
for(this_level in prefered_hierarchy){
   featured_sites <- rownames(annotate_sites[annotate_sites[,this_level]==1,])
   sites_per_feature[[this_level]] <- featured_sites[! featured_sites %in% already_seen]
   targetted_sites_per_feature[[this_level]] <- rownames(paired_sites[paired_sites$subjectHits %in% sites_per_feature[[this_level]],])
   already_seen <- c(already_seen,sites_per_feature[[this_level]])
}

# CHECKED


length(unlist(sites_per_feature))
any(duplicated(unlist(sites_per_feature)))

# targetted sites by feature 
pdf(file=paste(outDir,"distribution_of_sites.pdf",sep="/"),width=12,height=6)
annotation_colours <- brewer.pal(length(prefered_hierarchy),"Set2")
names(annotation_colours) <- prefered_hierarchy
#dev.new(width=12,height=6)
par(mfrow=c(1,2))
pie(sapply(sites_per_feature,length),col=annotation_colours[names(sites_per_feature)],labels=paste(names(sites_per_feature),sapply(sites_per_feature,length),sep="\n"),main="All AsiSI sites",cex=0.6)
pie(sapply(targetted_sites_per_feature,length),col=annotation_colours[names(targetted_sites_per_feature)],labels=paste(names(targetted_sites_per_feature),sapply(targetted_sites_per_feature,length),sep="\n"),main="Active AsiSI sites",cex=0.6)
dev.off()

# CHECKED



# Check from here
active_promoter_sites_GR <- allAsiSIGR[targetted_sites_per_feature[["promoter"]]]
other_promoter_sites_GR <- allAsiSIGR[sites_per_feature[["promoter"]][!sites_per_feature[["promoter"]] %in% targetted_sites_per_feature[["promoter"]]]] 
genic_active_sites_GR <- allAsiSIGR[rownames(paired_sites[! rownames(paired_sites) %in% targetted_sites_per_feature[["intergenic"]],])]
other_genic_sites_GR <- allAsiSIGR[(!names(allAsiSIGR) %in% names(genic_active_sites_GR)) & (!names(allAsiSIGR) %in% sites_per_feature[["intergenic"]])]

intergenic_active_sites_GR <- allAsiSIGR[targetted_sites_per_feature[["intergenic"]]] 

active_promoter_sites_GR$names <- paste("promoter_active_site",seq(1,length(active_promoter_sites_GR)),sep="_")
genic_active_sites_GR$names <- paste("genic_active_site",seq(1,length(genic_active_sites_GR)),sep="_")
intergenic_active_sites_GR$names <- paste("intergenic_active_site",seq(1,length(intergenic_active_sites_GR)),sep="_")


# These sets are correct

writeGRtoBED(paste(outDir,"active_promoter_sites.bed",sep="/"),active_promoter_sites_GR)
writeGRtoBED(paste(outDir,"other_promoter_sites.bed",sep="/"), other_promoter_sites_GR)
writeGRtoBED(paste(outDir,"active_genic_sites.bed",sep="/"), genic_active_sites_GR)
writeGRtoBED(paste(outDir,"other_genic_sites.bed",sep="/"), other_genic_sites_GR)
writeGRtoBED(paste(outDir,"active_intergenic_sites.bed",sep="/"), intergenic_active_sites_GR)

sessionInfo()

