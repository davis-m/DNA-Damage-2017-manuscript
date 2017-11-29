###########################################################################
#
# This script is part of a set of tools which aims to facilitate the
# classification of small RNAs in a streamlined and efficient manner.
# Copyright (C) 2017 EMBL - European Bioinformatics Institute
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

# This script takes the profiles determined for active AsiSI sites and other
# sites It will determine the difference between the profiles at these sites in
# damage and no damage conditions. The difference is determined according to
# the mean values of the two central matrix windows that flank the AsiSI site.
# Ultimately a Mann-Whitney test is performed to determine if the distributions
# of the differences in the profiles at the two sets of sites are significantly
# different. This script is specific to the publication analysis and will need
# editing to be appicable to wider sample sets.

args <- R.utils::commandArgs(asValue=TRUE)
if(! is.null(args$output)){
   outDir <- as.character(args$out)
}else{
   stop("Require an output directory")
}

if(! is.null(args$input)){
   inDir <- as.character(args$input)
}else{
   stop("Require an input directory")
}

if(! is.null(args$stats)){
   statfiles <- as.character(args$stats)
}else{
   stop("Require a table containing the stats files and their corresponding order file")
}

site_set_stats <- read.table(file=statfiles,sep="\t",as.is=TRUE,row.names=1)

if(!all(c("active_damage","active_no_damage","other_damage","other_no_damage") %in% rownames(site_set_stats))){
   stop("Require specific rownames in the table containing the stats files")
}

site_table <- list()
flanking_window_means <- list()
for(i in rownames(site_set_stats)){
   site_table[[i]] <- read.table(file=paste(inDir,"/",site_set_stats[i,1],sep=""),sep="",as.is=TRUE)
   active_order_table <- read.table(file=paste(inDir,"/",site_set_stats[i,2],sep=""),sep="\t",as.is=TRUE)
   if(!(nrow(active_order_table)==nrow(site_table[[i]]))){
      stop("Ordered regions and heatmap matrix not the same length")
   }
   rownames(site_table[[i]]) <- active_order_table[,4]
   if(ncol(site_table[[i]])!=4000){
      stop("Unexpected column number in computeMatrix")
   }
   # Take the mean value of the regions flanking the centre of the AsiSI site
   flanking_window_means[[i]] <- apply(site_table[[i]][,c(2000,2001)],1,mean)
   #flanking_window_means[[i]] <- apply(site_table[[i]][,c(1991,2010)],1,mean)
}

differences <- list()
differences[["active"]] <- flanking_window_means[["active_damage"]] - flanking_window_means[["active_no_damage"]][names(flanking_window_means[["active_damage"]])]
differences[["other"]] <- flanking_window_means[["other_damage"]] - flanking_window_means[["other_no_damage"]][names(flanking_window_means[["other_damage"]])]

pdf(file=paste(outDir,"/promoter.AsiSI_sites.profile_differences.pdf",sep=""))
boxplot(differences,ylab="Difference between damage and no damage profiles",xlab="Promoter AsiSI site sets",col="lightgrey",main="Differences between the damage and no damage profiles\nfor 2 sets of promoter associated AsiSI sites",pch=20)

# Mann–Whitney U test
test_result <- wilcox.test(differences[["active"]],differences[["other"]])
legend("topright",cex=0.7,legend=paste("Mann-Whitney test p-value:",signif(test_result$p.value,digits=3)))
dev.off()
