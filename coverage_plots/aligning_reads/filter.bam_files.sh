#!/usr/local/bin/bash

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

# This script will apply a series of filters to a directory of bam files, including the removal of mapped duplicates and low quality alignments.
# NOTE: all read manipulation is performed in tmp/

set -o pipefail
set -e
shopt -s nullglob          # Alternative to set environment variables (try shopt with no options to list variables)
                           # nullglob: when wildcard is not matched empty string returned (not pattern)

# Check that path to Picard .jar file is set
echo $PICARD_JAR
: "${PICARD_JAR:?env variable PICARD_JAR must be set as path to picard .jar file}"
                           
thetmp=/tmp/filtering.mpd.$$

# Check versions of tools
echo -e "\nJava version:"
eval "java -version"
echo -e "\nPicard MarkDuplicates version:"
eval "java -jar $PICARD_JAR MarkDuplicates --version" || true
echo -e "\nSamtools version:"
eval "samtools --version"

### Function to perform tmp cleanup if necesary

function rmit {
   #cp $thetmp/*redundant.fastq.gz $dest ### SHOULD INSERT || true clause?
   code=$?                 # Exit code of last executed command
   cd $dest
   rm -rf $thetmp
   exit $code
}

### Trap in case of epic fail (cleanup)

trap rmit SIGTERM EXIT     # trap will be initiated if a SIGTERM (kill) or EXIT (process) received

#export PATH=/nfs/research2/enright/lncRNAs/DATA/bin/:${PATH}

dest=no-such-dest
bamdir=no-base-supplied
qual=-1
rem=0
r_mem=no-memory-supplied
throw=-1
keep=-1

### Pass options

### Add options to remove or retain reads based on SAM bit flags k-eep and t-hrow

while getopts :o:k:t:b:q:m:rh opt  # :d:o:h -> kind of options (eg. argument required?). Leading : prevents printing of automated error when unrecognised option (? is set). opt -> loaded with each option letter in turn as processed
do
    case $opt in
    o)
      dest=$OPTARG              # This should be the complete path to the destination directory.
      ;;
    b)
      bamdir=$OPTARG            # This should be the directory containing the individual folders for each aligned sample.
      ;;
    q)
      qual=$OPTARG              # The quality cutoff score
      ;;
    m)
      r_mem=$OPTARG             # Memory for removal of duplicates
      ;;
    t)
      throw=$OPTARG             # Flags to be passed to -F of samtools for flag filtering
      ;;
    k)
      keep=$OPTARG              # Flags to be passed to -f of samtools for flag filtering
      ;;
    r)
      rem=1                     # Remove duplicate alignments with Picard
      ;;
    h)
      cat <<EOU
-b bam directory containing individual bam files to be filtered
-o output directory
-q quality score cutoff (must be greater than 0)
-r flag to specify whether picard tools should be used to remove duplicate alignments
-m the memory to be passed to java (-Xmx) for picard tools when removing duplicate alignments
-t bit-flag integer value to be passed to samtools view -f option
-k bit-flag integer value to be passed to samtools view -k option
EOU
      exit
      ;;
    :) echo "Flag $opt needs argument"
        exit 1;;
    ?) echo "Flag $opt unknown"
        exit 1;;
   esac
done

if ! [[ $qual =~ ^-?[0-9]+$ && $throw =~ ^-?[0-9]+$ && $keep =~ ^-?[0-9]+$ ]]
   then
   echo "-q, -k and -t take an integer"
   exit
fi


# At least one of -k, -t, -q or -r must be specified to apply a filter
if [[ $rem -lt 1 ]] && [[ $qual -lt 0 ]] && [[ $throw -lt 0 ]] && [[ $keep -lt 0 ]]
   then
   echo "-r, -k, -t or -q must be specified"
   exit
fi


# Require a memory specification for Java if r specified
if [[ $rem -eq 1 ]] && [[ $r_mem == "no-memory-supplied" ]]
   then
   echo "if -r specified then -m is required"
   exit
fi

# If $qual less than 0 exit -> require a quality score cutoff
if [[ ! -d $dest ]]
   then
   echo "-o must be a directory"
   exit
fi 

if [[ ! -d $bamdir ]]
   then
   echo "-b must be a directory"
   exit
fi 


# Add pattern match to check complete paths provided
if [[ ! $bamdir == /* ]]
   then
      echo "Require a complete path to files"
      exit
fi

if [[ ! $dest == /* ]]
   then
      echo "Require a complete path to destination directory"
      exit
fi
#exit

mkdir -p $thetmp
cd $thetmp

# HERE

for f in $bamdir/*.bam; do
   g=${f%%.bam}
   h=${g##*/}
#   echo $h
   
   # Sorting directory for the sample
   # Make a directory for the sample and move into it
   sample_tmp=$thetmp/sample
   
   echo -e "\n### Processing file $f\nMoving to directory $sample_tmp"
   mkdir -p $sample_tmp
   cd $sample_tmp

   #which samtools
   #cmmd="samtools view -bh -q50 -F12 $f | head > $h.filtered_hits.bam"
   # Alter so that q, t or k could be specified

   if [[ $qual -ge 0 || $throw -ge 0 || $keep -ge 0 ]]
      then

      # Set -q -f and -F depending upon whether option is supplied above otherwise will be empty
      # Add option variables (containing option or empty) to the samtools view command

      throwopt=""
      keepopt=""
      qualopt=""

      if [[ $qual -ge 0 ]]; then
         qualopt=" -q$qual " 
      fi
      
      if [[ $throw -ge 0 ]]; then
         throwopt=" -F$throw " 
      fi
      
      if [[ $keep -ge 0 ]]; then
         keepopt=" -f$keep " 
      fi
      
      # Change name of bam files here and below to .quality_or_flag_filtered.bam

      echo -e "\nBeginning samtools view filter:\nOptions:$qualopt$throwopt$keepopt\n# Commands:"

      # Filter based on quality first so that loci aren't removed because low quality alignment with high read quality doesn't get selected in the place of duplicates and then filtered
      cmmd="samtools view -bh$qualopt $throwopt $keepopt $f > $h.filtered_intermed.bam"
      
      #cmmd="samtools view -bh$qulopt$throwopt$keepopt $h.resorted_for_quality.bam > $h.quality_or_flag_filtered.bam"

      resortcmmd="samtools sort -o $h.quality_or_flag_filtered.bam -O bam -T $h.quality_sort $h.filtered_intermed.bam"
      
      #resortcmmd="samtools sort -o $h.resorted_for_quality.bam -O bam -T $h.quality_sort $f"
      
      echo -e "\n$cmmd"
      eval $cmmd
      echo -e "\n$resortcmmd"
      eval $resortcmmd
      f="$h.quality_or_flag_filtered.bam"
   fi
   
   # CHECKED
   
   if [[ $rem -eq 1 ]]
      then
      
      # Make a temp file for Picard overspill
      tmp_spill=$sample_tmp/spill
      mkdir -p $tmp_spill

      # NOTE: Reads must be FR format.

      outfile=""
      if [[ $qual -ge 0 || $throw -ge 0 || $keep -ge 0 ]]; then
         outfile="$h.quality_or_flag_filtered.rmdups.bam"
      else
         outfile="$h.rmdups.bam"
      fi

      echo -e "\nProceeding with duplicate removal\nCommands:"
      sortcmmd="samtools sort -n -o $h.resorted.bam -O bam -T $h.name_sort $f"
      fixcmmd="samtools fixmate $h.resorted.bam $h.mate_fixed.bam"
      
      # Samtools rmdup only removing one end of a duplicate alignment without correcting flags. Seems to be a known bug.

      #rmcmmd="samtools rmdup $h.mate_fixed.bam $h.rmdup.bam"
      resortcmmd="samtools sort -o $h.resorted_fixed.bam -O bam -T $h.coord_sort $h.mate_fixed.bam"
      # NOTE: Default picard duplicate selected on the basis of read quality.
      rmcmmd="java -Xmx$r_mem -jar $PICARD_JAR MarkDuplicates I=$h.resorted_fixed.bam O=$outfile M=$h.rmdup_stats.txt REMOVE_DUPLICATES=true READ_NAME_REGEX=null TMP_DIR=$tmp_spill"
   
      echo -e "\n$sortcmmd"
      eval $sortcmmd
      echo -e "\n$fixcmmd"
      eval $fixcmmd
      echo -e "\n$resortcmmd"
      eval $resortcmmd
      echo -e "\n$rmcmmd"
      eval $rmcmmd
      f=$outfile
   
      # Remove the tmp overspill
   fi

   
   # Also create an index for the new filtered file
   indexcmmd="samtools index $f" 

   echo -e "\nIndexing bam output:\n$indexcmmd"
   eval $indexcmmd

   # Collect and move all of the files needed 
   #touch $h.filtered_hits.bam
   echo -e "\nCopying files to $dest"
   
   for i in *.quality_or_flag_filtered.bam; do
      mv $i $dest/
   done
   
   for i in *.rmdups.bam; do
      mv $i $dest/
   done
   
   for i in *.rmdup_stats.txt; do
      mv $i $dest/
   done

   for i in *.bai; do
      mv $i $dest/
   done
   
   # Move out of the sample tmp directory and delete it
   echo -e "Leaving $sample_tmp and cleaning up\n"
   cd $thetmp
   rm -rf $sample_tmp
done


# CHECKED
