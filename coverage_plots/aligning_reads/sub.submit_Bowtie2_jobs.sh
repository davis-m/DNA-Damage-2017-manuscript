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

# This script organises alignments with Bowtie2 in tmp.

set -e
shopt -s nullglob
thetmp=/tmp/bowtie.mpd.$$


function rmit {
   errcode=$?
   rm -rf $thetmp
   exit $errcode
}

trap rmit SIGTERM  EXIT

dest=no-such-dest
file1=
file2=
threads=none
genomeBase=
options=none
samfile=
cmmd=

while getopts :f:s:o:t:a:g:e:h opt
do
    case "$opt" in
    o)
      dest=$OPTARG
      ;;
    f)
      file1=$OPTARG
      ;;
    s)
      file2=$OPTARG
    ;;
    t)
      threads=$OPTARG
    ;;
    a)
      options=$OPTARG
    ;;
    g)
      genomeBase=$OPTARG
    ;;
    e)
      samfile=$OPTARG
    ;;
    h)
      cat <<EOU
-f file of first end sequences
-s file of second end sequences
-o output directory
-t threads number
-g genome base name
-a alignment options
-e sam file output
EOU
      exit
      ;;
    :) echo "Flag $opt needs argument"
        exit 1;;
    ?) echo "Flag $opt unknown"
        exit 1;;
   esac
done

mkdir -p $thetmp
cd $thetmp

# Does the current Bowtie2 output file exist - if so exit

if [ -f $dest/$samfile.gz ]
   then
   echo "$samfile.gz exists."
   false   
fi

# Organise output file name

[[ -z "$file2" ]] && cmmd="bowtie2 $options -p $threads -x $genomeBase -U $file1 | gzip -c1 > $samfile.gz" || cmmd="bowtie2 $options -p $threads -x $genomeBase -1 $file1 -2 $file2 | gzip -c1 > $samfile.gz"

echo $cmmd
eval $cmmd

mkdir -p $dest


for f in *; do
   cp -r $f $dest
done


# CHECKED
