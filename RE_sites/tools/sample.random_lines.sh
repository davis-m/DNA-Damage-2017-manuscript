#!/usr/local/bin/bash

#######################################################################

# This script is part of a set of tools which aims to facilitate the
# classification of small RNAs in a streamlined and efficient manner.
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

#######################################################################

# This script will take a file, a number of lines to randomly sample and a
# number of samples to take. It will then produce random subsamples.

set -o pipefail
set -e
shopt -s nullglob  

while getopts :i:t:s:n:b:h opt  # :d:o:h -> kind of options (eg. argument required?). Leading : prevents printing of automated error when unrecognised optio
do
    case $opt in
    b)
      base=$OPTARG            # Base name for the samples to be produced
      ;;
    s)
      sample=$OPTARG          # Number of samples to take
      ;;
    t)
      suffix=$OPTARG          # File type to be appended to the output file name
      ;;
    n)
      number=$OPTARG          # Number of lines per sample
      ;;
    i)
      input=$OPTARG           # Input file to sample
      ;;
    h)
      cat <<EOU
-b base name for the output files
-s number of samples to take from the file
-n number of lines to sample
-t suffix to be tagged to the ouptut file names
-i input file to sample from
EOU
      exit
      ;;
    :) echo "Flag $opt needs argument"
        exit 1;;
    ?) echo "Flag $opt unknown"
        exit 1;;
   esac
done

if [[ -z $base || -z $sample || -z $suffix || -z $number || -z $input ]]
then
   echo "-b, -s, -t, -n and -i require options"
   exit 1
fi

if ! [[ $sample =~ ^[0-9]+$ && $number =~ ^[0-9]+$ ]]
   then
   echo "-s and -n take a positive integer"
   exit 1
fi


if [[ ! -s $input ]]
   then
   echo "-i must be a file of non-zero size"
   exit 1
fi

word_lines=$(gzip -cdf $input | wc -l)
if [[ $number -gt $word_lines ]]
then
   echo "Sampling more lines than there are in the file"
   exit 1
fi

for i in `seq 1 $sample`;
do
   outfile=$base.$i.$suffix;
   eval "gzip -cdf $input | shuf -n $number -o $outfile"
done
   

