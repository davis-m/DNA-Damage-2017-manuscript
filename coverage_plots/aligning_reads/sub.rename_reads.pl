#!/usr/local/bin/perl

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

# This script reformats read IDs

use strict;
use warnings;

my $additional_id = $ARGV[0];
my $fastq_file = $ARGV[1];

open(FASTQ,"gzip -cd $fastq_file |") || die "Could not open the fastq file ($fastq_file)\n";

my $line = 0;
my $id = 0;
while(<FASTQ>){
   $line++;
   chomp;
   my $current = $_;
   if(($line%4)==1){
      if(/@\w+_x(\d+)$/){
         $id++;
         print "@"."I$id"."_x$1"."$additional_id\n";
      }else{
         die "Could not parse the ID\n";
      }
   }else{
      print "$current\n";
   }
}
close(FASTQ) || die "Could not close the fastq file\n";

