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

# This script coordinates the conversion of sam files to bam files

use strict;
use warnings;

use Getopt::Long;

my $samDir = "";
my $outDir = "";
my $conversionScript = "$ENV{SRNA_PATH}/aligning_reads/sub.sam_to_bam_streamlined.pl";

if(! GetOptions(
   "sams=s" => \$samDir,
   "out=s"  => \$outDir
)){
   die "Can't process the command line\n";
}

die "Need a directory of sam files\n" if (length($samDir)==0);
die "Need a directory for the output\n" if (length($outDir)==0);

die "Can't find the sam directory\n" if (! -d $samDir);
die "Can't find the out directory\n" if (! -d $outDir);

my @all_sams = <$samDir/*sam*>;

#{
#   local $, = "\t";
#   print @all_sams;
#   print "\n";
#}

for my $a_sam (@all_sams){
   my $sample_id = "";

   #print STDERR "Converting $a_sam to BAM\n";

   if($a_sam=~/^.*\/([\w\-\.]+)\.sam.*$/){
      $sample_id = $1;
   }else{
      die "Could not find the sample ID\n";
   }

   my $convertCall = "$conversionScript --out=$outDir --sam=$a_sam --base=$sample_id";

   print STDERR "\n$sample_id:\n$convertCall\n";
   die "Could not run the SAM conversion call\n" if (system($convertCall));
}

print STDERR "Conversions complete!\n";
