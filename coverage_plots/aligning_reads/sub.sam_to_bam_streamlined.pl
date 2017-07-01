#!/usr/bin/env perl

###########################################################################
#
# This script is part of a lncRNA/small RNA annotation pipeline.
# Copyright (C) 2011 2012 2013 2014 EMBL - European Bioinformatics Institute 
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
###########################################################################

# 161112: sam_to_bam.pl
# Takes a sam file and uses samtools to convert this to a sorted and indexed bam file.
# Must be provided with a basename to organise the bam files.

use strict;
use warnings;
use Getopt::Long;

my $samFile = "";
my $fileBase = "";
my $outDir = "./";


if(!GetOptions( 
   "out=s"           => \$outDir,
   "sam=s"           => \$samFile,
   "base=s"          => \$fileBase
)){
   die "Option parsing has imploded spectacularly!";
}

die "Require a basename for output: --base\n" if (length ($fileBase) == 0);
die "Require a sam file as input: --sam\n" if (length ($samFile) == 0);
die "Output directory doesn't exist\n" if (! -d $outDir);
die "Sam file doesn't exist\n" if (! -s $samFile);

die "Could not check Samtools version\n" if (system("samtools --version 1>&2"));

print STDERR "# Converting $samFile to sorted and indexed BAM format #\n";
print STDERR "Placing output here: $outDir\n";

die "Output directory must not be empty\n" if (length($outDir)==0);
my $bamFile = "$outDir/$fileBase.bam";
#my $bamSort = "$outDir/$fileBase.sort";
my $indexIn = "$outDir/$fileBase.sort.bam";

my $viewCall = "gzip -cdf $samFile | samtools view -bS -o $bamFile -";
my $sortCall = "samtools sort $bamFile -T $outDir/$fileBase -o $indexIn";
my $indexCall = "samtools index $indexIn";

print STDERR "\nConverting SAM to BAM format\n";
print STDERR "$viewCall\n";
die "Failed to convert $samFile to bam format\n" unless system($viewCall)==0;

print STDERR "\nSorting BAM file\n";
print STDERR "$sortCall\n";
die "Failed to sort $bamFile\n" unless system($sortCall)==0;

print STDERR "\nIndexing BAM file\n";
print STDERR "$indexCall\n";
die "Failed to index $indexIn\n" unless system($indexCall)==0;

print STDERR "Cleaning up files\n";
my $BAM_delete = "rm -f $bamFile";
print STDERR "$BAM_delete\n";
die "Failed to delete BAM file\n" if system($BAM_delete);

print STDERR "\n# We have reached the end of the file conversion #\n";
