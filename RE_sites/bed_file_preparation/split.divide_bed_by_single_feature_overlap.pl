#!/usr/local/bin/perl

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

use strict;
use warnings;
use Getopt::Long;

# This script will take 2 bed files of features and separate the first BED file
# into elements that overlap the features in the second and those that don't. A
# base name for output is also required. Features in each output will also be
# merged to remove overlaps Merging is non-strand specific and collapses to 3
# columns.

my $baseBED   = "";
my $featuresBED = "";
my $basename = "";
my $mergeit = 0;

if (! GetOptions(
   "basebed=s"       =>    \$baseBED,
   "featuresbed=s"   =>    \$featuresBED,
   "outname=s"    =>    \$basename,
   "merge"        =>    \$mergeit            # Merging ouput is optional
)){
   die "Could not parse the command line options\n";
}

die "Require a base BED file from which elements will be seperated (--basebed)\n" if length($baseBED)==0;
die "Require a BED file containing elements on whiich the first BED file will be divided (--featuresbed)\n" if length($featuresBED)==0;
die "Require a base name for all the output (--outname)\n" if length($basename)==0;

die "Base BED file does not exist or is empty\n" if (! -s $featuresBED);
die "Feature BED file does not exist or is empty\n" if (! -s $featuresBED);

print STDERR "\nSeperating Features that overlap the restriction sites from non-overlapping sites\n\n";

my $overlapped_feature_file = "$basename.overlapping.bed";
my $not_overlapped_feature_file = "$basename.not_overlapping.bed";

my $sorted_overlapped_feature_file = "$basename.overlapping.sorted.bed";
my $sorted_not_overlapped_feature_file = "$basename.not_overlapping.sorted.bed";

my $merged_overlapped_feature_file = "$basename.overlapping.merged.bed";
my $merged_not_overlapped_feature_file = "$basename.not_overlapping.merged.bed";

my $overlapped_features = "bedtools intersect -wa -a $baseBED -b $featuresBED > $overlapped_feature_file";
my $non_overlapped_features = "bedtools intersect -wa -v -a $baseBED -b $featuresBED > $not_overlapped_feature_file";

my $sort_overlapped_features = "sort -k1,1 -k2,2n $overlapped_feature_file > $sorted_overlapped_feature_file";
my $sort_not_overlapped_features = "sort -k1,1 -k2,2n $not_overlapped_feature_file > $sorted_not_overlapped_feature_file";

my $merge_overlapped_features = "bedtools merge -i $sorted_overlapped_feature_file > $merged_overlapped_feature_file";
my $merge_non_overlapped_features = "bedtools merge -i $sorted_not_overlapped_feature_file > $merged_not_overlapped_feature_file";

print STDERR "Finding overlapping regions:\n$overlapped_features\n\nFinding non-overlapping regions:\n$non_overlapped_features\n";
die "Could not find regions overlapping with the features of interest\n" if (system($overlapped_features));
die "Could not find regions not overlapping with the features of interest\n" if (system($non_overlapped_features));

if($mergeit){
   print STDERR "Merging overlapping regions:\n$sort_overlapped_features\n$merge_overlapped_features\n\nMerging non-overlapping regions:\n$sort_not_overlapped_features\n$merge_non_overlapped_features\n";
   die "Could not sort regions overlapping with the features of interest\n" if (system($sort_overlapped_features));
   die "Could not merge regions overlapping with the features of interest\n" if (system($merge_overlapped_features));
   die "Could not sort regions not overlapping with the features of interest\n" if (system($sort_not_overlapped_features));
   die "Could not merge regions not overlapping with the features of interest\n" if (system($merge_non_overlapped_features));
}

print STDERR "\nComplete\n";
