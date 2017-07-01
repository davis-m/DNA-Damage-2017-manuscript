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

# This script will search for specified motifs in a fasta file.

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Perl;
use Data::Dumper;

my $fasta = "";
my $motif = "";

if(!(GetOptions(
   "fasta=s" => \$fasta,    
   "enzyme-cut=s" => \$motif
))){
   die "Could not get the command line options";
}

die "Need an input FASTA file: --fasta\n" if length($fasta)==0;
die "Input fasta file does not exist\n" if (! -s $fasta);
die "Need an input motif for an enzyme: --enzyme-cut\n" if length($motif)==0;

my $fasta_in = Bio::SeqIO->new(-file => "$fasta" , -format=>"Fasta", -alphabet => "dna");
if($motif=~/[^AGTC]/){
   die "Motif must only contain ATGC\n"
}

my $rev_motif = reverse_complement_as_string($motif);
print STDERR "Comparing motif to its reverse complement ($motif vs. $rev_motif)\n";


print STDERR "Beginning motif search\n";

my $motif_length=length($motif);
my $seq_count =0;
while ((my $seqobj = $fasta_in->next_seq())) {
   $seq_count++;
   my $sequence = $seqobj->seq();
   my $chromosome = $seqobj->display_id();
   print STDERR "Sequence number: $seq_count - chr $chromosome\n";
   
   #print "Sequence ",$seq->id, " first 10 bases ",
   #$seq->subseq(1,10), "\n";

   #print $sequence."\n";

   # Pattern match returns 0 based coordinates - returning BED 0-based coordinates

   if($rev_motif eq $motif){
      while($sequence=~/(?=$motif)/g){
         # NOTE: BED definition does not accept * in strand so all sites defined as undefined if palendrome
         print "$chromosome\t".($-[0])."\t".($-[0]+$motif_length)."\t$motif\t0\t.\n";
      }
   }else{
      while($sequence=~/(?=$motif)/g){
         print "$chromosome\t".($-[0])."\t".($-[0]+$motif_length)."\t$motif\t0\t+\n";
      }
      while($sequence=~/(?=$rev_motif)/g){
         print "$chromosome\t".($-[0])."\t".($-[0]+$motif_length)."\t$motif\t0\t-\n";
      }
   }
}
