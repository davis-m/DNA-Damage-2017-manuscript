#!/usr/local/bin/perl

###########################################################################

# This script is part of a set of tools which aims to facilitate the
# classification of small RNAs in a streamlined and efficient manner.
# Copyright (C) 2016 2017 EMBL - European Bioinformatics Institute

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

# This script will take a set of FASTQ files and run
# Bowtie2 to align each to a genome. Alignments are performed on a
# server farm.
# SRNA_PATH must be set as an environment variable.

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

print "\n##### Attempting to submit Bowtie2 jobs for all samples.\n\n";

my $dataDir="./";
my $outDir="./";
my $logDir="./";
my $bow2opts = "";
my $tallyopts = "";
my $desc = "";
my $genome_base ="";
my $subita="$ENV{SRNA_PATH}/aligning_reads/sub.submit_Bowtie2_jobs.sh";
my $renames_ids="$ENV{SRNA_PATH}/aligning_reads/sub.rename_reads.pl";
my $minlen =0;
my $threads=1;
my $RAM = 40000;
my $ROM = 20000;
my $tallyit=0;
my $tally_by_pos=0;

if (! GetOptions(
    "dataDir=s"        => \$dataDir,      # Directory containing the read files
    "outDir=s"         => \$outDir,       # Directory for the output
    "logDir=s"         => \$logDir,       # Directory to store the bowtie logs
    "bowtie2-genome=s" => \$genome_base,  # Index for the alignments
    "bowtie2-opts=s"   => \$bow2opts,     # Options to be passed to Bowtie2
    "description=s"    => \$desc,         # A file containing 3 columns: sample name, file1 and file2
    "min-len=i"        => \$minlen,       # A length filter applied at tally step
    "tally"            => \$tallyit,      # Run tally on the FASTQ files
    "tally-opts=s"     => \$tallyopts,    # Tally options file: required if files are to be tallied
    "threads=i"        => \$threads,
    "tally_pos"        => \$tally_by_pos, # Paired reads to be paired by position by tally.
    "ROM=i"            => \$ROM,
    "RAM=i"            => \$RAM
 )){
   die "Could not parse the option information\n";
}

die "Need a description files explaining each of the samples (--description=<FILE>)\n" if (length($desc) == 0);
die "Need a set of Bowtie2 options (--bowtie2-opts=<FILE>)\n" if (length($bow2opts)==0);
die "Need a Bowtie2 index base name (--bowtie2-genome=<DIR>)\n" if (length($genome_base)==0);
die "Need a set of tally options (--tally-opts=<FILE>)\n" if ($tallyit && (length($tallyopts)==0)); 
die "The tally options file is only required if --tally is specified\n" if (($tallyit==0) && (length($tallyopts)>0)); 
die "The tally_pos flag is only accepted if --tally is specified\n" if (($tallyit==0) && ($tally_by_pos>0)); 
die "ROM requested for Bowtie2 must be greater than 0\n" if (!($ROM > 0));
die "RAM requested for Bowtie2 must be greater than 0\n" if (!($RAM > 0));

die "Require full path for the data directory (--dataDir=<DIR>): $dataDir\n" if (!($dataDir =~ /^\/.*/));
die "Require full path for the log directory (--logDir=<DIR>): $logDir\n" if (!($logDir =~ /^\/.*/));
die "Require full path for the output directory (--outDir=<DIR>): $outDir\n" if (!($outDir =~ /^\/.*/));
die "Require full path for the Bowtie2 index base name (--bowtie2-genome=<DIR>): $genome_base\n" if (!($genome_base =~ /^\/.*/));

die "Minimum length only applicable if running tally (--min-len=<INTEGER>)\n" if (($minlen > 0)&&($tallyit==0));

die "Threads must be a positive integer\n" if ($threads < 1);
die "Could not track down the Bowtie2 submission script ($subita)\n" if (!(-s $subita));
die "The description file ($desc) seems to be missing or empty\n" if (!(-s $desc));
die "Can't find the data directory!\n" if (! -d $dataDir);
die "Can't find the output directory\n" if (! -d $outDir);
die "Can't find the log file directory\n" if (! -d $logDir);
die "Can't find the tally options file ($tallyopts)\n" if ($tallyit && (! (-e $tallyopts))); 
die "Can't find the Bowtie2 options file ($bow2opts)\n" if (! (-e $bow2opts)); 

my $bowtie2_index_1 = "$genome_base.1.bt2";

die "Can't idnetify the Bowtie2 index files ($bowtie2_index_1)\n" if (! (-s $bowtie2_index_1));

die "Could not request Tally version information\n" if system("tally --version");

# HERE
sub mysystem {
   my $cline = shift;
   my $retval = system("bash", "-c", qq{set -o pipefail; $cline});
   if ($retval) {
      print STDERR "Error while executing [$cline]\n";
   }
   return $retval;
}

sub read_desc {

   my %ChIP_info;

   die "Expecting a single argument\n" if (scalar(@_)!=1);
   my $description_file = shift;

   open(DESC, "< $description_file") || die "Could not open the description file\n";
   my $desc_header = <DESC>;
   chomp $desc_header;
   die "Header of description file not as expected\n" if(! (($desc_header =~ /^name\t/) && ($desc_header =~ /\tfile1(\t|$)/) && ($desc_header =~ /\tfile2(\t|$)/)) );
   my @headers = split "\t", $desc_header;

   while(<DESC>){
      chomp;
      my  @sample_line = split "\t";

      die "All lines in the file must be the same length\n" if (scalar(@headers)!=scalar(@sample_line));

      my $name = $sample_line[0];
      die "Multiple samples with the same name: $name\n" if(exists($ChIP_info{$name}));
      for (0..$#headers){
         $ChIP_info{$name}{$headers[$_]} = $sample_line[$_]; 
      }
   }
   close(DESC) || die "Could not close the description file\n";
   
   return(\%ChIP_info);
}

my $desc_REF = &read_desc($desc);
#print Dumper \$desc_REF;

sub read_option_file {
   die "Single argument needed\n" if (scalar(@_)!=1);
   my $optfile = shift;
   die "Options file $optfile does not exist\n" if (!-e $optfile);
   open(TALLYOPT, "< $optfile")|| die "Could not open the Tally options file: $optfile\n";
   my $optline = 0;
   my $read_options="";
   while(<TALLYOPT>){
      $optline++;
      chomp;
      $read_options = $_;
      die "Expecting only a single line in the options file\n" if $optline > 1;
   }
   close(TALLYOPT) || die "Could not close the Tally options file: $optfile\n";
   return($read_options);
}


sub fastq_check {
   die "Fastq test only expecting one option\n" if (scalar(@_)!=1);
   my $potential_fastq = shift;
   open(TEST,"gzip -cdf $potential_fastq | head -n4 |") || die "Could not open $potential_fastq\n"; 
   my $id = <TEST>;
   die "Not the expected ID line\n" if (!($id =~ /^@/));
   my $seq = <TEST>;
   chomp $seq;
   die "Not the expected seq line\n" if (!($seq =~ /^[ATGCN]*$/));
   my $second_id = <TEST>;
   die "Not the expected second ID line\n" if (!($second_id =~ /^\+/));
   my $qual = <TEST>;
   chomp $qual;
   die "Not the expected quality line\n" if (!($qual =~ /^\S*$/));
   #while(<TEST>){}
   close(TEST) || die "Could not close $potential_fastq\n";
}

sub bowtie2_submission {
   my $exp_arg = 8;
   die "Expecting $exp_arg arguments\n" if (scalar(@_)!=$exp_arg);
   my($fileA, $fileB, $genomeStuff, $theoptionsfile, $sampleID, $bowOut, $threadNo, $logstore) = @_;
   
   my $theoptions = &read_option_file($theoptionsfile);  
   
   my $anyfileB = ($fileB eq "-") ? "" : "-s $fileB";
   my $thread_flag = ($threadNo==1)? "" : "-n $threadNo";

   my @time_index = localtime(time);
   my $output_log = "$logstore/$sampleID.bowtie2_run.".$time_index[3]."-".($time_index[4]+1)."-".(1900+$time_index[5]).".log";
   die "Bowtie2 output log file already exists ($output_log)\n" if (-s $output_log);
   my $output_err = "$logstore/$sampleID.bowtie2_run.".$time_index[3]."-".($time_index[4]+1)."-".(1900+$time_index[5]).".err";
   die "Bowtie2 error log file already exists ($output_err)\n" if (-s $output_err);
   my $outFile = $sampleID.".bowtie2.sam";
   die "$outFile already exists\n" if (-s $outFile);

   my $bow2call = "bsub -o $output_log -e $output_err $thread_flag -M $RAM -R 'rusage[mem=$RAM]' -R 'select[tmp>$ROM]' -R 'rusage[tmp=$ROM]' \"$subita -o $bowOut -f $fileA $anyfileB -g $genomeStuff -a '$theoptions' -t $threadNo -e $outFile\"";
   print STDERR "Submitting a call for Bowtie2:\n$bow2call\n"; 
   die "Bowtie2 failed\n" if &mysystem($bow2call);
}

sub tallying {
   die "Expecting 7 arguments\n" if (scalar(@_)!=7);
   my ($sample, $outputDir, $file1, $file2, $options, $by_pos, $minimum) = @_;
   
   my $tallyoptions = &read_option_file($options);
  
   # Check tally options for conflicting arguments
   die "Found conflicting options in Tally options file ((-l)|(pair-by-offset)|(fastax-in)|(fastqx-in))\n" if($tallyoptions=~/(\-l)|(pair\-by\-offset)|(fastax\-in)|(fastqx\-in)/);

   my @all_out;

   if ($file2 ne "-"){
      
      if ($by_pos >0){
         die "In paired end tally fastq output format, read IDs must include x followed by a read's frequency at the end of the ID\n" if (! ($tallyoptions=~/_x%C%n/));

         my $first_temp_output = "$outputDir/$sample.1.temp.fq.gz";
         my $second_temp_output = "$outputDir/$sample.2.temp.fq.gz";

         my $outputfile = "$outputDir/$sample.1.tallied.min_$minimum.fq.gz";
         my $second_outputfile = "$outputDir/$sample.2.tallied.min_$minimum.fq.gz";

         die "Output files already exist:\n$outputfile\nOR\n$second_outputfile\n" if ((-s $outputfile)||(-s $second_outputfile));
         die "Temp output files already exist:\n$first_temp_output\nOR\n$second_temp_output\n" if ((-s $first_temp_output)||(-s $second_temp_output));
         
         my $tallyCall = "tally -i $file1 -j $file2 -l $minimum $tallyoptions --pair-by-offset -o $first_temp_output -p $second_temp_output";
         print STDERR "Submitting a call to Tally for $file1 and $file2:\n".$tallyCall."\n";
         die "Could not complete the call to Tally\n" if &mysystem($tallyCall);

         # REFORMAT IDS - they must match
         my $renameCall;
         $renameCall = "$renames_ids ' 1' $first_temp_output | gzip -c1 > $outputfile";
         print STDERR "\nRenaming reads in the first output file\n$renameCall\n";
         die "Couldn't reformat the IDs of the first output file\n" if &mysystem($renameCall);
         die "Could not find the reformatted output file from Tally $outputfile\n" if (! -s $outputfile);
         
         $renameCall = "$renames_ids ' 2' $second_temp_output | gzip -c1 > $second_outputfile";
         print STDERR "\nRenaming reads in the second output file\n$renameCall\n";
         die "Couldn't reformat the IDs of the second output file\n" if &mysystem($renameCall);
         die "Could not find the reformatted output file from Tally $second_outputfile\n" if (! -s $second_outputfile);

         # DELETE TEMP OUTPUT FILES

         my $delete_first = "rm -f $first_temp_output";
         die "Could not delete the first temporary output\n" if &mysystem($delete_first);
         my $delete_second = "rm -f $second_temp_output";
         die "Could not delete the first temporary output\n" if &mysystem($delete_second);

         @all_out = ($outputfile,$second_outputfile);
      
      }else{
         die "Currently paired files can only be tallied if paired read positions are consistent in the files. If this is the case set '--tally_pos'\n";
      }

   }else{
    
      die "--tally_pos should only be set if paired files are used" if($by_pos >0);

      my $outputfile = "$outputDir/$sample.tallied.min_$minimum.fq.gz";
      die "Output file already exists ($outputfile)\n" if (-s $outputfile);
      
      my $tallyCall = "tally -i $file1 -l $minimum $tallyoptions -o $outputfile";
      print STDERR "Submitting a call to Tally for $file1:\n".$tallyCall."\n";
      die "Could not complete the call to Tally\n" if &mysystem($tallyCall);
      die "Could not find the output file from Tally $outputfile\n" if (! -s $outputfile);
   
      my $second_outputfile = "-";
      @all_out = ($outputfile,$second_outputfile);
   }

   # RETURN an array of 2 output files if no file 2 second array member should be "-";
   return(\@all_out);
}

for my $sample (keys(%{$desc_REF})){
   print STDERR "\nProcessing $sample\n";
   my $bow2OutDest = $outDir."/$sample"."_bowtie2_aligned";
   print STDERR "\nWARNING: $bow2OutDest exists but just going to carry on and pretend this didn't happen\n\n" if (-e $bow2OutDest);

   my $first_file="";
   my $second_file="";
   print STDERR "Checking FASTQ input files\n";
   if(${$desc_REF}{$sample}{"file1"}ne"-"){
      $first_file = $dataDir."/".${$desc_REF}{$sample}{"file1"};
      die "The first file does not seem to exist: $first_file\n" if (! (-s $first_file));
      &fastq_check($first_file);
   }else{
      die "Require a file1 in sample description file\n"; 
   }
   
   if(${$desc_REF}{$sample}{"file2"} ne "-"){
      $second_file = $dataDir."/".${$desc_REF}{$sample}{"file2"};
      die "The second file does not seem to exist: $second_file\n" if (! (-s $second_file));
      &fastq_check($second_file);
   }else{
      $second_file="-";
   }

   print STDERR "First file:\t$first_file\n";
   print STDERR "Second file:\t$second_file\n";

   if($tallyit){
      print STDERR "\nBeginning to remove redundancy with tally\n";
      # Paired files can be tallied if pairing is positionally conserved. Therefore expect array of results. If second file is empty "-" will be returned as second output file.

      my $tallied_file_ref = &tallying($sample, $outDir, $first_file, $second_file, $tallyopts, $tally_by_pos, $minlen);
      $first_file  = $tallied_file_ref->[0];
      $second_file = $tallied_file_ref->[1];
   }

   print STDERR "\nBeginning to align with Bowtie2\n";
   
   # If tallied second file must be '-' and first will be uniquified, otherwise both may be files.
   &bowtie2_submission($first_file, $second_file, $genome_base, $bow2opts, $sample, $bow2OutDest, $threads, $logDir);
   print STDERR "\nCompleted submission for $sample\n\n";
}

# CHECKED








