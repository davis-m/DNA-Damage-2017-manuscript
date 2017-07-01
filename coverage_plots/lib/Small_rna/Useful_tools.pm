###########################################################################

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

#############################################################################

# A perl module containing useful subroutines.

package Small_rna::Useful_tools;
use strict;
use warnings;

use Exporter qw(import);                           # Use the Exporter import function.
our @EXPORT_OK = qw(pipesystem read_desc read_option_file exclude_options);         # List of functions that can be imported on demand.

# A system call subroutine that is pipe friendly.
sub pipesystem {
   my $cline = shift;
   my $retval = system("bash", "-c", qq{set -o pipefail; $cline});
   if ($retval) {
      print STDERR "Error while executing [$cline]\n";
   }
   return $retval;
}

# Read in a sample description file. Accepts file path, a set of expected headers as an array and the column used
# to specify the names of the sample in the resultant hash.

sub read_desc {

   my %sample_info;
   my $expected_args = 3;

   die "read_desc - Expected arguments: $expected_args\n" if (scalar(@_)!=$expected_args);
   my $description_file = shift;
   my $exp_header_ref = shift;
   
   if(!(ref($exp_header_ref) eq "ARRAY")){
      die "Expected header ($exp_header_ref) must be an array reference\n";
   }

   my @exp_headers = @{$exp_header_ref};
   my $which_name = shift;
   my %expected_headers = map {$_ => 1} @exp_headers;

   die "read_desc - Chosen header ($which_name) not in specified file headers (",join(", ",@exp_headers),")\n" if (!(exists($expected_headers{$which_name})));

   die "File to be read in is missing or empty\n" if (! -s $description_file);

   open(DESC, "< $description_file") || die "Could not open the description file\n";
   my $desc_header = <DESC>;
   chomp $desc_header;

   my @headers = split "\t", $desc_header;

   my %desc_head;
   my $header_position = 0;
   foreach (@headers) {
      $desc_head{$_} = $header_position;
      $header_position++;
   }

   foreach (@exp_headers){
      die "read_desc - In table ($description_file) a header is missing: $_" if (!exists($desc_head{$_}));
   }
   my $chosen_name_col = $desc_head{$which_name};

   my $lineCount=1;
   while(<DESC>){
      chomp;
      my  @sample_line = split "\t";
      $lineCount++;

      die "All lines in the file must be the same length: $lineCount\n" if (scalar(@headers)!=scalar(@sample_line));

      my $name = $sample_line[$chosen_name_col];
      die "Multiple samples with the same name: $name\n" if(exists($sample_info{$name}));
      for (0..$#headers){
         $sample_info{$name}{$headers[$_]} = $sample_line[$_];
      }
   }
   close(DESC) || die "Could not close the description file\n";

   return(\%sample_info);
}

# Read a single line option file for passing to command line commands.
# Accepts a path to the single lined options file.
# Returns options as a character object.

sub read_option_file {
   die "Single argument needed\n" if (scalar(@_)!=1);
   my $optfile = shift;
   die "Options file $optfile does not exist\n" if (!-e $optfile);
   open(COMMOPT, "< $optfile")|| die "Could not open the command line options file: $optfile\n";
   my $optline = 0;
   my $read_options="";
   while(<COMMOPT>){
      $optline++;
      chomp;
      $read_options = $_;
      die "Expecting only a single line in the options file\n" if $optline > 1;
   }
   close(COMMOPT) || die "Could not close the command line options file: $optfile\n";
   return($read_options);
}

# Provide an array ref of options which you want to ensure are absent in the options line read by read_option_file above
sub exclude_options{
   my $expected_args = 2;
   die "exclude_options - Expected arguments: $expected_args\n" if (scalar(@_)!=$expected_args);

   my $options_line = shift;
   my $to_exclude_array_ref = shift;

   my @matches;
   my @all_options = split " {1,}",$options_line;

   for my $checking (@{$to_exclude_array_ref}){
      for my $each_option (@all_options){
         if($each_option =~ /^$checking$|^$checking=/){
            push @matches, $checking;
         }
      }
   }
   die "Options identified that must be excluded to prevent clashes: ".join(", ",@matches)."\n" if scalar(@matches) >0;
}

"A true statement required by perl";
