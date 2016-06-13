#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

my $version='0.1';
my $license='';
my $help='';
my $inputMatrix='';
my $normOutput='';
my $protSizeMap='';
my %spe2protSize;
my $matrixHeader='';

#Get user input
GetOptions(
    'license|l'          => \$license,
    'help|h|?'           => \$help,
    'in_matrix|in=s'     => \$inputMatrix,
    'proteome_size|ps=s' => \$protSizeMap,
    'out_matrix|out=s'   => \$normOutput
);

#Check user input
if($help){
 &usage;
 exit(1);
}
if($license){
 &license;
 exit(1);
}
if(!-s $inputMatrix){
 print STDERR "FATAL: you must provide an input matrix with species IDs (rows) and protein counts (columns), with columns separated by tabulation (TAB).\n";
 &usage;
 exit(1);
}
if(!-s $protSizeMap){
 print STDERR "FATAL: you must provide an input mapping file with species IDs (column 1) and proteome size (column 2) separated by tabulation (TAB). Notice that the species identifiers in input files must be the the same.\n";
 &usage;
 exit(1);
}
if(!$normOutput){
 print STDERR "FATAL: You must provide an output file name.\n";
 &usage;
 exit(1);
}
if(-e $normOutput){
 print STDERR "FATAL: The output file you informed already exists. Please, change the name or delete the older file.\n";
 &usage;
 exit(1);
}

getMapInfo($protSizeMap);

my $countLinesInMatrix = `wc -l < $inputMatrix`;
chomp($countLinesInMatrix);
#print "The input matrix has $countLinesInMatrix lines\n";

my $countLinesInProtSizeMap = `wc -l < $protSizeMap`;
chomp($countLinesInProtSizeMap);
#print "The input proteome size mapping file has $countLinesInProtSizeMap lines\n";

if (!($countLinesInMatrix == ($countLinesInProtSizeMap + 1))) {
 die "FATAL: Input character matrix is expected to have the same number of lines in mapping file + header.
Your matrix has $countLinesInMatrix lines, and the mapping file has $countLinesInProtSizeMap.\n";
}

open(INMATRIX,$inputMatrix);

while(<INMATRIX>){
 chomp;
 if (/^Organism/) {
  $matrixHeader=$_;
  chomp($matrixHeader);
  open(OUTMATRIX,">>",$normOutput);
  print OUTMATRIX "$matrixHeader\n";
  close(OUTMATRIX);
 }
 next if (/^Organism/);
 my @countsForSpecies=split(/\t/,$_);
 my $normalizedCount;
 my $speciesID = shift(@countsForSpecies);
 open(OUTMATRIX,">>",$normOutput);
 print OUTMATRIX "$speciesID";
 foreach my $count (@countsForSpecies) {
 $normalizedCount = ($count/$spe2protSize{$speciesID})*1000;
 print OUTMATRIX "\t".sprintf("%.3f", $normalizedCount);
 }
 print OUTMATRIX "\n";
 close(OUTMATRIX);
}

close(INMATRIX);


# Get information from the mapping file and fill hash species2proteomeSize.
sub getMapInfo {

 my $map = shift;
 open(MAPFILE,$map);

 while(<MAPFILE>){
  chomp;
  my ($speID,$proteomeSize) = split(/\t/,$_);
  if(!$spe2protSize{$speID}){
   $spe2protSize{$speID}=$proteomeSize;
   #print "$spe2protSize{$speID}\t$speID\n";
  } else {
   die "FATAL: There is a repeated name for species (namely, $speID) in the mapping file.\n";
  }
 }

 close(MAPFILE);

}


sub usage{
    print STDERR "$0 version $version, Copyright (C) 2016 Renato Augusto Correa dos Santos\n";
    print STDERR "$0 comes with ABSOLUTELY NO WARRANTY; for details type `$0 -l'.\n";
    print STDERR "This is free software, and you are welcome to redistribute it under certain conditions;\n";
    print STDERR "type `$0 -l' for details.\n";
    print STDERR <<EOF;
NAME
    $0 Converts protein count data in a matrix to normalized (pseudo-continuous) counts, based on proteome size.

USAGE
    $0 

OPTIONS
    --in_matrix       -in   Matrix with non-normalized protein counts per species              REQUIRED
    --proteome_size   -ps   Mapping file with two columns with species and protein counts.     REQUIRED
    --out_matrix      -out  Output matrix with normalized counts.                              REQUIRED
    --help            -h    This help.
    --license         -l    License.

EOF
}

sub license {
print STDERR <<EOF;

Copyright (C) 2016 Renato Augusto Correa dos Santos
e-mail: renatoacsantos\@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
EOF
exit;
}
