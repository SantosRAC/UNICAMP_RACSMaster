#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

my $version='0.1';
my $license='';
my $help='';
my $inputMatrix='';
my $outputMatrix='';
my @taxa=();
my %family2countsPerTaxon;
my %filteredFamily2countsPerTaxon;
my @matrixHeader='';

#Get user input
GetOptions(
    'license|l'          => \$license,
    'help|h|?'           => \$help,
    'in_matrix|in=s'     => \$inputMatrix,
    'out_matrix|out=s'   => \$outputMatrix
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
 print STDERR "FATAL: you must provide an input matrix with species IDs (rows) and normalized protein counts (columns), with columns separated by tabulation (TAB).\n";
 &usage;
 exit(1);
}
if(!$outputMatrix){
 print STDERR "FATAL: You must provide an output file name.\n";
 &usage;
 exit(1);
}
if(-e $outputMatrix){
 print STDERR "FATAL: The output file you informed already exists. Please, change the name or delete the older file.\n";
 &usage;
 exit(1);
}

open(INMATRIX,$inputMatrix);

my $lineNumb=0;

while(<INMATRIX>){
 chomp;
 $lineNumb++;
 if (/^Organism/) {
  @matrixHeader=split(/\t/,$_);
  #print join("\t",@matrixHeader)."\n";
 }
 next if (/^Organism/);
 my @fields=split(/\t/,$_);
 push(@taxa,$fields[0]);
 foreach my $i (1 .. $#fields) {
  push(@{$family2countsPerTaxon{$matrixHeader[$i]}},$fields[$i])
  #print "$i - $fields[$i]";
 }
}

close(INMATRIX);

foreach my $family (keys %family2countsPerTaxon){
 if (keys %{{ map {$_, 1} @{$family2countsPerTaxon{$family}} }} == 1) {
  #print "ALL EQUAL: ".join("\t",@{$family2countsPerTaxon{$family}})."\n";
 } else {
  @{$filteredFamily2countsPerTaxon{$family}}=@{$family2countsPerTaxon{$family}};
  #print "NOT EQUAL: ".join("\t",@{$family2countsPerTaxon{$family}})."\n";
 }
}

open(OUTMATRIX,">>",$outputMatrix);
print OUTMATRIX "Organism\t".join("\t",keys %filteredFamily2countsPerTaxon)."\n";
close(OUTMATRIX);

foreach my $i (0 .. $#taxa) {
 open(OUTMATRIX,">>",$outputMatrix);
 print OUTMATRIX $taxa[$i];
 foreach my $f (keys %filteredFamily2countsPerTaxon) {
  print OUTMATRIX "\t".$filteredFamily2countsPerTaxon{$f}[$i];
 }
 print OUTMATRIX "\n";
 close(OUTMATRIX);
}

sub usage{
    print STDERR "$0 version $version, Copyright (C) 2016 Renato Augusto Correa dos Santos\n";
    print STDERR "$0 comes with ABSOLUTELY NO WARRANTY; for details type `$0 -l'.\n";
    print STDERR "This is free software, and you are welcome to redistribute it under certain conditions;\n";
    print STDERR "type `$0 -l' for details.\n";
    print STDERR <<EOF;
NAME
    $0 Filters out columns with with no variance.

USAGE
    $0 

OPTIONS
    --in_matrix       -in   Input matrix with normalized counts                       REQUIRED
    --out_matrix      -out  Output matrix with normalized counts (columns filtered).  REQUIRED
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
