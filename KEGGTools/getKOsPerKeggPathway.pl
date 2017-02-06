#!/usr/bin/perl

use warnings;
use strict;
use LWP::Simple;
use Getopt::Long;

my $version='0.1';
my $license='';
my $keggPathwaysFile='';
my $gmtFile='';
my $help='';
my %path2KO;

GetOptions(
  'help|h|?'    => \$help,
  'license|l'   => \$license,
  'infile|i=s'  => \$keggPathwaysFile,
  'outfile|o=s' => \$gmtFile,
);

if(!-s $keggPathwaysFile) {
  print "The input file with Kegg pathways is required.\n";
  &usage();
  exit(1);
}

if(-s $gmtFile) {
  print "The output file already exists.\n";
  &usage();
  exit(1);
}

if(!$gmtFile) {
  print "The output file name is required.\n";
  &usage();
  exit(1);
}

if($help) {
  &usage();
  exit(0);
}

if($license) {
 &license();
 exit(0);
}

open(KEGGPATHWAYS,$keggPathwaysFile);

while(<KEGGPATHWAYS>) {
 chomp;
 my ($keggPathwayID,$keggPathwayDesc)=split(/\t/,$_);
 $keggPathwayID =~ s/path://g;
 @{$path2KO{$keggPathwayID}}=();
 my $kegg_url = "http://rest.kegg.jp/link/ko/$keggPathwayID";
 my $content = get($kegg_url);
 if($content =~ /path:(map\d+)	ko:(K\d+)/) {
  my @contentLines = split(/\n/,$content);
  foreach my $line (@contentLines){
   my (undef,$ko)=split(/\t/,$line);
   $ko =~ s/ko://;
   if($ko ~~ @{$path2KO{$keggPathwayID}}){
    die "Duplicated KO for a pathway ($keggPathwayID): $ko\n";
   }
   push(@{$path2KO{$keggPathwayID}},$ko);
  }
  open(GMT,">>",$gmtFile);
  print GMT "$keggPathwayID\t$keggPathwayDesc\t".join("\t",@{$path2KO{$keggPathwayID}})."\n";
  close(GMT);
 }
 sleep 5;
}

close(KEGGPATHWAYS);

sub usage {
    print STDERR "$0 version $version, Renato Augusto Correa dos Santos\n";
    print STDERR <<EOF;

NAME
    $0 takes a list of Kegg Pathways resulting from using Kegg List API, and returns a GMT file (Pathway and corresponding KOs)

USAGE
    $0 --infile kegg.list.txt --outfile OUTFILE.gmt

OPTIONS
    --infile,   -i      Input file with Kegg Pathways from list of Kegg API.    REQUIRED
    --outfile,  -o      GMT file (pathways and corresponding KOs)               REQUIRED
    --help,     -h      This help.
    --license   -l      License.

EOF
}

sub license{
    print STDERR <<EOF;
Copyright (C) 2016 Renato Augusto Correa dos Santos
http://bce.bioetanol.org.br
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
