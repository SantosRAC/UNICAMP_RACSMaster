#!/usr/bin/perl

use warnings;
use strict;
use LWP::Simple;

my $listOfKeggCompounds = $ARGV[0];
my %keggCompoundFormulae;
my $compoundCount=1;

open(LISTOFCOMPOUNDS,$listOfKeggCompounds);

while(<LISTOFCOMPOUNDS>) {
 chomp;
 my $keggID='';
 my $keggIDformulae='';
 $keggID=$_;
 unless($keggID =~ /^C\d\d\d\d\d$/){ die "Compound identifier does not appear as expeced in line $_\n"; }
 if ($keggCompoundFormulae{$keggID}) {
  die "Something wrong with $keggID (ID appears to be repeated in input file!)\n";
 } else {
  my $kegg_url = "http://rest.kegg.jp/get/$keggID";
  my $content = get($kegg_url);
  if($content =~ /ENTRY\s+(C\d\d\d\d\d)\s+Compound/) {
   my @contentLines = split(/\n/,$content);
   foreach my $line (@contentLines){
    if ($line =~ /^FORMULA\s+(\S+)$/) {
     $keggIDformulae=$1;
     if ($keggCompoundFormulae{$keggID}) {
      die "It seems that you already included a formula for this compound ($keggID)\nPlease, check you input!\n";
     } else {
      $keggCompoundFormulae{$keggID}=$keggIDformulae;
      print "$keggID\t$keggIDformulae\n";
     }
    }
   }
  } else {
   die "It couldn't be possible to find this kegg identifier ($keggID) on REST page\n$content\n";
  }
 }
 sleep 2;
}

close(LISTOFCOMPOUNDS);
