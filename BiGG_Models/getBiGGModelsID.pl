#!/usr/bin/perl

use warnings;
use strict;

#TODO Command-line line arguments

my $metabolitesBiGGModels = $ARGV[0];
my $listKeggIdentifiersFile = $ARGV[1];

# Variables
my %biggID2KeggID;

# Read BiGG Models metabolite file
open(BIGGMODELSMET,$metabolitesBiGGModels);

while(<BIGGMODELSMET>){
 chomp;
 next if (/^bigg_id/);
 my ($biggId,$universalBiggId,$name,$modelList,$databaseLinks)=split(/\t/,$_);
 if($databaseLinks =~ /\{\"link\": \"http:\/\/identifiers\.org\/kegg\.compound\/C\d\d\d\d\d", \"id\": \"(C\d\d\d\d\d)\"\}/) {
  my @keggMatches = $databaseLinks =~ /\{\"link\": \"http:\/\/identifiers\.org\/kegg\.compound\/C\d\d\d\d\d", \"id\": \"(C\d\d\d\d\d)\"\}/g;
  foreach my $keggID (@keggMatches) {
   if(($biggID2KeggID{$keggID})){
    if($universalBiggId ~~ @{$biggID2KeggID{$keggID}}) {
     next;
    } else {
     push(@{$biggID2KeggID{$keggID}},$universalBiggId);
    }
   } else {
    @{$biggID2KeggID{$keggID}}=($universalBiggId);
   }
  }
 }
}

close(BIGGMODELSMET);

open(LISTKEGGIDS,$listKeggIdentifiersFile);

while(<LISTKEGGIDS>){
 chomp;
 my $keggID='';
 if(/^C\d\d\d\d\d$/){
  $keggID=$_;
 } else {
  die "Observed identifier is not presented as expected\n";
 }
 if ($biggID2KeggID{$keggID}) {
  print "$keggID\t".join(",",@{$biggID2KeggID{$keggID}})."\n";
 } else {
  print "$keggID\tNO\n";
 }
}

close(LISTKEGGIDS);
