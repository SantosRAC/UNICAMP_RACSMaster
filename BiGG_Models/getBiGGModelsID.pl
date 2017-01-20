#!/usr/bin/perl

use warnings;
use strict;

#TODO Command-line line arguments

my $metabolitesBiGGModels = $ARGV[0];
my $listModelIdentifiersFile = $ARGV[1];

# Variables
my %biggID2KeggID;
my %biggID2cycID;

# Read BiGG Models metabolite file
open(BIGGMODELSMET,$metabolitesBiGGModels);

while(<BIGGMODELSMET>){
 chomp;
 next if (/^bigg_id/);
 my ($biggId,$universalBiggId,$name,$modelList,$databaseLinks)=split(/\t/,$_);
 if($databaseLinks =~ /\{\"link\": \"http:\/\/identifiers\.org\/kegg\.compound\/C\d\d\d\d\d\", \"id\": \"(C\d\d\d\d\d)\"\}/) {
  my @keggMatches = $databaseLinks =~ /\{\"link\": \"http:\/\/identifiers\.org\/kegg\.compound\/C\d\d\d\d\d\", \"id\": \"(C\d\d\d\d\d)\"\}/g;
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
 if(/\{\"link\": \"http:\/\/identifiers\.org\/biocyc\/\S+\", \"id\": \"(\S+)\"\}/){
  my @cycMatches = $databaseLinks =~ /\{\"link\": \"http:\/\/identifiers\.org\/biocyc\/\S+\", \"id\": \"(\S+)\"\}/g;
  foreach my $cycID (@cycMatches) {
   if(($biggID2cycID{$cycID})){
    if($universalBiggId ~~ @{$biggID2cycID{$cycID}}) {
     next;
    } else {
     push(@{$biggID2cycID{$cycID}},$universalBiggId);
    }
   } else {
    @{$biggID2cycID{$cycID}}=($universalBiggId);
   }
  }
 }
}

close(BIGGMODELSMET);

open(LISTMODELIDS,$listModelIdentifiersFile);

while(<LISTMODELIDS>){
 chomp;
 my $keggID='';
 my $cycID='';
 if(($_ =~ /^C\d\d\d\d\d$/) or ($_ =~ /^G\d\d\d\d\d$/)){
  $keggID=$_;
  if ($biggID2KeggID{$keggID}) {
   print "$keggID\t".join(",",@{$biggID2KeggID{$keggID}})."\n";
  } else {
   print "$keggID\tNO_BiGG\n";
  }
 } else {
  $cycID=$_;
  if ($biggID2cycID{$cycID}) {
   print "$cycID\t".join(",",@{$biggID2KeggID{$keggID}})."\n";
  } else {
   print "$cycID\tNO_BiGG_CHECK_METACYCID\n";
  }
 }
}

close(LISTMODELIDS);
