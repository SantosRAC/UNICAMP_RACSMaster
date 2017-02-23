#!/usr/bin/perl

use warnings;
use strict;

my $fileCompoundsModel=$ARGV[0]; # File with BiGG identifiers for compounds
my $mappingKeggCompoundNames=$ARGV[1]; # File with KEGG compound ID and corresponding descriptive name

my %keggCompoundID2name;

open(MAPKEGGNAME,$mappingKeggCompoundNames);

while(<MAPKEGGNAME>){
 chomp;
 my ($keggID,$keggName) = split(/\t/,$_);
 if($keggCompoundID2name{$keggID}){
  die "Kegg Identifier for compound already in hash: $keggID\n";
 } else {
  $keggCompoundID2name{$keggID}=$keggName;
 }
}

close(MAPKEGGNAME);

print "!!SBTab\tSBtabVersion=\’1.0\’\tTableName=\’K. brasiliensis SBTab Compound Table\’\tTableType=\’Compound\’
!ID\t!Name\t!Identifiers:kegg.compound\t!Identifiers:bigg.metabolite\t!Identifiers:biocyc\n";

my $compoundSBTabID=1;

open(COMPOUNDSMODEL,$fileCompoundsModel);

while(<COMPOUNDSMODEL>){
 chomp;
 next if (/^IdentifierInModel/);
 my ($modelID,$oldBiGGId,$desc,$newBiGGId,$booleanRandomIDcreation)=split(/\t/,$_);
 if ($modelID =~ /C\d\d\d\d\d/){
  if ($desc eq 'blank_description') {
   if ($keggCompoundID2name{$modelID}) {
    $desc=$keggCompoundID2name{$modelID};
    print "sbtabcomp$compoundSBTabID\t$desc\t$modelID\t$newBiGGId\t\n";
   } else {
    die "It seems that $modelID has no description and no name is available either (mapping file)\n";
   }
  } else {
   print "sbtabcomp$compoundSBTabID\t$desc\t$modelID\t$newBiGGId\t\n";
  }
 } else {
  print "sbtabcomp$compoundSBTabID\t$desc\t\t$newBiGGId\t$modelID\n";
 }
 $compoundSBTabID++;
}

close(COMPOUNDSMODEL);
