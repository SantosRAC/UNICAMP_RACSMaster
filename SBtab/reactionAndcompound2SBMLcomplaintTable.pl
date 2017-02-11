#!/usr/bin/perl

use warnings;
use strict;

my $reactionTable=$ARGV[0]; # Reaction SBtab table
my $compoundTable=$ARGV[1]; # Compound SBtab table

open(REACTAB,$reactionTable);

while(<REACTAB>){
 chomp;
 if(/^!!SBtab/){
  print "$_\n";
 }
 elsif(/^!ID/){
  print "$_\t!SBML:reaction:id\n";
 } else {
  my ($ID,$ReactionFormula,$IdentifiersKeggReaction,$IdentifiersBiocyc,$GeneSymbol,$IsReversible) = split(/\t/,$_);
  if (($IdentifiersBiocyc) and ($IdentifiersKeggReaction)) {
   print "$ID,$ReactionFormula,$IdentifiersKeggReaction,$IdentifiersBiocyc,$GeneSymbol,$IsReversible\t$ID\n";
   die "Hey, something is wrong. Should not have both Kegg and Biocyc IDs\n";
  } elsif ($IdentifiersKeggReaction) {
   print "$ID\t$ReactionFormula\t$IdentifiersKeggReaction\t\t$GeneSymbol\t$IsReversible\t$ID\n";
  } elsif ($IdentifiersBiocyc) {
   print "$ID\t$ReactionFormula\t\t$IdentifiersBiocyc\t$GeneSymbol\t$IsReversible\t$ID\n";
  } else {
   die "Hey, something is wrong. None of the identifiers (Kegg or Biocyc) exist.\n";
  }
 }
}

close(REACTAB);

## Expected fields in compound SBtab table
## !ID     !Name   !Identifiers:kegg.compound      !Identifiers:bigg.metabolite    !Identifiers:biocyc     !StructureFormula

open(COMPTABLE,$compoundTable);

while(<COMPTABLE>){
 chomp;
 if (/^!!SBtab/){
  print "$_\n";
 } elsif (/^!ID/) {
  print "$_\t!SBML:species:id\n";
 } else {
 my ($ID,$Name,$IdentifiersKeggCompound,$IdentifiersBiggMetabolite,$IdentifiersBiocyc,$StructureFormula) = split(/\t/,$_);
 if (($IdentifiersKeggCompound) and ($IdentifiersBiocyc)) {
   die "Hey, something is wrong. Should not have both Kegg and Biocyc IDs: $IdentifiersKeggCompound and $IdentifiersBiocyc\n";
  } elsif ($IdentifiersKeggCompound) {
   if ($StructureFormula) {
    print "$ID\t$Name\t$IdentifiersKeggCompound\t$IdentifiersBiggMetabolite\t\t$StructureFormula\t$ID\n";
   } else {
    print "$ID\t$Name\t$IdentifiersKeggCompound\t$IdentifiersBiggMetabolite\t\t\t$ID\n";
   }
  } elsif ($IdentifiersBiocyc) {
   if ($StructureFormula) {
    print "$ID\t$Name\t\t$IdentifiersBiggMetabolite\t$IdentifiersBiocyc\t$StructureFormula\t$ID\n";
   } else {
    print "$ID\t$Name\t\t$IdentifiersBiggMetabolite\t$IdentifiersBiocyc\t\t$ID\n";
   }
  } else {
   die "Hey, something is wrong. None of the identifiers (Kegg or Biocyc) exist.\n";
  }
 }
}

close(COMPTABLE);
