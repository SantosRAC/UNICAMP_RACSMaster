#!usr/bin/perl

use warnings;
use strict;
use diagnostics;
use XML::Twig;

my $infile=$ARGV[0]; # SBML3

open(INFILE,$infile);

my $twig=XML::Twig->new();
$twig->parsefile($infile) or die "cannot parse $infile";
my $root = $twig->root;
my @models = $root->children('model');

print "!!SBtab\tSBtabVersion=\’1.0\’\tTableType=\’Reaction\’\tTableName=\’K. brasiliensis SBTab Reaction Table\’
!ID\t!ReactionFormula\t!Identifiers:kegg.reaction\t!Identifiers:biocyc\t!Gene:Symbol\t!IsReversible\n";

my $reactionSBTabID=1;
my $geneNumberTest=1; #TODO: Remove this "geneNumberTest" and include real genes (more complex situations, in which it is necessary to create Enzyme/ Gene tables)

# Parsing reaction list
foreach my $mod (@models) {
 my @reactionLists = $mod->children('listOfReactions');
 foreach my $real (@reactionLists) {
  my @allReactions = $real->children('reaction');
  foreach my $rea (@allReactions) {
   my @reactantsReaction=();
   my @productsReaction=();
   my $reactID=$rea->att('id');
   my $reactMetaID=$rea->att('metaid');
   my $reversibilityReaction=$rea->att('reversible');
   $reversibilityReaction = ucfirst($reversibilityReaction);
   # Parsing reactant list in reaction 
   my @reactantList = $rea->children('listOfReactants');
   foreach my $reactl (@reactantList) {
    my @allReactants = $reactl->children('speciesReference');
    foreach my $reaspe (@allReactants) {
     my $spID = $reaspe->att('species');
     unless ($spID ~~ @reactantsReaction) {
      push(@reactantsReaction,$spID)
     }
    }
   }
   # Parsing product list in reaction 
   my @productList = $rea->children('listOfProducts');
   foreach my $productl (@productList) {
    my @allProducts = $productl->children('speciesReference');
    foreach my $prospe (@allProducts) {
     my $spID = $prospe->att('species');
     unless ($spID ~~ @productsReaction) {
      push(@productsReaction,$spID);
     }
    }
   }
   my $sbtabReactID='sbtabreact'.$reactionSBTabID;
   my $geneTestID='testGeneSymbol'.$geneNumberTest; #TODO: Remove this line in the future (Gene and Enzyme tables)
   if($reactID =~ /R\d\d\d\d\d/){
    print "$sbtabReactID\t".join(" + ",@reactantsReaction)." <=> ".join(" + ",@productsReaction)."\t$reactID\t\t$geneTestID\t$reversibilityReaction\n";
   } else {
    print "$sbtabReactID\t".join(" + ",@reactantsReaction)." <=> ".join(" + ",@productsReaction)."\t\t$reactID\t$geneTestID\t$reversibilityReaction\n";
   }
   $reactionSBTabID++;
   $geneNumberTest++; # TODO: Remove this line in the future (Gene and Enzyme tables)
  }
 }
}

close(INFILE);
