#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use XML::Twig;

my @genesExclusiveUmaydis=();
my %geneAssociation;

my $metabolicModel = $ARGV[0]; # U. maydis metabolic model provided by Christian Lieven
my $listOfGenesExclusivelyInUmaydis = $ARGV[1]; # List of genes exclusively found in U. maydis compared to PRIAM
my $geneAssociationFile = $ARGV[2]; # List with gene association between K. brasiliensis and Ustilago maydis, output from the script filterBlast.py (repo)

# Open list of genes exclusively found in U. maydis (compared to PRIAM) and store in a list
open(LISTOFEXCGENES,$listOfGenesExclusivelyInUmaydis);
while(<LISTOFEXCGENES>) {
 chomp;
 my $gene = $_;
 unless($gene ~~ @genesExclusiveUmaydis){
  push(@genesExclusiveUmaydis,$gene);
 }
}
close(LISTOFEXCGENES);

# Read the file with association betwenn proteins in K. brasiliensis and U. maydis and store in a hash
open(GENEASSOCFILE,$geneAssociationFile);
while(<GENEASSOCFILE>){
 chomp;
 my (undef,$maydisGene,$kbrGene)=split(/\t/,$_);
 $maydisGene = lc($maydisGene);
 $geneAssociation{$maydisGene}=$kbrGene;
}
close(GENEASSOCFILE);

print "Reaction_Name\tReaction_Genes\tReaction_Boolean_Model\tReaction_Reversibility\tGenes_Kalmanozyma\tReaction_Reactants\tReaction_Products\tEC_number\tBiocyc\tSubsystem\n";

# Getting models in SBML file
# This script was designed to work with the SBML2 metabolic model of U. maydis provided by MSc Christian Lieven
my $twig=XML::Twig->new();
$twig->parsefile($metabolicModel) or die "cannot parse $metabolicModel";
my $root = $twig->root;
my @models = $root->children('model');

my %speciesInfo;
my %speciesKeggID;

foreach my $mod (@models){
 # Getting information about each species in model
 my @listOfSpecies = $mod->children('listOfSpecies');
 foreach my $allSpecies (@listOfSpecies) {
  my @speciesThemselves = $allSpecies->children('species');
  foreach my $species (@speciesThemselves) {
   my $speciesID = $species->att('id');
   my $speciesName = $species->att('name');
   if((not $speciesID) or (not $speciesName)) { die "Species name or identifier don't exist. Check model.\n" }
   $speciesInfo{$speciesID}=$speciesName;
   my @notes = $species->children('notes');
   foreach my $note (@notes) {
    my @bodies = $note->children('body');
    foreach my $body (@bodies) {
     my @paragraphs = $body->children('p');
     foreach my $p (@paragraphs) {
      my $textParagraph = $p->text();
      if ($textParagraph =~ /KEGG:/) {
       my $speKeggID = $textParagraph;
       $speKeggID =~ s/KEGG: //g;
       #$speciesID =~ s/__45__/-/g;
       #$speciesID =~ s/__43__/+/g;
       #$speciesID =~ s/_c$//g;
       #$speciesID =~ s/_e$//g;
       #$speciesID =~ s/_p$//g;
       #$speciesID =~ s/_m$//g;
       $speciesKeggID{$speciesID}=$speKeggID;
      }
     }
    }
   }
  }
 }
 # Getting information about each reaction in U. maydis metabolic model
 my @listOfReactions = $mod->children('listOfReactions');
 foreach my $allReactions (@listOfReactions) {
  my @reactionsThemselves = $allReactions->children('reaction');
  foreach my $reaction (@reactionsThemselves) {
   my $reactionName = $reaction->att('id');
   my $reactionReversibility = $reaction->att('reversible');
   my @reacGenes=();
   my $reacBoolean='';
   my @reacReactants=();
   my @reacProducts=();
   my $reacECnumber='';
   my $reactBiocyc='';
   my $reacSubsystem='';
   my @genesKbr=();
   my @notes = $reaction->children('notes');
   foreach my $note (@notes) {
    my @bodies = $note->children('body');
    foreach my $body (@bodies) {
     my @paragraphs = $body->children('p');
     foreach my $p (@paragraphs) {
      my $textParagraph = $p->text();
      if ($textParagraph =~ /BIOCYC:/) {
       $reactBiocyc = $textParagraph;
       $reactBiocyc =~ s/BIOCYC: //g;
      } elsif ($textParagraph =~ /SUBSYSTEM:/) {
       $reacSubsystem = $textParagraph;
       $reacSubsystem =~ s/SUBSYSTEM: //g;
      } elsif ($textParagraph =~ /GENE_ASSOCIATION:/) {
       @reacGenes = $textParagraph =~ /\((um\d\d\d\d\d)\)/g;
       if ($textParagraph =~ /and/) {
        $reacBoolean='and';
       } elsif ($textParagraph =~ /or/) {
        $reacBoolean='or';
       } else {
        $reacBoolean='no';
       }
      } elsif ($textParagraph =~ /EC Number:/) {
       $reacECnumber = $textParagraph;
       $reacECnumber =~ s/EC Number: //g;
      } elsif ($textParagraph =~ /Confidence level:/) {
       #TODO
      } elsif ($textParagraph =~ /Note:/) {
       #TODO
      } elsif ($textParagraph =~ /Notes:/) {
       #TODO
      } else {
       #TODO
      }
     }
    }
   }
   my @reactants = $reaction->children('listOfReactants');
   foreach my $react (@reactants) {
    my @species = $react->children('speciesReference');
    foreach my $spe (@species) {
     my $speID=$spe->att('species');
     my $speStoic=$spe->att('stoichiometry');
     if($speciesKeggID{$speID}) {
      my $finalSpecies = "$speciesKeggID{$speID}".":$speStoic";
      push(@reacReactants,$finalSpecies);
     } else {
      my $finalSpecies = "$speID".":$speStoic";
      push(@reacReactants,$finalSpecies);
     }
    }
   }
   my @products = $reaction->children('listOfProducts');
   foreach my $prod (@products) {
    my @species = $prod->children('speciesReference');
    foreach my $spe (@species) {
     my $speID=$spe->att('species');
     my $speStoic=$spe->att('stoichiometry');
     if($speciesKeggID{$speID}) {
      my $finalSpecies = "$speciesKeggID{$speID}".":$speStoic";
      push(@reacProducts,$finalSpecies);
     } else {
      my $finalSpecies = "$speID".":$speStoic";
      push(@reacProducts,$finalSpecies);
     }
    }
   }
   my $booleanFoundUmaydisExclusive='';
   foreach my $geneThisReaction (@reacGenes){
    my $kbrGene = $geneAssociation{$geneThisReaction};
    if ($kbrGene ~~ @genesExclusiveUmaydis) {
     $booleanFoundUmaydisExclusive=1;
     push(@genesKbr,$kbrGene);
    }
   }
   if($booleanFoundUmaydisExclusive) {
    print "$reactionName\t".join(",",@reacGenes)."\t$reacBoolean\t$reactionReversibility\t".join(",",@genesKbr)."\t".join(",",@reacReactants)."\t".join(",",@reacProducts)."\t$reacECnumber\t$reactBiocyc\t$reacSubsystem\n";
   }
  }
 }
}
