#!/usr/bin/perl
#TODO: Add validation of SBML model using Online validator API: http://sbml.org/Facilities/Validator/ and LWP (Perl)
require v5.10.1;
use strict;
use warnings;
use diagnostics;
#use LWP::UserAgent; # This will be used when the Perl code to connect to SBML validator is available
use Getopt::Long;
use XML::Twig;
use XML::Writer; # helper module for Perl programs that write an XML document.
use HTML::Entities; # Encode or decode strings with HTML entities

my $version='0.1';
my $sbmlLevel='';
my $sbmlVersion='';
my $infile='';
my $outfile='';
my $license='';
my $help='';
my $debug=0;
my $reactionTable='';
#my $ua = LWP::UserAgent->new; # This will be used only when the online validator is available

##############################################
#Variables related to info in existing model
##############################################
my @compartmentsInModel=();
my @reactionsInModel=();
my @metabolitesInModel=();
my @genesInModel=();
# Units in model
my %unitsModelInfo;
my %unitDefinitionModelInfo;
# Objectives in model
my %objectivesInModelInfo;
my %fluxObjectivesInModelInfo;
# Species, reactions
my %speciesModelInfo;
my %reactionsModelInfo;
my @allowedCompartments=('UNK_COMP','cytoplasm','Cytosol','mitochondrion','Mitochondria','peroxisome','Golgi','Golgi_Apparatus','vacuole','ER','Endoplasmic_Reticulum','plasma_membrane','nucleus','Periplasm','extracellular','Vacuole','Nucleus','Peroxisome','Extra_organism','Extraorganism');
my %infoSBMLLine;
my %infoModelLine;

##############################################
#Variables for reactions in input table
##############################################
my %inputTableReactionProductsStoic;
my %inputTableReactionReactantsStoic;
my %inputTableReactionProductsComp;
my %inputTableReactionReactantsComp;
my %inputTableReactionProductsUpperBoundFlux;
my %inputTableReactionReactantsUpperBoundFlux;
my %inputTableReactionProductsLowerBoundFlux;
my %inputTableReactionReactantsLowerBoundFlux;
my %inputTableReactionReactionType;
my @reactionInTable=();
my %reaction2productsInTable;
my %reaction2reactantsInTable;

GetOptions(
    'license|l'       => \$license,
    'help|h|?'        => \$help,
    'infile|i=s'      => \$infile,
    'input_reaction_table|t=s' => \$reactionTable,
);

if($help){
        &usage();
    exit 1;
}
if($license){
        &license();
    exit 1;
}
if(!-s $infile){
    print STDERR "FATAL: you must provide a SBML file (level 3 version 1) as input\n";
    &usage();
    exit 0;     
}
if(!-s $reactionTable){
    print STDERR "FATAL: you must provide a table with information about the reaction(s) and metabolite(s)\n";
    &usage();
    exit 0;
}

parseXML();
getInputTable();
checkInputTable();

##############################################
#Check input table with information about
#reactions and metabolites
##############################################
sub getInputTable {
 open(INPUTTABLE,$reactionTable);
 while(<INPUTTABLE>){
  chomp;
  next if (/^#/);
  my ($reactionBiGGID,$reactionBiGGName,$reactionType,$reactionECorTCDB,$transportBoolean,$reversible,$metaboliteBiGGID,$metaboliteType,$metStoichiometry,$metCompartment,$geneAssociation,$upperBoundFlux,$lowerBoundFlux) = split(/\t/,$_);

  #Checks BiGG standards for reaction and pseudo-reaction in input table
  $reactionType=lc($reactionType);
  if($reactionType eq 'reaction'){
   $inputTableReactionReactionType{$reactionBiGGID}=$reactionType;
   if($reactionBiGGID =~ /R_(\S+)?/){
    print "Reaction \"$reactionBiGGID\" is OK!\n";
    unless ($reactionBiGGID ~~ @reactionInTable) {push(@reactionInTable,$reactionBiGGID);}
   } else {
    print "Reaction \"$reactionBiGGID\" does not follow BiGG standard\n";
   }
  }
  elsif(($reactionType eq 'pseudoreaction') or ($reactionType eq 'pseudo-reaction')){
   $inputTableReactionReactionType{$reactionBiGGID}=$reactionType;
   if(($reactionBiGGID =~ /EX_(\S+)?/) or ($reactionBiGGID =~ /DM_(\S+)?/) or ($reactionBiGGID =~ /SK_(\S+)?/) or ($reactionBiGGID =~ /(R_)?ATPM/)){
    print "Pseudo-reaction \"$reactionBiGGID\" is OK!\n";
    unless ($reactionBiGGID ~~ @reactionInTable) {push(@reactionInTable,$reactionBiGGID);}
   } else {
    print "Pseudo-reaction \"$reactionBiGGID\" does not follow BiGG standard\n";
   }
  } else {
   die "Reaction \"$reactionBiGGID\" is neither a reaction nor a pseudo-reaction (ATPM, biomass, exchange, sink, or demand reactions)\n";
  }

  #Checks BiGG standards for metabolite in input table
  if($metaboliteBiGGID =~ /M_(\S+)(_\S+)?/){
   if($metaboliteBiGGID ~~ @metabolitesInModel){print "Metabolite ($metaboliteBiGGID) already exists in model\nInformation for new reaction ($reactionBiGGID) will be retrieved from existing metabolite.\n";} # only information specific for metabolite, not those specific for reactions.
   $metaboliteType=lc($metaboliteType);
   if($metaboliteType eq 'product'){
    unless($metaboliteBiGGID ~~ @{$reaction2productsInTable{$reactionBiGGID}}){
     push(@{$reaction2productsInTable{$reactionBiGGID}},$metaboliteBiGGID);
    }
    $inputTableReactionProductsStoic{$reactionBiGGID}{$metaboliteBiGGID}=$metStoichiometry;
    $inputTableReactionProductsComp{$reactionBiGGID}{$metaboliteBiGGID}=$metCompartment;
    $inputTableReactionProductsUpperBoundFlux{$reactionBiGGID}{$metaboliteBiGGID}=$upperBoundFlux;
    $inputTableReactionProductsLowerBoundFlux{$reactionBiGGID}{$metaboliteBiGGID}=$lowerBoundFlux;
   }
   elsif($metaboliteType eq 'reactant'){
    unless($metaboliteBiGGID ~~ @{$reaction2reactantsInTable{$reactionBiGGID}}){
     push(@{$reaction2reactantsInTable{$reactionBiGGID}},$metaboliteBiGGID);
    }
    $inputTableReactionReactantsStoic{$reactionBiGGID}{$metaboliteBiGGID}=$metStoichiometry;
    $inputTableReactionReactantsComp{$reactionBiGGID}{$metaboliteBiGGID}=$metCompartment;
    $inputTableReactionReactantsUpperBoundFlux{$reactionBiGGID}{$metaboliteBiGGID}=$upperBoundFlux;
    $inputTableReactionReactantsLowerBoundFlux{$reactionBiGGID}{$metaboliteBiGGID}=$lowerBoundFlux;
   } else {
    die "Metabolite BiGG ID \"$metaboliteBiGGID\" is neither a \'reactant\' nor a \'product\'\n";
   }
  } else {
   die "Metabolite in input table \"$metaboliteBiGGID\" does not follow the BiGG standards\n";
  }

 }
 close(INPUTTABLE);
}

##############################################
#Main function that parses the XML file
##############################################

sub parseXML {

 # Create a Twig object to represent the XML
 my $twig=XML::Twig->new();
 $twig->parsefile($infile) or die "cannot parse $infile";
 my $root = $twig->root;

 # Lines below get information from the line with 'sbml' tag
 my $sbmlLevel='';
 if($root->att('level')){
  $sbmlLevel=$root->att('level');
  $infoSBMLLine{'level'}=$sbmlLevel;
 }

 my $sbmlVersion='';
 if($root->att('version')){
  $sbmlVersion=$root->att('version');
  $infoSBMLLine{'version'}=$sbmlVersion;
 }

 my $xmlnsFBC='';
 if($root->att('xmlns:fbc')){
  $xmlnsFBC = $root->att('xmlns:fbc');
  $infoSBMLLine{'xmlns:fbc'}=$xmlnsFBC;
 }

 my $sboTerm='';
 if($root->att('sboTerm')){
  $sboTerm = $root->att('sboTerm');
  $infoSBMLLine{'sboTerm'}=$sboTerm;
 }

 my $xmlns='';
 if($root->att('xmlns')){
  $xmlns = $root->att('xmlns');
  $infoSBMLLine{'xmlns'}=$xmlns;
 }

 my $fbcRequired='';
 if($root->att('fbc:required')){
  $fbcRequired = $root->att('fbc:required');
  $infoSBMLLine{'fbc:required'}=$fbcRequired;
 }

 # Check if user input is a SBML Level 3, Version 1
 if (($sbmlLevel == 3) and ($sbmlVersion == 1)) {
  if ($debug > 0) { print "SBML SPEC: Level $sbmlLevel, version $sbmlVersion\n";}
 } else {
  die "You must use SBML Level 3, version 1\nCheck specitifications at www.sbml.org\n";
 }

 # Check if there is a model (required inside sbml)
 if (!$root->children('model')){ die "SBML Level 3 Version 1 must have a model object inside <sbml> ... </sbml>!\n"; }
 
 my @models = $root->children('model');

 foreach my $mod (@models){

  my $ModelID='';
  if($mod->att('id')){
   $ModelID = $mod->att('id');
   $infoModelLine{'id'}=$ModelID;
  }

  my $ModelFBCStrict='';
  if($mod->att('fbc:strict')){
   $ModelFBCStrict = $mod->att('fbc:strict');
   $infoModelLine{'fbc:strict'}=$ModelID;
  }
  
  #Checking list of unit definitions
  my @listOfUnitDefinitions = $mod->children('listOfUnitDefinitions');
  foreach my $unitDefinition (@listOfUnitDefinitions){
   my @listOfUnits = $unitDefinition->children('unitDefinition');
   my $countUnitsDef=0;
   foreach my $unit (@listOfUnits){
    $unitDefinitionModelInfo{$countUnitsDef}=$unit->att('id');
    $countUnitsDef++;
    my @listOfUnits2 = $unit->children('listOfUnits');
    foreach my $unit2 (@listOfUnits2){
     my @listOfUnits3 = $unit2->children('unit');
     my $countUnits=0;
     foreach my $unit3 (@listOfUnits3){
      $unitsModelInfo{$countUnits}{'exponent'}=$unit3->att('exponent');
      $unitsModelInfo{$countUnits}{'kind'}=$unit3->att('kind');
      $unitsModelInfo{$countUnits}{'multiplier'}=$unit3->att('multiplier');
      $unitsModelInfo{$countUnits}{'scale'}=$unit3->att('scale');
      $countUnits++;
     }
    }
   }
  }

  #Checking list of objectives in model
  my @listOfObjectives = $mod->children('fbc:listOfObjectives');
  foreach my $objectiveListUnit (@listOfObjectives){
   if($objectiveListUnit->att('fbc:activeObjective') eq 'obj'){
    print "Objective fbc:activeObjective: \'".$objectiveListUnit->att('fbc:activeObjective')."\'\n";
   } else {
    die "Objective fbc:activeObjective attribute is not as expected.\n";
   }
   my @objectives = $objectiveListUnit->children('fbc:objective');
   my $objectiveCount=1;
   foreach my $objective (@objectives){
    $objectivesInModelInfo{$objectiveCount}{'fbc:id'}=$objective->att('fbc:id');
    $objectivesInModelInfo{$objectiveCount}{'fbc:type'}=$objective->att('fbc:type');
    my @listOfFluxObjectives = $objective->children('fbc:listOfFluxObjectives');
    foreach my $fluxObjective (@listOfFluxObjectives){
     my @fluxObjectives = $fluxObjective->children('fbc:fluxObjective');
     my $foCount=1;
     foreach my $fo (@fluxObjectives){
      $fluxObjectivesInModelInfo{$objectiveCount}{$foCount}{'fbc:coefficient'}=$fo->att('fbc:coefficient');
      $fluxObjectivesInModelInfo{$objectiveCount}{$foCount}{'fbc:reaction'}=$fo->att('fbc:reaction');
      #TODO: Check if the reaction is in model. If it is not, complain!!
      print "$objectivesInModelInfo{$objectiveCount}{'fbc:id'}\t$objectivesInModelInfo{$objectiveCount}{'fbc:type'}\t$fluxObjectivesInModelInfo{$objectiveCount}{$foCount}{'fbc:coefficient'}\t$fluxObjectivesInModelInfo{$objectiveCount}{$foCount}{'fbc:reaction'}\n";
     }
    }
    $objectiveCount++;
   }
  }

  #Checking FBC list of gene products
  #The model must have at least one, according to the fbc package specification
  my @GeneProductLists = $mod->children('fbc:listOfGeneProducts');
  foreach my $geneProduct (@GeneProductLists){
   my @geneProducts=$geneProduct->children('fbc:geneProduct');
   foreach my $gp (@geneProducts){
    my $gpID='';
    my $gpLabel='';
    my $gpName='';
    $gpID = $gp->att('fbc:id');
    $gpLabel = $gp->att('fbc:label');
    $gpName = $gp->att('fbc:name');
    # check if gene identifier starts with 'G_' (required by BiGG, paper of 2015)
    if($gpID =~ /^G_(\S+)$/){
     if($gpID ~~ @genesInModel){
      #TODO
     } else {
      push(@genesInModel,$gpID);
     }
    } else {
     die "Gene $gpName does not follow required regular expression in BiGG paper (check King et al, 2015).\n";
    }
   }
  }

  # Parsing compartments
  my @compartmentLists = $mod->children('listOfCompartments');
  foreach my $compl (@compartmentLists) {
   my @AllCompartments = $compl->children('compartment');
   foreach my $comp (@AllCompartments) {
    my $compId='';
    my $compName='';
    $compId = $comp->att('id');
    $compName = $comp->att('name');
    if (($compId ~~ @allowedCompartments) or ($compName ~~ @allowedCompartments)){
     push(@compartmentsInModel,$compId);
    } else {
     if($compId and $compName){
      die "Compartment (ID: $compId, name: $compName) present in model is not among allowed compartments!\n";
     } else {
      die "Compartment (ID: $compId) present in model is not among allowed compartments!\n";
     }
    }
   }
  }
 
  # Parsing species in the model
  my @speciesList = $mod->children('listOfSpecies');
  foreach my $spemodell (@speciesList) {
   my @AllSpeciesModel = $spemodell->children('species');
   foreach my $spemodel (@AllSpeciesModel) {

    #######################################
    # Variables for each species
    #######################################
    my $modSpeciesID = $spemodel->att('id');
    if ($modSpeciesID =~ /M_(\S+)(_\S+)?/){
     if($modSpeciesID ~~ @metabolitesInModel){
      die "There is a duplicated metabolite in existing model: $modSpeciesID. Please check this!\n";
     } else {
      push(@metabolitesInModel,$modSpeciesID);
     }
    } else {
     die "Identifier for metabolite in model is not regular BiGG ID: $modSpeciesID\n";
    }
    my $modSpeciescharge = $spemodel->att('fbc:charge');
    my $modSpecieshasOnlySubstanceUnits = $spemodel->att('hasOnlySubstanceUnits');
    my $modSpeciesConstant = $spemodel->att('constant');
    my $modSpeciesboundaryCondition = $spemodel->att('boundaryCondition');
    my $modSpeciesName = $spemodel->att('name') ? $spemodel->att('name') : $spemodel->att('id');
    my $modSpeciesCompartment = $spemodel->att('compartment');
    my $modeSpeciesMetaID = $spemodel->att('metaid');
    my $annotResource='';

    # Parsing annotations for a given species in the model
    my @spInModelAnnotations = $spemodel->children('annotation');
    foreach my $annot (@spInModelAnnotations) {
     my $firstAnnotChild = $annot->first_child();
     my $secondAnnotChild = $firstAnnotChild->first_child();
     my $thirdAnnotChild = $secondAnnotChild->first_child();
     my $fourthAnnotChild = $thirdAnnotChild->first_child();
     my $fifthAnnotChild = $fourthAnnotChild->first_child();
     if($fifthAnnotChild->name =~ /rdf:li/) {
      $annotResource=$fifthAnnotChild->att('rdf:resource');
     }
    }

    # Adding information about model species in hash
    $speciesModelInfo{$modSpeciesID}{'overall'}{'id'}=$modSpeciesID;
    $speciesModelInfo{$modSpeciesID}{'overall'}{'metaid'}=$modeSpeciesMetaID;
    $speciesModelInfo{$modSpeciesID}{'overall'}{'charge'}=$modSpeciescharge;
    $speciesModelInfo{$modSpeciesID}{'overall'}{'hasOnlySubstanceUnits'}=$modSpecieshasOnlySubstanceUnits;
    $speciesModelInfo{$modSpeciesID}{'overall'}{'boundaryCondition'}=$modSpeciesboundaryCondition;
    $speciesModelInfo{$modSpeciesID}{'overall'}{'constant'}=$modSpeciesConstant;
    $speciesModelInfo{$modSpeciesID}{'overall'}{'name'}=$modSpeciesName;
    $speciesModelInfo{$modSpeciesID}{'overall'}{'compartment'}=$modSpeciesCompartment;
    $speciesModelInfo{$modSpeciesID}{'annotation'}{'resource'}=$annotResource;
   }
  }

  # Parsing reactions
  my @reactionLists = $mod->children('listOfReactions');
  foreach my $real (@reactionLists) {

   my @AllReactions = $real->children('reaction');
   foreach my $rea (@AllReactions) {
    my $modReactID='';
    my $modReactMetaID='';
    my $modReactName='';
    my $modReactFast='';
    my $modReactReversibility='';
    my $modReactLowerBoundary='';
    my $modReactUpperBoundary='';
    my $reactAnnotResource='';
    my @reactantsReaction=();
    my @productsReaction=();

    $modReactID=$rea->att('id');
    if ($modReactID =~ /R_(\S+)?/){
     if($modReactID ~~ @reactionsInModel){
      die "There is a duplicated reaction in existing model: $modReactID. Please check this!\n";
     } else {
      push(@reactionsInModel,$modReactID); 
     }
    } else {
     die "Identifier for reaction in model is not regular BiGG ID: \n";
    }
    $modReactName=$rea->att('name');
    $modReactFast=$rea->att('fast');
    $modReactMetaID=$rea->att('metaid');
    $modReactReversibility=$rea->att('reversible');
    $modReactLowerBoundary=$rea->att('fbc:lowerFluxBound');
    $modReactUpperBoundary=$rea->att('fbc:upperFluxBound');

    # Parsing list of reactants of a given reaction
    my @reactantList = $rea->children('listOfReactants');
    foreach my $reactl (@reactantList) {
     my $spID='';
     my @AllReactants = $reactl->children('speciesReference');
     foreach my $reaspe (@AllReactants) {
      $spID = $reaspe->att('species');
     }
    }

    # Parsing list of products of a given reaction
    my @productList = $rea->children('listOfProducts');
    foreach my $productl (@productList) {
     my $spID='';
     my @AllProducts = $productl->children('speciesReference');
     foreach my $productspe (@AllProducts) {
      $spID = $productspe->att('species');
     }
    }

    # Parsing annotations in reactions
    my @reacInModelAnnotations = $rea->children('annotation');
    foreach my $annot (@reacInModelAnnotations) {
     my $firstAnnotChild = $annot->first_child();
     my $secondAnnotChild = $firstAnnotChild->first_child();
     my $thirdAnnotChild = $secondAnnotChild->first_child();
     my $fourthAnnotChild = $thirdAnnotChild->first_child();
     my $fifthAnnotChild = $fourthAnnotChild->first_child();
     if($fifthAnnotChild->name =~ /rdf:li/) {
      $reactAnnotResource=$fifthAnnotChild->att('rdf:resource');
     }
    }

    # Adding information about model reaction in hash
    $reactionsModelInfo{$modReactID}{'overall'}{'id'}=$modReactID;
    $reactionsModelInfo{$modReactID}{'overall'}{'name'}=$modReactName;
    $reactionsModelInfo{$modReactID}{'overall'}{'fast'}=$modReactFast;
    $reactionsModelInfo{$modReactID}{'overall'}{'metaid'}=$modReactMetaID;
    $reactionsModelInfo{$modReactID}{'overall'}{'resource'}=$reactAnnotResource;
    $reactionsModelInfo{$modReactID}{'overall'}{'reversible'}=$modReactReversibility;
    $reactionsModelInfo{$modReactID}{'overall'}{'lowerboundflux'}=$modReactLowerBoundary;
    $reactionsModelInfo{$modReactID}{'overall'}{'upperboundflux'}=$modReactUpperBoundary;
    @{ $reactionsModelInfo{$modReactID}{'reaction'}{'listOfReactants'} }=@reactantsReaction;
    @{ $reactionsModelInfo{$modReactID}{'reaction'}{'listOfProducts'} }=@productsReaction;

   }
  }
 }

}

#############################################
#Function that checks some things from the input reactions (reaction table)
##############################################

sub checkInputTable {
 foreach my $rct (@reactionInTable){
  print "Reaction: $rct\tReaction type: $inputTableReactionReactionType{$rct}\n";
  my $lbFlux='';
  my $ubFlux='';
 
  #TODO 
  # Control of pseudo-reactions. Defined variables related to pseudo-reactions
  # Exchange reactions have just one metabolite
  # Demand and sink reaction have just one metabolite
  # Biomass reactions have many metaobolite
  # ATP maintenance reaction is a special case

  foreach my $met (@{$reaction2reactantsInTable{$rct}}){
   #Check the lower bound flux
   if($ubFlux){
    if($ubFlux!=$inputTableReactionReactantsUpperBoundFlux{$rct}{$met}){
     die "Something wrong. Different lines for the same reaction with different upper bound flux set\n";
    }
   } else {
    $ubFlux=$inputTableReactionReactantsUpperBoundFlux{$rct}{$met};
   }
   if($lbFlux){
    if($lbFlux!=$inputTableReactionReactantsLowerBoundFlux{$rct}{$met}){
     die "Something wrong. Different lines for the same reaction with different lower bound flux set\n";
    }
   } else {
    $lbFlux=$inputTableReactionReactantsLowerBoundFlux{$rct}{$met};
   }
   print "----|Reactant| Metabolite: $met\tstoichiometry: $inputTableReactionReactantsStoic{$rct}{$met}\tCompartment: $inputTableReactionReactantsComp{$rct}{$met}\tReaction upper bound flux: $inputTableReactionReactantsUpperBoundFlux{$rct}{$met}\tReaction lower bound flux: $inputTableReactionProductsLowerBoundFlux{$rct}{$met}\n";
  }
  foreach my $met (@{$reaction2productsInTable{$rct}}){
   print "----|Product| Metabolite: $met\tstoichiometry: $inputTableReactionProductsStoic{$rct}{$met}\tCompartment: $inputTableReactionProductsComp{$rct}{$met}\tReaction upper bound flux: $inputTableReactionProductsUpperBoundFlux{$rct}{$met}\tReaction lower bound flux: $inputTableReactionProductsLowerBoundFlux{$rct}{$met}\n";
  }
 }
}

#############################################
# Main function that includes reactions in model
# Should use module XML::Writer
# Should the module IO::File to output the resulting model
##############################################

#
my $LISTOFUNITSDEFwriter = XML::Writer->new(OUTPUT => 'self', NEWLINES => 0);
$LISTOFUNITSDEFwriter->startTag('listOfUnitDefinitions');
$LISTOFUNITSDEFwriter->endTag('listOfUnitDefinitions');

my $UNITDEFwriter = XML::Writer->new(OUTPUT => 'self', NEWLINES => 0);
$UNITDEFwriter->startTag('unitDefinition',%unitDefinitionModelInfo);
$UNITDEFwriter->endTag('unitDefinition');

my $LISTOFUNITSwriter = XML::Writer->new(OUTPUT => 'self', NEWLINES => 0);
$LISTOFUNITSwriter->startTag('listOfUnits');
$LISTOFUNITSwriter->endTag('listOfUnits');

my $allStringsOfUnitsInModel='';
foreach my $key (keys %unitsModelInfo){
 my $stringOfUnitsInModel="\<unit";
 foreach my $key2 (keys $unitsModelInfo{$key}){
  $stringOfUnitsInModel=$stringOfUnitsInModel." $key2=\"$unitsModelInfo{$key}{$key2}\"";
 }
 $stringOfUnitsInModel=$stringOfUnitsInModel." \/\>\n";
 $allStringsOfUnitsInModel=$allStringsOfUnitsInModel.$stringOfUnitsInModel;
 #print "STRING UNIT: $stringOfUnitsInModel";
}
print "$allStringsOfUnitsInModel\n";
#TODO Allow the incorporation of new units?

#
my $MODELwriter = XML::Writer->new(OUTPUT => 'self', NEWLINES => 0);
$MODELwriter->startTag('model', %infoModelLine);
$MODELwriter->characters("Hello, world!");
$MODELwriter->endTag('model');

#
my $SBMLwriter = XML::Writer->new(OUTPUT => 'self', NEWLINES => 0);
$SBMLwriter->startTag('sbml', %infoSBMLLine);
my $ModelTagString=$MODELwriter->to_string;
$SBMLwriter->characters($ModelTagString);
$SBMLwriter->endTag('sbml');

#
my $finalResult=$SBMLwriter->to_string;
$finalResult =~ s/&lt;/</g;
$finalResult =~ s/&gt;/>/g;
print $finalResult."\n";

##############################################
#Usage
##############################################
sub usage{
    print STDERR "$0 version $version, Copyright (C) 2016 Renato Augusto Correa dos Santos, Diego Mauricio Riano Pachon\n";
    print STDERR "$0 comes with ABSOLUTELY NO WARRANTY; for details type `$0 -l'.\n";
    print STDERR "This is free software, and you are welcome to redistribute it under certain conditions;\n";
    print STDERR "type `$0 -l' for details.\n";
    print STDERR <<EOF;
NAME
    $0 Manipulates SBML file (Level 3 Version 1) reconstruction file using BiGG rules
    
USAGE
    $0 -i infile.xml --input_table table.txt
i
OPTIONS
    --infile                   -i       Input file (SBML3 format)                        REQUIRED
    --input_reaction_table     -t       Input table (separated by TAB)                   REQUIRED
                                        Lines in table must present information about metabolite and reaction association,
                                        and the following fields (in this order) are required:

                                        a) Reaction BiGG ID,
                                        b) Reaction name
                                        c) Reaction type (reaction or pseudoreaction)
                                        d) Reaction EC/TCDB code
                                        e) A boolean field whether it is a transport reaction or not
                                        f) A boolean field indicating whether the reaction is reversible or not
                                        g) Metabolite BiGG ID
                                        ih) Metabolite type in reaction (reactant or product)
                                        i) Metabolite stoichiometry coefficient in reaction
                                        j) Metabolite compartment
                                        k) Gene association: list of gene names in BiGG format, separated by comma
                                        m) Upper bound flux
                                        n) Lower bound flux

    --help                     -h       This help.
    --license                  -l       License.

    --k                                 Skip online validator (not implemented yet)

EOF
}

##############################################
#License
##############################################
sub license{
    print STDERR <<EOF;

Copyright (C) 2016 Renato Augusto Correa dos Santos, Diego Mauricio Riano Pachon
e-mail: renatoacsantos\@gmail.com, diriano\@gmail.com

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
