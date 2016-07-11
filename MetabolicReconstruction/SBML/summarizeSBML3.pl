#!/usr/bin/perl
#TODO: Add validation of SBML model using Online validator API: http://sbml.org/Facilities/Validator/ and LWP (Perl)
require v5.10.1;
use strict;
use warnings;
use diagnostics;
use LWP::UserAgent;
use Getopt::Long;
use XML::Twig;

my $version='0.1';
my $sbmlLevel='';
my $sbmlVersion='';
my $infile='';
my $outfile='';
my $license='';
my $help='';
my $debug=0;

######################
# Variables used for summarizing SBML elements 
######################

my @compartmentsInModel=();
my @reactantsInReactions=();
my @productsInReactions=();
my @speciesInModel=();
my %speciesModelInfo;
my @reactionsInModel=();
my %reactionsModelInfo;
my $summary='';
my $summary_species='';
my $summary_reaction='';
my $ua = LWP::UserAgent->new;

#####################
# Optional variables
#####################
my $valitateBoolean='';

GetOptions(
    'license|l'    => \$license,
    'help|h|?'     => \$help,
    'infile|i=s'   => \$infile,

# Getting summary statistics from model
    'summary'                   =>  \$summary,
    'summary_species=s'         =>  \$summary_species,
    'summary_reaction=s'        =>  \$summary_reaction,

# Optionals
    'v'                         =>  \$valitateBoolean,

# Debug option (for developers)
    'debug|d:i'    => \$debug,
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
if(!$summary_species && !$summary_reaction){
 $summary=1;
}

unless ($valitateBoolean) {
 
 $ua->agent("CTBE-summarizeSBML3/0.1");
 $ua->from("bce.ctbe\@gmail.com");
 $ua->env_proxy;

 my $htmlValidateFile = getHTML($infile);
 my $twigValidator=XML::Twig->new();
 $twigValidator->parsefile($htmlValidateFile);
 my $rootValitator = $twigValidator->root;
 #print $rootValitator->name."\n";

}

parseXML();

sub parseXML {

 my $twig=XML::Twig->new();
 $twig->parsefile($infile) or die "cannot parse $infile";
 my $root = $twig->root;
 $sbmlLevel = $root->att('level');
 $sbmlVersion = $root->att('version');

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
 
  # Parsing compartments
  my @compartmentLists = $mod->children('listOfCompartments');
  foreach my $compl (@compartmentLists) {
   my @AllCompartments = $compl->children('compartment');
   foreach my $comp (@AllCompartments) {
   my $compId='';
   my $compName='';
   $compId = $comp->att('id');
   $compName = $comp->att('name');
   push(@compartmentsInModel,$compId);
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
    my $modSpeciesName = $spemodel->att('name') ? $spemodel->att('name') : $spemodel->att('id');
    my $modSpeciesCompartment = $spemodel->att('compartment');
    my $modeSpeciesMetaID = $spemodel->att('metaid');
    my $annotResource='';
    my $Notes_COMPOUND_ID='';

    next if($modSpeciesID ~~ @speciesInModel);
    push(@speciesInModel,$modSpeciesID);

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

    # Parsing notes for a given species in the model
    my @spNotes = $spemodel->children('notes');
    foreach my $notes (@spNotes) {
     my $COMPOUND_ID_HTML = $notes->first_child();
     my $COMPOUND_ID_TEXT = $COMPOUND_ID_HTML->text;
#     $Notes_COMPOUND_ID = $COMPOUND_ID_TEXT =~ s/COMPOUND_ID: //r;
     $Notes_COMPOUND_ID = $COMPOUND_ID_TEXT ;
     $Notes_COMPOUND_ID =~ s/COMPOUND_ID: //;
    }

    # Adding information about model species in hash
    $speciesModelInfo{$modSpeciesID}{'overall'}{'id'}=$modSpeciesID;
    $speciesModelInfo{$modSpeciesID}{'overall'}{'metaid'}=$modeSpeciesMetaID;
    $speciesModelInfo{$modSpeciesID}{'overall'}{'name'}=$modSpeciesName;
    $speciesModelInfo{$modSpeciesID}{'overall'}{'compartment'}=$modSpeciesCompartment;
    $speciesModelInfo{$modSpeciesID}{'notes'}{'compound_id'}=$Notes_COMPOUND_ID;
    $speciesModelInfo{$modSpeciesID}{'annotation'}{'resource'}=$annotResource;
   }
  }

  # Parsing reactions
  my @reactionLists = $mod->children('listOfReactions');
  foreach my $real (@reactionLists) {

   my @AllReactions = $real->children('reaction');
   foreach my $rea (@AllReactions) {
    my $reactID='';
    my $reactMetaID='';
    my $reactAnnotResource='';
    my @reactantsReaction=();
    my @productsReaction=();

    $reactID=$rea->att('id');
    $reactMetaID=$rea->att('metaid');
    push(@reactionsInModel,$reactID);

    # Parsing list of reactants of a given reaction
    my @reactantList = $rea->children('listOfReactants');
    foreach my $reactl (@reactantList) {
     my $spID='';
     my @AllReactants = $reactl->children('speciesReference');
     foreach my $reaspe (@AllReactants) {
      $spID = $reaspe->att('species');
      #print "$spID\n"; # print observed unsorted reactants
      push(@reactantsReaction,$spID);
      next if ($spID ~~ @reactantsInReactions);
      push(@reactantsInReactions,$spID);
     }
    }

    # Parsing list of products of a given reaction
    my @productList = $rea->children('listOfProducts');
    foreach my $productl (@productList) {
     my $spID='';
     my @AllProducts = $productl->children('speciesReference');
     foreach my $productspe (@AllProducts) {
      $spID = $productspe->att('species');
      #print "$spID\n"; # print observed unsorted products
      push(@productsReaction,$spID);
      next if ($spID ~~ @productsInReactions);
      push(@productsInReactions,$spID);
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
    $reactionsModelInfo{$reactID}{'overall'}{'id'}=$reactID;
    $reactionsModelInfo{$reactID}{'overall'}{'metaid'}=$reactMetaID;
    $reactionsModelInfo{$reactID}{'overall'}{'resource'}=$reactAnnotResource;
    @{ $reactionsModelInfo{$reactID}{'reaction'}{'listOfReactants'} }=@reactantsReaction;
    @{ $reactionsModelInfo{$reactID}{'reaction'}{'listOfProducts'} }=@productsReaction;

   }
  }
 }

 # Provide summary statistics
 if ($summary) {
  print "######### SUMMARY STATISTICS OF THE METABOLIC RECONSTRUCTION #########
There is/are ".scalar(@compartmentsInModel)." involved compartment(s) in the model
There is/are ".scalar(@reactionsInModel)." involved reaction(s) in the model
There is/are ".scalar(@speciesInModel)." involved metabolite(s)/species in the model
There is/are ".scalar(@reactantsInReactions)." involved reactants in reactions of the model
There is/are ".scalar(@productsInReactions)." involved products in reactions of the model\n";
 }

 # Provide summary statistics for a given species
 if ($summary_species) {
  if ($speciesModelInfo{$summary_species}{'overall'}{'id'}) {
   if ($speciesModelInfo{$summary_species}{'overall'}{'id'} =~ /(G\d{5})/) {
    print "######### SUMMARY STATISTICS FOR COMPOUND $summary_species #########
    Identifier: $speciesModelInfo{$summary_species}{'overall'}{'id'},
    Meta ID: $speciesModelInfo{$summary_species}{'overall'}{'metaid'},
    Compartment: $speciesModelInfo{$summary_species}{'overall'}{'compartment'},
    Compound Identifier (added as note in the model): $speciesModelInfo{$summary_species}{'notes'}{'compound_id'},
    Resource: $speciesModelInfo{$summary_species}{'annotation'}{'resource'}\n"
   } else {
    print "######### SUMMARY STATISTICS FOR COMPOUND $summary_species #########
    Identifier: $speciesModelInfo{$summary_species}{'overall'}{'id'},
    Meta ID: $speciesModelInfo{$summary_species}{'overall'}{'metaid'},
    Name: $speciesModelInfo{$summary_species}{'overall'}{'name'},
    Compartment: $speciesModelInfo{$summary_species}{'overall'}{'compartment'},
    Compound Identifier (added as note in the model): $speciesModelInfo{$summary_species}{'notes'}{'compound_id'},
    Resource: $speciesModelInfo{$summary_species}{'annotation'}{'resource'}\n"
   }
  } else {
   die "Compound $summary_species does not exist in the model or no identifier information is available\n";
  }
 }

 if ($summary_reaction) {
  print "######### SUMMARY STATISTICS FOR COMPOUND $summary_reaction #########
    Identifier: $reactionsModelInfo{$summary_reaction}{'overall'}{'id'}
    Meta ID: $reactionsModelInfo{$summary_reaction}{'overall'}{'metaid'}
    Resource: $reactionsModelInfo{$summary_reaction}{'overall'}{'resource'}
    Reactants: ".join(',',@{ $reactionsModelInfo{$summary_reaction}{'reaction'}{'listOfReactants'} })."
    Products: ".join(',',@{ $reactionsModelInfo{$summary_reaction}{'reaction'}{'listOfProducts'} })."\n";
 }

}

sub getHTML{
 my $file = shift;
 
 return $html;
}

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
    $0 Summarizes SBML file Level 3 Version 1 reconstruction file with KEGG identifiers
    
USAGE
    $0 -i infile.xml (same as using --summary)

    $0 -i infile.xml --summary

    $0 --infile network.sbml --summary_reaction R02108

    $0 --infile network.sbml --summary_species C00364

    $0 -i infile.xml --add_reaction XXXXXX --reactants_new_reaction XXXXX --products_new_reaction XXXXX #TODO

OPTIONS
    --infile          -i       Input file (XML)                                                         REQUIRED
    --summary                  Shows overall statistics (reactions, metabolites, compartments, etc.)    OPTIONAL
    --summary_species          Shows overall details for a given species in the model                   OPTIONAL
    --summary_reaction         Shows overall details for a given reaction in the model                  OPTIONAL
    --help            -h       This help.
    --license         -l       License.

    --k                        Skip online validator

ADVANCED OPTIONS (developers) #TODO
    --debug           -d    debug (INT).

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
