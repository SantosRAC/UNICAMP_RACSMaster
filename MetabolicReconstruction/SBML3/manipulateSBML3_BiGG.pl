#!/usr/bin/perl
#TODO: Add validation of SBML model using Online validator API: http://sbml.org/Facilities/Validator/ and LWP (Perl)
require v5.10.1;
use strict;
use warnings;
use diagnostics;
#use LWP::UserAgent; # This will be used when the Perl code to connect to SBML validator is available
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
my $reactionTable='';
#my $ua = LWP::UserAgent->new; # This will be used only when the online validator is available

##############################################
#Variables related to the novel reactions
##############################################
my %reaction2associatedLinks;
my %reaction2name;
my %reaction2id;
my %reaction2associatedGenes;
my %reaction2reversibility;
my %reaction2type;
my %reaction2products;
my %reaction2reactants;
my %reaction2compartment;
my %reaction2EC;
my @allowedCompartments=('UNK_COMP','cytoplasm','mitochondrion','peroxisome','Golgi','vacuole','ER','plasma_membrane','nucleus','extracellular');

##############################################
#Variables related to info existent in model
##############################################
my @compartmentsInModel=();

GetOptions(
    'license|l'       => \$license,
    'help|h|?'        => \$help,
    'infile|i=s'      => \$infile,
    'input_table|t=s' => \$reactionTable,
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

checkInputTable();
parseXML();

##############################################
#Check input table with information about
#reactions and metabolites
##############################################

sub checkInputTable {

 open(INPUTTABLE,$reactionTable);

 while(<INPUTTABLE>){
  chomp;
  next if (/^#/);
  my ($reactionBiGG,$reactionType,$reversible,$metaboliteBiGG,$metaboliteType,$stoichiometry,$compartment) = split(/\t/,$_);
  # Check if reaction ID is a BiGG identifier
  if($reactionBiGG =~ /^R\d{5}$/) {
   $reactionType = lc($reactionType);
   # Check if reaction type is allowed (conversion or transport)
   if(($reactionType =~ /conversion/) or ($reactionType =~ /transport/)){
    # Check if reaction type was previously set. If it exists, check if next lines for the reaction have the same type.
    if ($reaction2type{$reactionBiGG}){
     if($reaction2type{$reactionBiGG} eq $reactionType) {
      #print "$reactionBiGG\t$reactionType (check reaction type)\n";
     } else {
      die "FATAL: Ambiguous reaction type for reaction: $reactionBiGG\n";
     }
    } else {
     $reaction2type{$reactionBiGG}=$reactionType;
    }
    $reversible = lc($reversible);
    if(($reversible =~ /true/) or ($reversible =~ /false/)){
     # Check if reaction reversibility was previously set. If it exists, check if next lines for the reaction have the same reversibility.
     if($reaction2reversibility{$reactionBiGG}){
      if($reaction2reversibility{$reactionBiGG} eq $reversible){
       #print "$reactionBiGG\t$reactionType\t$reversible (check reaction reversibility)\n";
      } else{
       die "FATAL: Ambiguous reaction reversibility for reaction: $reactionBiGG\n";
      }
     } else {
      $reaction2reversibility{$reactionBiGG}=$reversible;
     }
    } else {
     die "FATAL: The reaction reversibility must be 'true' or 'false'\n";
    }
    # Check metabolite ID (must be a BiGG identifier)
    if($metaboliteBiGG =~ /C\d{5}/) {
     #print "$reactionBiGG\t$reactionType\t$reversible\t$metaboliteBiGG(Check metabolite for BiGG reactions)\n";
     $metaboliteType=lc($metaboliteType);
     if(($metaboliteType eq "product") or ($metaboliteType eq "reactant")){
      # Addes metabolite
      if($metaboliteType eq "product"){
       push(@{$reaction2products{$reactionBiGG}},$metaboliteBiGG);
      }
      if($metaboliteType eq "reactant"){
       push(@{$reaction2reactants{$reactionBiGG}},$metaboliteBiGG);
      }
      if ($compartment ~~ @allowedCompartments) {
       if($reaction2compartment{$reactionBiGG}) {
        #print "";
       } else {
        $reaction2compartment{$reactionBiGG}=$compartment;
       }
      } else {
       die "FATAL: compartment ($compartment) for metabolite $metaboliteBiGG in reaction $reactionBiGG is not allowed.\nOptions are: UNK_COMP (if compartment is unknown), cytoplasm, mitochondrion, peroxisome, Golgi, vacuole, ER, plasma_membrane, nucleus, extracellular\n";
      }
     } else {
      die "FATAL: a metabolite in your table does not have an allowed type (product or reactant)\n";
     }
    } else {
     die "FATAL: a metabolite in your table does not satisfy BiGG rules for compound ID.\n Check reaction $reactionBiGG, metabolite $metaboliteBiGG\n";
    }
   } else {
    die "FATAL: You must provide an allowed reaction type: 'Transport' or 'Conversion'\n";
   }
  } elsif ($reactionBiGG =~ /^BIOMASS$/) {
   #TODO Include case in which the BIOMASS reaction is being added
  } else {
    die "FATAL: Reaction type is not as expected ($reactionBiGG).\nIt must be a reaction in BiGG REACTION (e.g. R02108) or 'BIOMASS' (biomass pseudo-reaction)\n";
  }
 }

 #Check if a reaction has product but no reactant(s)
 foreach my $rea (keys %reaction2products){
  if(!$reaction2reactants{$rea}) {
   die "$rea has products but does not have reactants!\n";
  }
 }

 #Check if a reaction has reactant but no product(s)
 foreach my $rea (keys %reaction2reactants){
  if(!$reaction2products{$rea}) {
   die "$rea has reactants but does not have products!\n";
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
 
 my @rootChildren=$root->children('model');

 foreach my $rootChild (@rootChildren) {
  my @modelChildren=$rootChild->children();
  foreach my $modChild (@modelChildren) {

   #Get species in the model
   my @internModelChildrenSpe=$modChild->children('species');
   foreach my $intModChildSpe (@internModelChildrenSpe) {
    my $idCheckSpeciesInModel=$intModChildSpe->att('id');
    #Add species informed in the input table to the model
   }

   #Get compartments in the model
   my @internModelChildrenComp=$modChild->children('compartment');
   foreach my $intModChildComp (@internModelChildrenComp) {
    my $idCheckCompInModel=$intModChildComp->att('id');
    if($idCheckCompInModel ~~ @compartmentsInModel){
     #do nothing!
    } else {
     push(@compartmentsInModel,$idCheckCompInModel);
    }
   }
   #Check if compartment being added (for a given reaction) was existent in the model;
   #The script must add the compartment to the model if is doesn't exist.
   foreach my $comp (values %reaction2compartment){
    if($comp ~~ @compartmentsInModel){
     #print "Compartment $comp is in the model\n";
    } else {
     #Include here code to add compartment to the model
     print "Compartment \"$comp\" is not in the model!\n";
    }
   }

   #Get reactions in the model
   my @internModelChildrenRea=$modChild->children('reaction');
   foreach my $intModChildRea (@internModelChildrenRea) {
    my $idCheckReactionInModel=$intModChildRea->att('id');
    # Check if reaction being added was existent in the model; die if it already exists
    foreach my $rea (keys %reaction2type){
     if($rea eq $idCheckReactionInModel) {die "Reaction ($rea) already exists in model\n";}
    }
   }
  }
 }

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
    $0 Manipulates SBML file (Level 3 Version 1) reconstruction file using BiGG rules
    
USAGE
    $0 -i infile.xml --input_table table.txt

OPTIONS
    --infile          -i       Input file (XML,SBML3)                           REQUIRED
    --input_table     -t       Input table (separated by TAB)                   REQUIRED
                               Lines in table must present information about metabolite and reaction association,
                               and the following fields (in this order) are required:

                               a) Reaction (BiGG identifier),
                               b) Type of reaction (transport or conversion)
                               c) A boolean field indicating whether the reaction is reversible (yes/no)
                               d) Metabolite/ compound/ species BiGG identifier
                               e) Type of metabolite/ compound/ species (reactant or product)
                               f) Species Stoichiometry
                               g) Metabolite Compartment

                               Notes:
                               1) If the compounds involved in reaction are present in different compartments,
                               the reaction is assumed to be of the type 'transport'
                               2) If the reaction is of type 'BIOMASS' #TODO

    --help            -h       This help.
    --license         -l       License.

    --k                        Skip online validator (not implemented yet)

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
