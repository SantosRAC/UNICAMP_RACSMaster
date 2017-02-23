#!/usr/bin/perl

use warnings;
use strict;

my $mappingTableCompoundFormulae=$ARGV[0];
my $sbtabCompoundTable=$ARGV[1];

my %metabolite2formula;

open(MAPFILE,$mappingTableCompoundFormulae);

while(<MAPFILE>){
 chomp;
 my ($metaboliteID,$metaboliteFormula,undef)=split(/\t/,$_);
 if($metabolite2formula{$metaboliteID}){
  die "There is a previous formula included in hash for this metabolite: $metaboliteID\n";
 } else {
  $metabolite2formula{$metaboliteID}=$metaboliteFormula;
 }
}

close(MAPFILE);

open(SBTABCT,$sbtabCompoundTable);

## Expected header of input table
##!!SBtab SBtabVersion=’1.0’      TableName=’K. brasiliensis SBTab Compound Table’        TableType=’Compound’
##!ID     !Name   !Identifiers:kegg.compound      !Identifiers:bigg.metabolite    !Identifiers:biocyc

print "!!SBtab SBtabVersion=’1.0’\tTableName=’K. brasiliensis SBTab Compound Table’\tTableType=’Compound’\n!ID\t!Name\t!Identifiers:kegg.compound\t!Identifiers:bigg.metabolite\t!Identifiers:biocyc\t!StructureFormula\n";

while(<SBTABCT>){
 chomp;
 next if (/SBtabVersion/);
 next if (/Identifiers:bigg.metabolite/);
 my ($ID,$Name,$IdentifiersKeggCompound,$IdentifiersBiggMetabolite,$IdentifiersBiocyc)=split(/\t/,$_);
 if (($IdentifiersKeggCompound) and ($IdentifiersBiocyc)) {
  die "Both KEGG and BioCyc compound identifiers are present. Check line with $IdentifiersBiocyc\n";
 } elsif ($IdentifiersKeggCompound) {
  if($metabolite2formula{$IdentifiersKeggCompound}){
   print "$ID\t$Name\t$IdentifiersKeggCompound\t$IdentifiersBiggMetabolite\t\t$metabolite2formula{$IdentifiersKeggCompound}\n";
  } else {
   print "$ID\t$Name\t$IdentifiersKeggCompound\t$IdentifiersBiggMetabolite\t\tNO_FORMULA_AVAILABLE\n";
  }
 } elsif ($IdentifiersBiocyc) {
  $IdentifiersBiocyc =~ s/-/__45__/g;
  if ($metabolite2formula{$IdentifiersBiocyc}) {
   print "$ID\t$Name\t\t$IdentifiersBiggMetabolite\t$IdentifiersBiocyc\t$metabolite2formula{$IdentifiersBiocyc}\n";
  } else {
   print "$ID\t$Name\t\t$IdentifiersBiggMetabolite\t$IdentifiersBiocyc\tNO_FORMULA_AVAILABLE\n";
  }
 } else {
  die "Neither of the identifiers (KEGG or BioCyc) is present in line:\n$_\n";
 }
}

close(SBTABCT);
