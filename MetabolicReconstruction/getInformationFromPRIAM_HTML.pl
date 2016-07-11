#!/usr/bin/perl

use warnings;
use strict;

my @infiles = <ma*.html>;

foreach my $infile (@infiles) {
 open(INFILE,$infile);
 while(<INFILE>) {
  chomp;
  if(/^<area.+onmouseover=\"showInfo\((.+)\);\"\/>/) {
   my $EC='';
   my $typeEC='';
   my @fields=split(/<br><br>/,$1);
   my $gene='';
   my @genesThisEC=();
   my $prob='';
   my @probabilitiesGenesThisEC=();

   # Get EC code
   if (/EC:(\d+\.\d+\.\d+\.\d+)<br><br>/ or /EC:(\d+\.\d+\.\d+\.\-)<br><br>/ or /EC:(\d+\.\d+\.\-\.\-)<br><br>/ or /EC:(\d+\.\-\.\-\.\-)<br><br>/) {
    my $tmpEC=$1;
    if(!$EC) {
     $EC = $tmpEC;
    } else {
     die "Something wrong with EC $tmpEC (EC $EC already exists)\n";
    }
   }

   # Enzyme predicted using PRIAM signature
   if(/This activity is predicted with a probability of (\d+\.\d+)/) {
    $typeEC="DETECTED_GENE";
    for my $f (@fields) {
     if($f =~ /^Predicted enzyme/){
      while($f =~ /(evm\.model\.KI\d+\.\d+\.\d+)( | \(incomplete\) )\(proba=(\d+\.\d+)\)/g) {
       $gene=$1;
       $prob=$3;
       push(@genesThisEC,$gene);
       push(@probabilitiesGenesThisEC,$prob);
      }
     }
    }
   }

   # Predictable, but not detected
   if(/None of the sequences of your proteome has matched signatures for that EC numberThus this activity is predicted as beeing absent from your proteome/){
    if($typeEC){ die "Weird: $EC already set as \'$typeEC\', but trying to set as \'ABSENT_IN_PROTEOME\'"; }
    if(/However, in the context of the network, this reaction might be present/) {
     $typeEC="ABSENT_IN_PROTEOME_GAPFILLED";
    } else {
     $typeEC="ABSENT_IN_PROTEOME"; # In this case there aren't probabilities nor proteins
    }
   }

   # Enzyme with no PRIAM signature; there aren't genes nor probs
   if(/This EC does not have any signature in PRIAM yet/) {
   if($typeEC){ die "Weird: $EC already set as \'$typeEC\', but trying to set as \'UNPREDICTABLE\'"; }
    if(/However, in the context of the network, this reaction is supposed to be present/){
     $typeEC='UNPREDICTABLE_GAPFILLED';
    } else {
     $typeEC='UNPREDICTABLE';
    }
   }

   # Enzyme is predictable, but not found. There is a gene and probability.
   if(/PRIAM predict this activity as beeing absent from your proteome/){
    if($typeEC){
     die "Weird: $EC already set as \'$typeEC\', but trying to set as \'FOUND_LOW_PROB\'";
    }
    if(/However, in the context of the network, this reaction is supposed to be present/) {
     $typeEC='FOUND_LOW_PROB_GAPFILLED';
    } else {
     $typeEC='FOUND_LOW_PROB';
    }
    for my $f (@fields) {
     if($f =~ /^Best enzyme candidate/){
      while($f =~ /(evm\.model\.KI\d+\.\d+\.\d+)( | \(incomplete\) )\(proba=(\d+\.\d+)\)/g) {
       $gene=$1;
       $prob=$3;
       push(@genesThisEC,$gene);
       push(@probabilitiesGenesThisEC,$prob);
      }
     }
    }
   }

   # Ambiguous cases
   if(/However, this EC number is associated with multiple possible reactions and this is not one of those retained in the automatically reconstructed metabolic network/){
    if($typeEC eq "DETECTED_GENE"){
     $typeEC=$typeEC."_AMBIGUOUS";
    } else {
     die "Weird: $EC already set as \'$typeEC\' (should be \'DETECTED_GENE\'), but trying to set as \'AMBIGUOUS\'";
    }
   }

   # Check whether gene and probability arrays are empty (length = zero)
   if(scalar(@genesThisEC) == 0) {
    push(@genesThisEC,"NO_GENES");
   }
   if(scalar(@probabilitiesGenesThisEC) == 0) {
    push(@probabilitiesGenesThisEC,"NO_PROBS");
   }

   # Modify variable storing html file name and print result for that enzyme
   $infile =~ s/\.html//g;
   print "$infile\t$EC\t$typeEC\t".join(",",@genesThisEC)."\t".join(",",@probabilitiesGenesThisEC)."\n";

  }
 }
 close(INFILE);
}
