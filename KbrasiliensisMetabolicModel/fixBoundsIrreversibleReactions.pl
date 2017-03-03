#!/usr/bin/perl

use warnings;
use strict;

my $sbmlModel=$ARGV[0];

open(SBMLM,$sbmlModel);

while(<SBMLM>){
 chomp;
 if(/^<reaction .+ reversible="(\S+)"/){
  my $newLine='';
  my $revers=$1;
  #print "$1\n";
  if($revers eq 'false'){
   $newLine = $_;
   $newLine =~ s/fbc:lowerFluxBound=\"cobra_default_lb\"/fbc:lowerFluxBound=\"0\"/g;
   print "$newLine\n";
  } else {
   print "$_\n";
  }
 }
 else {
  print "$_\n";
 }
}

close(SBMLM);
