#!/usr/bin/perl

use warnings;
use strict;
use LWP::Simple;

my $genesKeggPriam = $ARGV[0];
my %keggReactionsFormulae;
my $reactionCount=1;

open(GENEKEGGPRIAM,$genesKeggPriam);

while(<GENEKEGGPRIAM>) {
 next if (/^TypeOfAssign/);
 chomp;
 my ($TypeOfAssign,$BooleanAND,$NumberOfKeggReactions,$KeggReactions,$GeneAndPriamProb)=split(/\t/,$_);
 my @KeggReactionIDs = split(",",$KeggReactions);
 next if ($NumberOfKeggReactions == 0);
 my @geneMatches = $GeneAndPriamProb =~ /(evm\.model\.KI\d+\.\d+\.\d+)/g;
 foreach my $keggID (@KeggReactionIDs) {
  unless($keggReactionsFormulae{$keggID}) {
   my $kegg_url = "http://rest.kegg.jp/find/reaction/$keggID";
   my $content = get($kegg_url);
   if($content =~ /rn:(R\d+)	(.+$)$/) {
    my @contentLines = split(/\n/,$content);
    if(scalar(@contentLines) > 1){
     die "Something wrong with reaction $keggID. There is more than one formulae for this reaction!\n";
    }
    foreach my $line (@contentLines){
     my ($reactionIDrest,$formulae)=split(/\t/,$line);
     $reactionIDrest =~ s/rn://g;
     if($reactionIDrest eq $keggID){
      $keggReactionsFormulae{$keggID}=$formulae;
     } else {
      die "Something wrong in reaction $keggID\n";
     }
    }
   }
   sleep 2;
  }
 }
 print "R$reactionCount\t";
 foreach my $kid (@KeggReactionIDs){
  print "$keggReactionsFormulae{$kid};;;";
 }
 print "\t";
 foreach my $kid (@KeggReactionIDs){
  print "$kid;;;";
 }
 print "\t";
 if($BooleanAND eq 'True'){
  print join(" and ",@geneMatches);
 } else {
  print join(" or ",@geneMatches);
 }
 print "\n";
 $reactionCount++;
}

close(GENEKEGGPRIAM);
