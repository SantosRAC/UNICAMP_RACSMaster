#!/usr/bin/perl

use warnings;
use strict;
use LWP::Simple;

# Get a list of KEGG Compounds and get their corresponding names

my $KeggCompoundFile = $ARGV[0];
my %keggCompoundsNames;

open(KEGGCOMPIDS,$KeggCompoundFile);

while(<KEGGCOMPIDS>) {
 chomp;
 unless(/^C\d\d\d\d\d$/) { die "KEGG identifier does not appear as expected\n"; }
 my $keggCompID=$_;
 my $kegg_url = "http://rest.kegg.jp/find/compound/$keggCompID";
 my $content = get($kegg_url);
 if($content =~ /cpd:(C\d\d\d\d\d)	(.+)$/) {
  my @contentLines = split(/\n/,$content);
  if(scalar(@contentLines) > 1){
   die "Something wrong with compound $keggCompID. There is more than one line of names for this compound!\n";
  }
  foreach my $line (@contentLines){
   my ($compoundIDrest,$compNames)=split(/\t/,$line);
   my @compoundNames=split(/;/,$compNames);
   $compoundIDrest =~ s/cpd://g;
   if($keggCompoundsNames{$compoundIDrest}) {
    die "$compoundIDrest already exist in hash!\n";
   } else {
    $keggCompoundsNames{$compoundIDrest}=$compoundNames[0];
    print "$compoundIDrest\t$keggCompoundsNames{$compoundIDrest}\n";
   }
  }
 }
 sleep(2);
}

close(KEGGCOMPIDS);
