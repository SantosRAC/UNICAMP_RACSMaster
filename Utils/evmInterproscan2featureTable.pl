#!/usr/bin/perl

use warnings;
use strict;

my $evmGffFile = $ARGV[0];
my $interproScanFile = $ARGV[1];

# Information InterProScan

my %interProScanAnnotation;

open(INTERPROSCAN,$interproScanFile);

while(<INTERPROSCAN>){
 chomp;
 #evm.model.KI545862.1.400        c1932f5ac3f6f4274218e5d6b2427aa4        728     Phobius TRANSMEMBRANE   Region of a membrane-bound protein predicted to be embedded in the membrane.    367     384     -       T       07-04-2016
 my @fields=split(/\t/,$_);
 # Gene identifiers have 'TU', not 'model'
 my($gene,undef,undef,undef,undef,$desc,undef,undef,undef,undef,undef)=split(/\t/,$_);
 $gene =~ s/\.model\./\.TU\./g;
 unless ($desc eq /^$/){
  if($interProScanAnnotation{$gene}){
   unless($desc ~~ @{$interProScanAnnotation{$gene}}){
    push(@{$interProScanAnnotation{$gene}},$desc);
   }
  }else{
   @{$interProScanAnnotation{$gene}}=($desc);
  }
 }
}

close(INTERPROSCAN);

# Information GFF3

#KI545851.1      EVM     gene    695     2101    .       -       .       ID=evm.TU.KI545851.1.1;Name=EVM_prediction_KI545851.1.1
#KI545851.1      EVM     mRNA    695     2101    .       -       .       ID=evm.model.KI545851.1.1;Parent=evm.TU.KI545851.1.1;Name=EVM_prediction_KI545851.1.1
#KI545851.1      EVM     exon    695     2101    .       -       .       ID=evm.model.KI545851.1.1.exon1;Parent=evm.model.KI545851.1.1
#KI545851.1      EVM     CDS     695     2101    .       -       0       ID=cds.evm.model.KI545851.1.1;Parent=evm.model.KI545851.1.1

my @sequences=();
my @features=();
my %featuresInfo;
my %CDScount;
my %featuresSeqs;

open(EVMGFFFILE,$evmGffFile);

while(<EVMGFFFILE>){
 chomp;
 next if(/^#/);
 my ($seq,$source,$feattype,$init_pos,$end_pos,undef,$strand,$codon_start,$additionalInfo)=split(/\t/,$_);
 die unless(scalar(split(/\t/,$_)) == 9);
 unless($seq ~~ @sequences){
  push(@sequences,$seq);
 }
 my @addInfoFields=split(/;/,$additionalInfo);

 my $featIdentifier='';
 my $featIdentParent='';
 my $countMatchID=0;
 my $countMatchParent=0; 

 foreach my $addInfoField (@addInfoFields){
  if($addInfoField =~ /ID=/){
   $countMatchID++;
   $featIdentifier=$addInfoField;
   if($countMatchID == 1){
    $featIdentifier =~ s/ID=//g;
    if(($featIdentifier ~~ @features) and ($feattype ne 'CDS')){
     die "Something wrong: $featIdentifier\n";
    } elsif($feattype eq 'CDS') {
     if($featuresInfo{$featIdentifier}){
      $CDScount{$featIdentifier}++;
      $featuresInfo{$featIdentifier}{'CDS'}{$CDScount{$featIdentifier}}{'codon_start'}=$codon_start+1;
      $featuresInfo{$featIdentifier}{'CDS'}{$CDScount{$featIdentifier}}{'init'}=$init_pos;
      $featuresInfo{$featIdentifier}{'CDS'}{$CDScount{$featIdentifier}}{'end'}=$end_pos;
      $featuresInfo{$featIdentifier}{'CDS'}{$CDScount{$featIdentifier}}{'strand'}=$strand;
     }else{
      push(@features,$featIdentifier);
      $CDScount{$featIdentifier}=1;
      $featuresInfo{$featIdentifier}{'CDS'}{1}{'codon_start'}=$codon_start+1;
      $featuresInfo{$featIdentifier}{'CDS'}{1}{'init'}=$init_pos;
      $featuresInfo{$featIdentifier}{'CDS'}{1}{'end'}=$end_pos;
      $featuresInfo{$featIdentifier}{'CDS'}{1}{'strand'}=$strand;
      $featuresInfo{$featIdentifier}{'feattype'}=$feattype;
     }
    }else{
     push(@features,$featIdentifier);
     $featuresInfo{$featIdentifier}{'init'}=$init_pos;
     $featuresInfo{$featIdentifier}{'end'}=$end_pos;
     $featuresInfo{$featIdentifier}{'feattype'}=$feattype;
     $featuresInfo{$featIdentifier}{'strand'}=$strand;
    }
   } else {
    die "Something wrong\n";
   }
  }
  if($addInfoField =~ /Parent=/){
   unless($feattype eq 'gene'){
    if($feattype eq 'exon'){
     $featIdentParent=$addInfoField;
     $featIdentParent =~ s/Parent=//g;
     $featuresInfo{$featIdentifier}{'parent'}=$featIdentParent;
    } else {
     $featIdentParent=$addInfoField;
     $featIdentParent =~ s/Parent=//g;
     $featuresInfo{$featIdentifier}{'parent'}=$featIdentParent;
    }
   }
  }
 }

 if($featuresSeqs{$seq}){
  unless($featIdentifier ~~ @{$featuresSeqs{$seq}}){
   push(@{$featuresSeqs{$seq}},$featIdentifier);
  }
 } else {
  @{$featuresSeqs{$seq}}=($featIdentifier);
 }

}

close(EVMGFFFILE);

my %parent2exon;

foreach my $seq (@sequences){
 print ">$seq\n";
 foreach my $feat (@{$featuresSeqs{$seq}}){
  if($featuresInfo{$feat}{'feattype'} eq 'CDS'){
   foreach my $cdsNum (keys $featuresInfo{$feat}{'CDS'}){
    if($featuresInfo{$feat}{'CDS'}{$cdsNum}{'strand'} eq '-'){
     print "$featuresInfo{$feat}{'CDS'}{$cdsNum}{'end'}\t$featuresInfo{$feat}{'CDS'}{$cdsNum}{'init'}\t$featuresInfo{$feat}{'feattype'}\n";
    } else {
     print "$featuresInfo{$feat}{'CDS'}{$cdsNum}{'init'}\t$featuresInfo{$feat}{'CDS'}{$cdsNum}{'end'}\t$featuresInfo{$feat}{'feattype'}\n";
    }
    print "			codon_start	$featuresInfo{$feat}{'CDS'}{$cdsNum}{'codon_start'}\n";
    my $feat2=$feat;
    $feat2 =~ s/cds\.//g;
    print "			gene	$featuresInfo{$feat2}{'parent'}\n";
    if($interProScanAnnotation{$featuresInfo{$feat2}{'parent'}}){
     foreach my $note (@{$interProScanAnnotation{$featuresInfo{$feat2}{'parent'}}}){
      print "			product	$note\n";
     }
    } else {
     print "			note	hypothetical protein\n";
    }
   }
  } else {
   if ($featuresInfo{$feat}{'strand'} eq '-'){
    print "$featuresInfo{$feat}{'end'}\t$featuresInfo{$feat}{'init'}\t$featuresInfo{$feat}{'feattype'}\n";
   } else {
    print "$featuresInfo{$feat}{'init'}\t$featuresInfo{$feat}{'end'}\t$featuresInfo{$feat}{'feattype'}\n";
   }
   if($featuresInfo{$feat}{'feattype'} eq 'exon'){
    my $feat2=$featuresInfo{$feat}{'parent'};
    print "			gene	$featuresInfo{$feat2}{'parent'}\n";
    if($interProScanAnnotation{$featuresInfo{$feat2}{'parent'}}){
     foreach my $note (@{$interProScanAnnotation{$featuresInfo{$feat2}{'parent'}}}){
      print "			product $note\n";
     }
    } else {
     print "			note	hypothetical protein\n";
    }
    if($parent2exon{$featuresInfo{$feat}{'parent'}}){
     $parent2exon{$featuresInfo{$feat}{'parent'}}++;
     print "			number	$parent2exon{$featuresInfo{$feat}{'parent'}}\n";
    } else {
     $parent2exon{$featuresInfo{$feat}{'parent'}}=1;
     print "			number	$parent2exon{$featuresInfo{$feat}{'parent'}}\n";
    }
   } elsif($featuresInfo{$feat}{'feattype'} eq 'gene'){
    print "			gene	$feat\n";

   } elsif($featuresInfo{$feat}{'feattype'} eq 'mRNA') {

   } else {
    die "Something wrong...\n";
   }
  }
 }
}

