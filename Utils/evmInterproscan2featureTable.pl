#!/usr/bin/perl

use warnings;
use strict;

my $evmGffFile = $ARGV[0];
my $interproScanFile = $ARGV[1];
my $tblOut = $ARGV[2];

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
 $gene =~ s/\./\_/g;

 # Replace any problematic description (reported in disc.report)
 # tbl2ans run: ~/Software/tbl2asn/linux64.tbl2asn -p . -j "[organism=Kalmanozyma brasiliensis] [strain=GHG001]" -M n -y "Comment" -i GCA_000497045.1_PSEUBRA1_genomic.fna -Z disc.report -t template.sbt -V b

 #FATAL: Remove organism from product name
 if($desc =~ /Animal heme peroxidase superfamily profile/){
  print "$gene: Replaced \'Animal heme peroxidase superfamily profile\' by \'Heme peroxidase superfamily\' in:\n$desc\n";
  $desc =~ s/Animal heme peroxidase superfamily profile/Heme peroxidase superfamily/g;
 }
 if($desc =~ /Haem peroxidase, animal/){
  print "$gene: Replaced \'Haem peroxidase, animal\' by \'Haem peroxidase\' in:\n$desc\n";
  $desc =~ s/Haem peroxidase, animal/Haem peroxidase/g;
 }
 if($desc =~ /Domain found in NIK1-like kinases, mouse citron and yeast ROM1, ROM2/){
  print "$gene: Replaced \'Domain found in NIK1-like kinases, mouse citron and yeast ROM1, ROM2\' by \'hypothetical protein\' in:\n$desc\n";
  $desc =~ s/Domain found in NIK1-like kinases, mouse citron and yeast ROM1, ROM2/hypothetical protein/g;
 }
 if($desc =~ /Putative DNA-binding domain in centromere protein B, mouse jerky and transposases\./){
  print "$gene: Replaced \'Putative DNA-binding domain in centromere protein B, mouse jerky and transposases.\' by \'hypothetical protein\' in:\n$desc\n";
  $desc =~ s/Putative DNA-binding domain in centromere protein B, mouse jerky and transposases\./hypothetical protein/g;
 }
 if(($desc =~ /(Staphylococcal nuclease homologues)/) or ($desc =~ /(Staphylococcal nuclease homologue)/)){
  print "$gene: Replaced \'$1\' by \'Nuclease\' in:\n$desc\n";
  $desc = 'Nuclease';
 }

 #FATAL: Possible parsing error or incorrect formatting; remove inappropriate symbols
 # Features starting with '
 if(lc($desc) =~ /homeobox/){
  print "$gene: Renamed \'$desc\' to \'Homeobox\'\n";
  $desc = 'Homeobox';
 }
 # Features ends with '.'
 if($desc =~ /\.$/){
  print "$gene: Removed \'\.\' at the end of the feature in\n$desc\n";
  $desc =~ s/\.$//g;
 }
 # Features with '. '
 if($desc =~ /\. /){
  print "$gene: Found '\. '. Replaced \'$desc\' by 'hypothetical protein'\n";
  $desc = 'hypothetical protein';
 }
 # features contains '@'
 if($desc =~ /\@/){
  print "$gene: Found '\@'. Replaced \'$desc\' by 'hypothetical protein'\n";
  $desc = 'hypothetical protein';
 }
 #features starts with 'hypothetical protein' but not equals 'hypothetical protein'
 if(lc($desc) =~ /hypothetical protein/){
  print "$gene: Replaced \'$desc\' by \'hypothetical protein\' in:\n$desc\n";
  $desc = 'hypothetical protein';
 }

 # Use American spelling
 if($desc =~ /organisation/){
  print "$gene: Replaced \'organisation\' by \'organization\' in:\n$desc\n";
  $desc =~ s/organisation/organization/g;
 }
 if($desc =~ /dimerisation/){
  print "$gene: Replaced \'dimerisation\' by \'dimerization\' in:\n$desc\n";
  $desc =~ s/dimerisation/dimerization/g;
 }
 if($desc =~ /characteris/){
  print "$gene: Replaced \'characteris\' by \'characteriz\' in:\n$desc\n";
  $desc =~ s/characteris/characteriz/g;
 }
 if($desc =~ /disulphide/){
  print "$gene: Replaced \'disulphide\' by \'disulfide\' in:\n$desc\n";
  $desc =~ s/disulphide/disulfide/g;
 }

 # Better 'hypothetical protein' than 'conserved protein'
 if(lc($desc) =~ /conserved protein/){
  print "$gene: Renamed \'$desc\' to \'hypothetical protein\'\n";
  $desc='hypothetical protein';
 }
 if($desc =~ /Region of a membrane-bound protein/){
  print "$gene: Renamed \'$desc\' to \'Membrane-bound protein\'\n";
  $desc='Membrane-bound protein';
 }

 if($desc =~ /domain of unknown function/){
  print "$gene: Replaced \'domain of unknown function\' by \'protein of unknown function\' in:\n$desc\n";
  $desc =~ s/domain of unknown function/protein of unknown function/g;
 }
 if($desc =~ /SET and RING finger associated domain. Domain of unknown function in SET domain containing proteins and in Deinococcus radiodurans DRA1533\./){
  print "$gene: Replaced \'SET and RING finger associated domain. Domain of unknown function in SET domain containing proteins and in Deinococcus radiodurans DRA1533\.\' by \'Protein of unknown function\' in:\n$desc\n";
  $desc =~ s/SET and RING finger associated domain. Domain of unknown function in SET domain containing proteins and in Deinococcus radiodurans DRA1533\./Protein of unknown function/;
 }
 if(/Domain of unknown function in PX-proteins/){
  print "$gene: Replaced \'Domain of unknown function in PX-proteins\' by \'Protein of unknown function\' in:\n$desc\n";
  $desc =~ s/Domain of unknown function in PX-proteins/Protein of unknown function/;
 }
 if(/Domain of unknown function in Sec63p, Brr2p and other proteins\./){
  print "$gene: Replaced \'Domain of unknown function in Sec63p, Brr2p and other proteins\' by \'Protein of unknown function\' in:\n$desc\n";
  $desc =~ s/Domain of unknown function in Sec63p, Brr2p and other proteins\./Protein of unknown function/;
 }
 if($desc =~ /Domain of unknown function/){
  print "$gene: Replaced \'Domain of unknown function\' by \'Protein of unknown function\' in:\n$desc\n";
  $desc =~ s/Domain of unknown function/Protein of unknown function/g;
 }
 if($desc =~ /Predicted/){
  print "$gene: Replaced \'Predicted\' by \'Putative\' in:\n$desc\n";
  $desc =~ s/Predicted/Putative/g;
 }
 if($desc =~ /Uncharacterized/){
  print "$gene: Replaced \'Uncharacterized\' by \'Putative\' in:\n$desc\n";
  $desc =~ s/Uncharacterized/Putative/g;
 }
 if(lc($desc) eq 'zinc finger'){
  print "Renamed ($gene): \'$desc\' to \'Zinc finger protein\'\n";
  $desc='Zinc finger protein';
 }
 
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
    $featIdentifier =~ s/\./\_/g;
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
   $addInfoField =~ s/\./\_/g;
   unless($feattype eq 'gene'){
    $featIdentParent=$addInfoField;
    $featIdentParent =~ s/Parent=//g;
    $featuresInfo{$featIdentifier}{'parent'}=$featIdentParent;
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
my $locusCount=0;
my %gene2locusTag;

# Open tbl file
open(TBLFILE,">",$tblOut);

print TBLFILE ">Features	SeqID	table_name\n";

foreach my $seq (@sequences){
 print TBLFILE ">Feature	$seq	Table1\n";
 foreach my $feat (@{$featuresSeqs{$seq}}){
  $feat=~s/\./\_/g;
  if($featuresInfo{$feat}{'feattype'} eq 'CDS'){
   foreach my $cdsNum (keys $featuresInfo{$feat}{'CDS'}){
    if($featuresInfo{$feat}{'CDS'}{$cdsNum}{'strand'} eq '-'){
     print TBLFILE "$featuresInfo{$feat}{'CDS'}{$cdsNum}{'end'}\t$featuresInfo{$feat}{'CDS'}{$cdsNum}{'init'}\t$featuresInfo{$feat}{'feattype'}\n";
    } else {
     print TBLFILE "$featuresInfo{$feat}{'CDS'}{$cdsNum}{'init'}\t$featuresInfo{$feat}{'CDS'}{$cdsNum}{'end'}\t$featuresInfo{$feat}{'feattype'}\n";
    }
    print TBLFILE "			codon_start	$featuresInfo{$feat}{'CDS'}{$cdsNum}{'codon_start'}\n";
    my $feat2=$feat;
    $feat2 =~ s/cds\_//g;
    my $feat_preRemoved=$featuresInfo{$feat2}{'parent'};
    $feat_preRemoved=~s/evm_TU_//g;
    print TBLFILE "			gene	$feat_preRemoved\n";
    print TBLFILE "			protein_id	$feat_preRemoved\n";
    if($gene2locusTag{$featuresInfo{$feat2}{'parent'}}){
     print TBLFILE "			locus_tag	$gene2locusTag{$featuresInfo{$feat2}{'parent'}}\n";
    } else {
     die "There is no LOCUS TAG for CDS: $feat\n";
    }
    if($interProScanAnnotation{$featuresInfo{$feat2}{'parent'}}){
     foreach my $note (@{$interProScanAnnotation{$featuresInfo{$feat2}{'parent'}}}){
      print TBLFILE "			product	$note\n";
     }
    } else {
     print TBLFILE "			note	hypothetical protein\n";
    }
   }
  } else {
   if ($featuresInfo{$feat}{'strand'} eq '-'){
    print TBLFILE "$featuresInfo{$feat}{'end'}\t$featuresInfo{$feat}{'init'}\t$featuresInfo{$feat}{'feattype'}\n";
   } else {
    print TBLFILE "$featuresInfo{$feat}{'init'}\t$featuresInfo{$feat}{'end'}\t$featuresInfo{$feat}{'feattype'}\n";
   }
   if($featuresInfo{$feat}{'feattype'} eq 'exon'){
    my $feat2=$featuresInfo{$feat}{'parent'};
    my $feat_preRemoved=$featuresInfo{$feat2}{'parent'};
    $feat_preRemoved=~s/evm_TU_//g;
    print TBLFILE "			gene	$feat_preRemoved\n";
    if($gene2locusTag{$featuresInfo{$feat2}{'parent'}}){
     print TBLFILE "			locus_tag	$gene2locusTag{$featuresInfo{$feat2}{'parent'}}\n";
    } else {
     die "There is no LOCUS TAG for exon: $feat\n";
    }
    if($interProScanAnnotation{$featuresInfo{$feat2}{'parent'}}){
     foreach my $note (@{$interProScanAnnotation{$featuresInfo{$feat2}{'parent'}}}){
      print TBLFILE "			product $note\n";
     }
    } else {
     print TBLFILE "			note	hypothetical protein\n";
    }
    if($parent2exon{$featuresInfo{$feat}{'parent'}}){
     $parent2exon{$featuresInfo{$feat}{'parent'}}++;
     print TBLFILE "			number	$parent2exon{$featuresInfo{$feat}{'parent'}}\n";
    } else {
     $parent2exon{$featuresInfo{$feat}{'parent'}}=1;
     print TBLFILE "			number	$parent2exon{$featuresInfo{$feat}{'parent'}}\n";
    }
   } elsif($featuresInfo{$feat}{'feattype'} eq 'gene'){
    my $feat_preRemoved=$feat;
    $feat_preRemoved=~s/evm_TU_//g;
    print TBLFILE "			gene	$feat_preRemoved\n";
    $locusCount++;
    $gene2locusTag{$feat}="KALMBRA_".$locusCount;
    print TBLFILE "			locus_tag	$gene2locusTag{$feat}\n";
   } elsif($featuresInfo{$feat}{'feattype'} eq 'mRNA') {
    my $feat_preRemoved=$featuresInfo{$feat}{'parent'};
    $feat_preRemoved=~s/evm_TU_//g;
    print TBLFILE "			gene	$feat_preRemoved\n";
    print TBLFILE "			locus_tag	$gene2locusTag{$featuresInfo{$feat}{'parent'}}\n";
   } else {
    die "Something wrong...\n";
   }
  }
 }
}

close(TBLFILE);
