#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

## Versions
# 0.1
# 0.2 changes:
## - included command-line arguments (Getopt long module)
## -

my $version='0.2';
my $help='';
my $license='';
my $evmGffFile = '';
my $interproScanFile = '';
my $scafLengthsFile = '';
my $gapNsFile = '';
my $tblOut = '';

GetOptions(
  'help|h|?'            => \$help,
  'license|l'           => \$license,
  'evm_gff|e=s'         => \$evmGffFile,
  'interpro_tsv|i=s'    => \$interproScanFile,
  'scaf_lengths|sl=s'   => \$scafLengthsFile,
  'gaps=s'              => \$gapNsFile,
  'tbl_out|o=s'         => \$tblOut,
);

if(-s $tblOut) {
  print "The output file (feature table) already exists.\nPlease, delete this file before running the script again.\n";
  &usage();
  exit(1);
}

if(!$evmGffFile) {
  print "A GFF from EVM is required.\n";
  &usage();
  exit(1);
}

if(!$interproScanFile) {
  print "A TSV (tab-separated) file from InterProScan5 is required.\n";
  &usage();
  exit(1);
}

if($help) {
  &usage();
  exit(0);
}

if($license) {
 &license();
 exit(0);
}

# Information about gaps in sequences
my %gapsInSeqs;
my %numGapsInSeqs;

open(GAPSFILE,$gapNsFile);

while(<GAPSFILE>){
 chomp;
 my($seq,$init,$end)=split(/\t/,$_);
 if($gapsInSeqs{$seq}){
  $numGapsInSeqs{$seq}++;
  $gapsInSeqs{$seq}{$numGapsInSeqs{$seq}}{'init'}=$init;
  $gapsInSeqs{$seq}{$numGapsInSeqs{$seq}}{'end'}=$end;
 } else {
  $numGapsInSeqs{$seq}=1;
  $gapsInSeqs{$seq}{$numGapsInSeqs{$seq}}{'init'}=$init;
  $gapsInSeqs{$seq}{$numGapsInSeqs{$seq}}{'end'}=$end;
 }
}

close(GAPSFILE);

# Information InterProScan

my %interProScanAnnotation;
my %interProScanAnnotDbXref;
my %interProScanECcode;

open(INTERPROSCAN,$interproScanFile);

while(<INTERPROSCAN>){
 chomp;
 #evm.model.KI545862.1.400        c1932f5ac3f6f4274218e5d6b2427aa4        728     Phobius TRANSMEMBRANE   Region of a membrane-bound protein predicted to be embedded in the membrane.    367     384     -       T       07-04-2016
 my @fields=split(/\t/,$_);
 # Gene identifiers have 'TU', not 'model'
 my($gene,undef,undef,$db,$dbID,$desc,undef,undef,undef,undef,undef)=split(/\t/,$_);
 $gene =~ s/\.model\./\.TU\./g;
 $gene =~ s/\./\_/g;

 # Replace any problematic description (reported in disc.report)
 # tbl2asn run: ~/Software/tbl2asn/linux64.tbl2asn -p . -j "[organism=Kalmanozyma brasiliensis] [strain=GHG001]" -M n -y "Re-annotation of K. brasiliensis GHG001 including RNAseq experimental data" -i GCA_000497045.1_PSEUBRA1_genomic.fna -Z disc.report -t template.sbt -V b

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

 # NCBI does not like features with more than 100 characters
 if(length($desc) >= 100){
  print "$gene: Feature with more than 100 characters. Renamed from \'$desc\' to \'hypothetical protein\'\n";
  $desc = 'hypothetical protein';
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
  print "$gene: Found '\@'. Renamed from \'$desc\' by 'hypothetical protein'\n";
  $desc = 'hypothetical protein';
 }
 #features starts with 'hypothetical protein' but not equals 'hypothetical protein'
 if(lc($desc) =~ /hypothetical protein/){
  print "$gene: Renamed from \'$desc\' to \'hypothetical protein\' in:\n$desc\n";
  $desc = 'hypothetical protein';
 }
 # Features with ';'
 if($desc eq 'ARF-like small GTPases; ARF, ADP-ribosylation factor'){
  print "$gene: Renamed from \'$desc\' to \'ARF, ADP-ribosylation factor\'\n";
  $desc = 'ARF, ADP-ribosylation factor';
 }

 # Use American spelling
 if($desc =~ /organisation/){
  print "$gene: Replaced \'organisation\' by \'organization\' in:\n$desc\n";
  $desc =~ s/organisation/organization/g;
 }
 if($desc =~ /tetramerisation/){
  print "$gene: Replaced \'tetramerisation\' by \'tetramerization\' in:\n$desc\n";
  $desc =~ s/tetramerisation/tetramerization/;
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
  print "$gene: Renamed from \'$desc\' to \'hypothetical protein\'\n";
  $desc='hypothetical protein';
 }
 # NCBI does not like feature starting with 'region'
 if($desc =~ /Region of a membrane-bound protein/){
  print "$gene: Found 'region of a membrane-bound' at the beginning of the feature. Renamed from \'$desc\' to \'Membrane-bound protein\'\n";
  $desc='Membrane-bound protein';
 }
 if($desc =~ / gene /){
  print "$gene: Found \'gene\'. Renamed from \'$desc\' to \'hypothetical protein\'\n";
  $desc='hypothetical protein';
 }
 if($desc =~ /Kinetochore CENP-C fungal homologue, Mif2, N-terminal/){
  print "$gene: Renamed from \'Kinetochore CENP-C fungal homologue, Mif2, N-terminal\' to \'Kinetochore CENP-C, Mif2\'\n";
  $desc='Kinetochore CENP-C, Mif2';
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

 # Add description to gene
 # If available, include E.C. number 
 unless ($desc eq /^$/){
  if(($desc =~ /(EC \d+\.\d+\.\d+\.\d+)/) or ($desc =~ /(EC \d+\.\d+\.\d+\.\-)/)){
   $desc = $1;
   $desc =~ s/EC //g; 
   if($interProScanECcode{$gene}){
    unless($desc ~~ @{$interProScanECcode{$gene}}){
     push(@{$interProScanECcode{$gene}},$desc);
    }
   } else {
    @{$interProScanECcode{$gene}}=($desc);
   }
  } else {
   if($interProScanAnnotation{$gene}){
    unless($desc ~~ @{$interProScanAnnotation{$gene}}){
     push(@{$interProScanAnnotation{$gene}},$desc);
    }
   }else{
    @{$interProScanAnnotation{$gene}}=($desc);
   }
  }
 }

 # Databases must be listed here: https://www.ncbi.nlm.nih.gov/genbank/collab/db_xref/
 # Currently (August 2017), the following DBs are allowed: PFAM, TIGRFAM
 # In the .tbl file the db_ref must be capitalized
 unless (($dbID eq /^$/) or ($db eq /^$/)){
  if((uc($db) eq 'PFAM') or (uc($db) eq 'TIGRFAM')){
   $db = uc($db);
   my $tmpDBXRef="$db"."\t"."$dbID";
   if($interProScanAnnotDbXref{$gene}){
    unless($tmpDBXRef ~~ @{$interProScanAnnotDbXref{$gene}}){
     push(@{$interProScanAnnotDbXref{$gene}},$tmpDBXRef);
    }
   }else{
    @{$interProScanAnnotDbXref{$gene}}=($tmpDBXRef);
   }
  }
 }
}

close(INTERPROSCAN);

# Information GFF3

#KI545851.1      EVM     gene    695     2101    .       -       .       ID=evm.TU.KI545851.1.1;Name=EVM_prediction_KI545851.1.1
#KI545851.1      EVM     mRNA    695     2101    .       -       .       ID=evm.model.KI545851.1.1;Parent=evm.TU.KI545851.1.1;Name=EVM_prediction_KI545851.1.1
#KI545851.1      EVM     exon    695     2101    .       -       .       ID=evm.model.KI545851.1.1.exon1;Parent=evm.model.KI545851.1.1
#KI545851.1      EVM     CDS     695     2101    .       -       0       ID=cds.evm.model.KI545851.1.1;Parent=evm.model.KI545851.1.1

# Information about sequences (identifier in FASTA and length)
my @sequences=();
my %seqLengthInit;
my %seqLengthEnd;

# 
my @features=();
my %featuresInfo;
my %CDScount;
my %featuresSeqs;

open(SEQFILE,$scafLengthsFile);

while(<SEQFILE>){
 chomp;
 my ($seq,$seqInit,$seqEnd)=split(/\t/,$_);
 die "There are two sequences ($seq) with the same identifier in length file!\n" if ($seq ~~ @sequences);
 push(@sequences,$seq);
 $seqLengthInit{$seq}=$seqInit;
 $seqLengthEnd{$seq}=$seqEnd;
}

close(SEQFILE);

open(EVMGFFFILE,$evmGffFile);

while(<EVMGFFFILE>){
 chomp;
 next if(/^#/);
 my ($seq,$source,$feattype,$init_pos,$end_pos,undef,$strand,$codon_start,$additionalInfo)=split(/\t/,$_);
 die unless(scalar(split(/\t/,$_)) == 9);
 die unless($seq ~~ @sequences);
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
    die "Something wrong: there is more than one \"ID\" in the last field of the GFF3 file.\n";
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
 print TBLFILE "$seqLengthInit{$seq}	 $seqLengthEnd{$seq}	source\n";
 print TBLFILE "			mol_type	genomic DNA\n";
 print TBLFILE "			organism	Kalmanozyma brasiliensis\n";

 # Print gaps (assembly_gap)
 if($numGapsInSeqs{$seq}){
  my $numGaps = $numGapsInSeqs{$seq};
  foreach my $gapnum (1 .. $numGaps) {
   print TBLFILE "$gapsInSeqs{$seq}{$gapnum}{'init'}\t$gapsInSeqs{$seq}{$gapnum}{'end'}\tassembly_gap\n";
   print TBLFILE "			gap_type	within scaffold\n";
   print TBLFILE "			estimated_length	unknown\n";
   print TBLFILE "			linkage_evidence	paired-ends\n";
  }
 }

 foreach my $feat (@{$featuresSeqs{$seq}}){
  $feat=~s/\./\_/g;

  if($featuresInfo{$feat}{'feattype'} eq 'CDS'){

   my @initPositionsCDS=();
   my @endPositionsCDS=();
   my $startPosCDS=0;

   foreach my $cdsNum (keys $featuresInfo{$feat}{'CDS'}){
    push(@initPositionsCDS,$featuresInfo{$feat}{'CDS'}{$cdsNum}{'init'});
   }
   
   foreach my $cdsNum (keys $featuresInfo{$feat}{'CDS'}){
    push(@endPositionsCDS,$featuresInfo{$feat}{'CDS'}{$cdsNum}{'end'});
   }

   my @sortEndPositionsCDS=();
   my @sortInitPositionsCDS=();

   if($featuresInfo{$feat}{'CDS'}{1}{'strand'} eq '-') {
    @sortEndPositionsCDS = sort {$b <=> $a} @initPositionsCDS;
    @sortInitPositionsCDS = sort {$b <=> $a} @endPositionsCDS;
   } else {
    @sortInitPositionsCDS = sort {$a <=> $b} @initPositionsCDS;
    @sortEndPositionsCDS = sort {$a <=> $b} @endPositionsCDS;
   }

   foreach my $cdsLinePrintPositions (0 .. scalar(@sortInitPositionsCDS)-1){
    if($featuresInfo{$feat}{'CDS'}{1}{'strand'} eq '-'){
     if($startPosCDS==0){
      $startPosCDS=1;
      print TBLFILE "$sortInitPositionsCDS[$cdsLinePrintPositions]\t$sortEndPositionsCDS[$cdsLinePrintPositions]\tCDS\n";
     }else{
      print TBLFILE "$sortInitPositionsCDS[$cdsLinePrintPositions]\t$sortEndPositionsCDS[$cdsLinePrintPositions]\n";
     }
    } else {
     if($startPosCDS==0){
      $startPosCDS=1;
      print TBLFILE "$sortInitPositionsCDS[$cdsLinePrintPositions]\t$sortEndPositionsCDS[$cdsLinePrintPositions]\tCDS\n";
     }else{
      print TBLFILE "$sortInitPositionsCDS[$cdsLinePrintPositions]\t$sortEndPositionsCDS[$cdsLinePrintPositions]\n";
     }
    }
   }
       
   my $feat2=$feat;
   $feat2 =~ s/cds\_//g;
   print TBLFILE "			codon_start	$featuresInfo{$feat}{'CDS'}{1}{'codon_start'}\n";
   print TBLFILE "			protein_id	gnl|BCE_CTBE|$gene2locusTag{$featuresInfo{$feat2}{'parent'}}\n";
   print TBLFILE "			transcript_id	gnl|BCE_CTBE|mrna.$gene2locusTag{$featuresInfo{$feat2}{'parent'}}\n";
   if($gene2locusTag{$featuresInfo{$feat2}{'parent'}}){
    print TBLFILE "			locus_tag	$gene2locusTag{$featuresInfo{$feat2}{'parent'}}\n";
    #print TBLFILE "			gene	$gene2locusTag{$featuresInfo{$feat2}{'parent'}}\n";
   } else {
    die "There is no LOCUS TAG for CDS: $feat\n";
   }
   if($interProScanAnnotation{$featuresInfo{$feat2}{'parent'}}){
    foreach my $note (@{$interProScanAnnotation{$featuresInfo{$feat2}{'parent'}}}){
     # If there is a EC associated with this protein, it cannot be annotated as 'hypothetical'
     unless (($interProScanECcode{$featuresInfo{$feat2}{'parent'}}) and ($note eq "hypothetical protein")) {
      print TBLFILE "			product	$note\n";
     }
    }
   } else {
    print TBLFILE "			note	hypothetical protein\n";
   }
   # Print (if any) EC number associated with this CDS
   if($interProScanECcode{$featuresInfo{$feat2}{'parent'}}){
    foreach my $ec (@{$interProScanECcode{$featuresInfo{$feat2}{'parent'}}}){
     print TBLFILE "			EC_number	$ec\n";
    }
   }
   # Parsing db_xref information
   if($interProScanAnnotDbXref{$featuresInfo{$feat2}{'parent'}}){
    foreach my $db_xref (@{$interProScanAnnotDbXref{$featuresInfo{$feat2}{'parent'}}}){
     my ($db,$dbID)=split(/\t/,$db_xref);
     print TBLFILE "			db_xref	$db:$dbID\n";
    }
   }
  } # Closing CDS

  elsif ($featuresInfo{$feat}{'feattype'} eq 'gene') {
   $locusCount++;
   my $final_locusTagCount='';
   #old locus_tag has five numbers: '00098'
   if (length($locusCount) == 1){
    $final_locusTagCount='00000'.$locusCount;
   } elsif (length($locusCount) == 2) {
    $final_locusTagCount='0000'.$locusCount;
   } elsif (length($locusCount) == 3) {
    $final_locusTagCount='000'.$locusCount;
   } elsif (length($locusCount) == 4) {
    $final_locusTagCount='00'.$locusCount;
   } elsif (length($locusCount) == 5) {
    $final_locusTagCount='0'.$locusCount;
   } else {
    die "LOCUS TAG number (\$locusCount) seems to be above the allowed.\n";
   }
   $gene2locusTag{$feat}="PSEUBRA_".$final_locusTagCount;
   if ($featuresInfo{$feat}{'strand'} eq '-'){
    print TBLFILE "$featuresInfo{$feat}{'end'}\t$featuresInfo{$feat}{'init'}\n";
   } else {
    print TBLFILE "$featuresInfo{$feat}{'init'}\t$featuresInfo{$feat}{'end'}\n";
   }
   #print TBLFILE "			gene	$gene2locusTag{$feat}\n";   
   print TBLFILE "			locus_tag	$gene2locusTag{$feat}\n";
  } # Closing gene

  elsif ($featuresInfo{$feat}{'feattype'} eq 'exon') {
   if ($featuresInfo{$feat}{'strand'} eq '-'){
    print TBLFILE "$featuresInfo{$feat}{'end'}\t$featuresInfo{$feat}{'init'}\t$featuresInfo{$feat}{'feattype'}\n";
   } else {
    print TBLFILE "$featuresInfo{$feat}{'init'}\t$featuresInfo{$feat}{'end'}\t$featuresInfo{$feat}{'feattype'}\n";
   }
   my $feat2=$featuresInfo{$feat}{'parent'};
   if($gene2locusTag{$featuresInfo{$feat2}{'parent'}}){
    print TBLFILE "			locus_tag	$gene2locusTag{$featuresInfo{$feat2}{'parent'}}\n";
    #print TBLFILE "			gene	$gene2locusTag{$featuresInfo{$feat2}{'parent'}}\n";
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
  } # Closing exon


  elsif ($featuresInfo{$feat}{'feattype'} eq 'mRNA') {
   if ($featuresInfo{$feat}{'strand'} eq '-'){
    print TBLFILE "$featuresInfo{$feat}{'end'}\t$featuresInfo{$feat}{'init'}\t$featuresInfo{$feat}{'feattype'}\n";
   } else {
    print TBLFILE "$featuresInfo{$feat}{'init'}\t$featuresInfo{$feat}{'end'}\t$featuresInfo{$feat}{'feattype'}\n";
   }
   #print TBLFILE "			gene	$gene2locusTag{$featuresInfo{$feat}{'parent'}}\n";
   print TBLFILE "			protein_id	gnl|BCE_CTBE|$gene2locusTag{$featuresInfo{$feat}{'parent'}}\n";
   print TBLFILE "			transcript_id	gnl|BCE_CTBE|mrna.$gene2locusTag{$featuresInfo{$feat}{'parent'}}\n";
   print TBLFILE "			locus_tag	$gene2locusTag{$featuresInfo{$feat}{'parent'}}\n";
  }

  else { # If it is not gene, mRNA, exon, or CDS: something is wrong!
   die "Something wrong...\n";
  }

 } # Close the foreach feature
}

close(TBLFILE);

sub usage {
    print STDERR <<EOF;

NAME
    $0 version $version
    $0 takes an EVM GFF file and InterProScan5 results, and generates a tbl for NCBI annotation submission (tbl2asn input)

BASIC USAGE
    $0 --evm_gff evmannotation.gff --interpro_tsv interproannot.tsv --tbl_out organismannot.tbl

OPTIONS
    --evm_gff        -e      EVM input file in the GFF format                                 REQUIRED
    --interpro_tsv   -i      InterProScan5 results output file in the TSV (tab-separed)       REQUIRED
    --scaf_lengths   -sl     Scaffold lengths used as input                                   REQUIRED ?
    --gaps                   Position of gaps in scaffolds                                    REQUIRED ?
    --tbl_out        -o      Output file in the file, as required by tbl2asn (feature table)
    --help,          -h      This help.
    --license        -l      License.

EOF
}

# Subroutine that prints 
sub license{
    print STDERR <<EOF;
Copyright (C) 2017 Renato Augusto Correa dos Santos
e-mail: renatoacsantos\@gmail.com

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
