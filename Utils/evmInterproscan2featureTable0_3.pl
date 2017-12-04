#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

## Versions
# 0.1
# 0.2 changes:
## - included command-line arguments (Getopt long module)
# 0.3 changes
## - get information from the output of seqkit ($ seqkit locate --ignore-case --only-positive-strand --pattern "N+" GCA_000497045.1_PSEUBRA1_genomic.fna| sed 's/ +/       /g'| cut -f1,5,6)
## - assume the file with positions of Ns is required
## - get information about products from the IPR description (column besides IPR identifier)
## - fix product names
## - how tbl2asn (version available in NCBI on Nov-27-17) is run: 
### $ linux64.tbl2asn -j "[organism=Kalmanozyma brasiliensis] [strain=GHG001]" -M n -y "Re-annotation of K. brasiliensis GHG001 including RNAseq experimental data" -i GCA_000497045.1_PSEUBRA1_genomic.fna -Z disc.report -t template.sbt -V b -a r10u -l paired-ends

my $version='0.2';
my $help='';
my $license='';
my $evmGffFile = '';
my $interproScanFile = '';
my $scafLengthsFile = '';
my $gapNsFile = '';
my $tblOut = '';
my $logFile = '';
my $locusTag = '';

# Information about gaps in sequences
my %gapsInSeqs;
my %numGapsInSeqs;

# Information about sequences (identifier in FASTA and lengths)
my @scaf_sequences=();
my %seqLengthInit;
my %seqLengthEnd;

GetOptions(
  'help|h|?'            => \$help,
  'license|l'           => \$license,
  'evm_gff|e=s'         => \$evmGffFile,
  'interpro_tsv|i=s'    => \$interproScanFile,
  'scaf_lengths|sl=s'   => \$scafLengthsFile,
  'log_file=s'          => \$logFile,
  'locus_tag=s'         => \$locusTag,
  'gaps=s'              => \$gapNsFile,
  'tbl_out|o=s'         => \$tblOut,
);

if(!$tblOut) {
  print "A name for the output file (.tbl, a feature table) is required.\n";
  &usage();
  exit(1);
}

if(-s $tblOut) {
  print "The output file ($tblOut) already exists.\nPlease, delete this file before running the script again.\n";
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

if(!$scafLengthsFile){
 print "A file with lengths of scaffolds in the EVM GFF is required.\n";
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

if(!$logFile){
 $logFile='log_evm2tbl.txt';
}

if(-s $logFile) {
  print "The log file ($logFile) already exists.\nPlease, delete this file before running the script again.\n";
  &usage();
  exit(1);
}

if(!$locusTag){
 print "User must provide a locus tag.\n";
 &usage();
 exit(1);
}

if(!$gapNsFile){
 print "User must provide gap positions (runs of Ns).\n";
 &usage();
 exit(1);
}

# If user passes file with gap positions available in an external file
# This part opens this file, storing 
open(GAPSFILE,$gapNsFile);

while(<GAPSFILE>){
 chomp;
 next if (/^seqID/);
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
open(LOGFILE,'>',$logFile);

while(<INTERPROSCAN>){
 chomp;
 # evm.model.KI545862.1.400	c1932f5ac3f6f4274218e5d6b2427aa4	728	SUPERFAMILY	SSF74788		418727	4.92E-86	T	07-04-2016	IPR016159	Cullin repeat-like-containing domain		

 my @fields=split(/\t/,$_);
 my($gene,undef,undef,$db,$dbID,undef,undef,undef,undef,undef,undef,$ipr,$desc)=split(/\t/,$_);
 if ($ipr){
  if ($ipr =~ /^IPR/) {
   $gene =~ s/\.model\./\.TU\./g;
   $gene =~ s/\./\_/g;

   #Fixing description
   ## Remove organism from product name
   ### feature contains 'fungi', metazoa/fungi, fungi/archaea, 'eukaryotes'
   ### feature contains 'staphylococcal'
   ## Typo
   ### feature contains heam- (should be replaced by hem-)
   ### 'homologue' implies evolutionary relationship/ change to -like
   ### Plurals 'purines', 'domains', 'synthetases'

   if(($desc =~ /Adenine nucleotide alpha hydrolase-like domains/) and ($ipr eq 'IPR023382')){
    print LOGFILE "$gene: Substituted \'domains\' by \'domain\' in:\n$desc ($ipr)\n";
    $desc =~ s/domains/domain/g;
   }
   if(($desc =~ /purines/) and ($ipr eq "IPR001248")){
    print LOGFILE "$gene: Substituted \'purines\' by \'purine\' in:\n$desc ($ipr)\n";
    $desc =~ s/purines/purine/g;
   }
   if(($desc =~ /Aspartyl-tRNA synthetases/) and ($ipr eq 'IPR004523')){
    print LOGFILE "$gene: Substituted \'synthetases\' by \'synthetase\' in:\n$desc ($ipr)\n";
    $desc =~ s/synthetases/synthetase/g;
   }
   if(($desc =~ /SLC41 divalent cation transporters, integral membrane domain/) and ($ipr eq 'IPR006667')){
    print LOGFILE "$gene: Substituted \'transporters\' by \'transporter\' in:\n$desc ($ipr)\n";
    $desc =~ s/transporters/transporter/g;
   }
   if($desc =~ /, eukaryotes$/){
    print LOGFILE "$gene: Removed \', eukaryotes\' at the end of:\n$desc ($ipr)\n";
    $desc =~ s/, eukaryotes//g;
   }
   if($desc =~ /, fungi$/){
    print LOGFILE "$gene: Removed \', fungi\' at the end of:\n$desc ($ipr)\n";
    $desc =~ s/, fungi//g;
   }
   if(($desc =~ /Staphylococcal nuclease/) and ($ipr eq 'IPR016071')){
    print LOGFILE "$gene: Substituted \'Staphylococcal nuclease\' by \'Nuclease\' in:\n$desc ($ipr)\n";
    $desc =~ s/Staphylococcal n/N/g;
   }
   if($desc =~ /, metazoa\/fungi/){
    print LOGFILE "$gene: Removed \', metazoa\/fungi\' at the end of:\n$desc\n";
    $desc =~ s/, metazoa\/fungi//g;
   }
   if($desc =~ /, fungi\/archaea/){
    print LOGFILE "$gene: Removed \', fungi\/archaea\' at the end of:\n$desc\n";
    $desc =~ s/, fungi\/archaea//g;
   }
   if($desc =~ / fungi$/){
    print LOGFILE "$gene: Removed \', , fungi\/archaea\' at the end of:\n$desc\n";
    $desc =~ s/ fungi//g;
   }
   if($desc =~ /heam-/){
    print LOGFILE "$gene: Substituted \'heam-\' by \'hem-\' in\n$desc\n";
    $desc =~ s/heam-/hem-/g;
   }
   if(($desc =~ /^Partial /) and ($ipr eq 'IPR006693')){
    print LOGFILE "$gene ($ipr): Removed \'Partial \' at the beginning of\n$desc ($ipr)\n";
    $desc =~ s/Partial //g;
   }
   if(($desc =~ /Activator of Hsp90 ATPase homologue 1-like/) and ($ipr eq 'IPR013538')){
    print LOGFILE "$gene ($ipr): Substituted \'ATPase homologue 1-like\' by \'ATPase 1-like\' in:\n$desc ($ipr)\n";
    $desc =~ s/Activator of Hsp90 ATPase homologue 1-like/Activator of Hsp90 ATPase 1-like/g;
   }
   if(($desc =~ /DnaJ homologue, subfamily C, member 28, conserved domain/) and ($ipr eq 'IPR018961')){
    print LOGFILE "$gene ($ipr): Substituted \'DnaJ homologue\' by \'DnaJ-like\' in:\n$desc ($ipr)\n";
    $desc =~ s/DnaJ homologue, subfamily C, member 28, conserved domain/DnaJ-like, subfamily C, member 28, conserved domain/g;
   }
   if(($desc =~ /DNA mismatch repair protein MutS-homologue MSH6/) and ($ipr eq 'IPR015536')){
    print LOGFILE "$gene ($ipr): Substituted \'MutS-homologue\' by \'MutS-like\' in:\n$desc ($ipr)\n";
    $desc =~ s/DNA mismatch repair protein MutS-homologue MSH6/DNA mismatch repair protein MutS-like MSH6/g;
   }
   if(($desc =~ /Protein notum homologue/) and ($ipr eq 'IPR004963')){
    print LOGFILE "$gene ($ipr): Substituted \'notum homologue\' by \'notum-like\' in:\n$desc ($ipr)\n";
    $desc =~ s/Protein notum homologue/Protein notum-like/g;
   }

   if($interProScanAnnotation{$gene}){
    unless($desc ~~ @{$interProScanAnnotation{$gene}}){
     push(@{$interProScanAnnotation{$gene}},$desc);
    }
   } else{
    @{$interProScanAnnotation{$gene}}=($desc);
   }
  }
 }

}

close(INTERPROSCAN);
close(LOGFILE);

my @features=();
my %featuresInfo;
my %CDScount;
my %featuresSeqs;

open(SEQFILE,$scafLengthsFile);

while(<SEQFILE>){
 chomp;
 my ($seq,$seqInit,$seqEnd)=split(/\t/,$_);
 die "There are two sequences ($seq) with the same identifier in length file!\n" if ($seq ~~ @scaf_sequences);
 push(@scaf_sequences,$seq);
 $seqLengthInit{$seq}=$seqInit;
 $seqLengthEnd{$seq}=$seqEnd;
}

close(SEQFILE);

open(EVMGFFFILE,$evmGffFile);

while(<EVMGFFFILE>){
 chomp;
 next if(/^#/);
 my ($seq,$source,$feattype,$init_pos,$end_pos,undef,$strand,$codon_start,$additionalInfo)=split(/\t/,$_);
 unless(scalar(split(/\t/,$_)) == 9){
  die scalar(split(/\t/,$_))." is the number of fields in line\n$_\n";
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

foreach my $seq (@scaf_sequences){
 print TBLFILE ">Feature	$seq	Table1\n";
 print TBLFILE "$seqLengthInit{$seq}	 $seqLengthEnd{$seq}	source\n";
 print TBLFILE "			mol_type	genomic DNA\n";
 print TBLFILE "			organism	Kalmanozyma brasiliensis\n";

=cut
 # Print gaps (assembly_gap), if external file is available
 if($gapNsFile){
  if($numGapsInSeqs{$seq}){
   my $numGaps = $numGapsInSeqs{$seq};
   foreach my $gapnum (1 .. $numGaps) {
    print TBLFILE "$gapsInSeqs{$seq}{$gapnum}{'init'}\t$gapsInSeqs{$seq}{$gapnum}{'end'}\tassembly_gap\n";
    print TBLFILE "			gap_type	within scaffold\n";
    print TBLFILE "			estimated_length	unknown\n";
    print TBLFILE "			linkage_evidence	paired-ends\n";
   }
  }
 }
=cut

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
   $gene2locusTag{$feat}=$locusTag."_".$final_locusTagCount;
   if ($featuresInfo{$feat}{'strand'} eq '-'){
    print TBLFILE "$featuresInfo{$feat}{'end'}\t$featuresInfo{$feat}{'init'}	gene\n";
   } else {
    print TBLFILE "$featuresInfo{$feat}{'init'}\t$featuresInfo{$feat}{'end'}	gene\n";
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
     print TBLFILE "			product	$note\n";
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
    $0 --evm_gff evmannotation.gff --interpro_tsv interproannot.tsv --tbl_out organismannot.tbl --scaf_lengths scaffolds_lengths.txt --locus_tag PSEUBRA --gaps runs_of_Ns.txt

OPTIONS
    --evm_gff        -e      EVM input file in the GFF format                                 REQUIRED
    --interpro_tsv   -i      InterProScan5 results output file in the TSV (tab-separed)       REQUIRED
    --scaf_lengths   -sl     Scaffold lengths used as input                                   REQUIRED
    --gaps                   Position of gaps in scaffolds                                    
      The script assumes a modified file after running seqkit, starting with
      seqID   start   end                                                                     REQUIRED
    --log_file               
      If the name of a log file is not set, it will be automatically named "log_evm2tbl.txt"  OPTIONAL
    --tbl_out        -o      Output file in the file, as required by tbl2asn (feature table)  REQUIRED
    --locus_tag              Locus tag                                                        REQUIRED
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
