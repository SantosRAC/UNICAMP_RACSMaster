#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

my %familyTPs; # Identified as positive, from the positive dataset
my %familyFPs; # Identified as positive, from the non-positive dataset
my %familyFNs; # Identified as negative, from the positive dataset
my %familyMembers;

my $license='';
my $table='';
my $mapfile='';
my $outfile='';
my $help='';
my $debug='';
my $version='0.1';

#Get input from user
GetOptions(
    'license|l'        => \$license,
    'intable|i=s'      => \$table,
    'mapfile|m=s'      => \$mapfile,
    'outfile|o=s'      => \$outfile,
    'help|h|?'         => \$help,
    'debug|d=i'        => \$debug
);

#Check input from user
if ($help){
    &usage();
    exit(1);
}
if ($license){
    &license();
    exit(1);
}
if (!-e $table){
    print STDERR "FATAL: You must provide an input table generated after parsing hmmsearch and mapping files.\n";
    &usage();
    exit(1);
}
if (!$outfile){
    print STDERR "FATAL: You must provide a name for the output file.\n";
    &usage();
    exit(1);
}
if (-e $outfile){
    print STDERR "FATAL: The output file already exists. Delete the older one, or use another name.\n";
    &usage();
    exit(1);
}
if (!-e $mapfile){
    print STDERR "FATAL: You must provide valid mapping file with TPs (proteins and corresponding families in TCDB).\n";
    &usage();
    exit(1);
}

########################################
# Mapping table with proteins in TCDB and corresponding family or subfamily in which it was classified
########################################
open(MAPFILE,$mapfile);

while(<MAPFILE>){
 chomp;
 my($protein,$family) = split(/\t/,$_);
 unless($protein ~~ @{$familyMembers{$family}}){
  push(@{$familyMembers{$family}},$protein);
 }
}

close(MAPFILE);

########################################
# Resulting table, processed after running HMMER
########################################
open(TABLE,$table);

while(<TABLE>){
 chomp;
 next if (/^Status/);
 # HEADER OF THE TABLE: Status  TCDB_Protein    TCDB_Family     HMMs    Sequence_Evalue
 my ($Status,$TCDB_Protein,$TCDB_Family,$HMMs,$Sequence_Evalue) = split(/\t/,$_);
 if($Status eq 'OK_merged'){
  my @familiesMerge=split(/,/,$HMMs);
  foreach my $f (@familiesMerge){
   if($TCDB_Protein ~~ @{$familyMembers{$f}}){
    next if($TCDB_Protein ~~ @{$familyTPs{$f}});
    push(@{$familyTPs{$f}},$TCDB_Protein);
   } else {
    die "Something wrong with the script used previously. Protein $TCDB_Protein with the status \'OK_merged\' is not in TPs (family: $f)\n";
   }
  }
 }
 elsif ($Status eq 'OK'){
  if($TCDB_Protein ~~ @{$familyMembers{$TCDB_Family}}){
   next if($TCDB_Protein ~~ @{$familyTPs{$TCDB_Family}});
   push(@{$familyTPs{$TCDB_Family}},$TCDB_Protein);
  } else {
   die "Something wrong with the script used previously. Protein $TCDB_Protein with the status \'OK\' is not in TPs (family: $TCDB_Family)\n";
  }
 }
 elsif ($Status eq 'NO_HMM'){
  if($TCDB_Protein ~~ @{$familyMembers{$TCDB_Family}}){
   next if($TCDB_Protein ~~ @{$familyFNs{$TCDB_Family}});
   push(@{$familyFNs{$TCDB_Family}},$TCDB_Protein);
  } else {
   die "Something wrong with the script used previously. Protein $TCDB_Protein with the status \'OK\' is not in TPs\n";
  }
 }
 elsif ($Status eq 'Other') {
  if($TCDB_Protein ~~ @{$familyMembers{$TCDB_Family}}){
   if($HMMs =~ /clstr/){
    my @HMMsCorName=split(/_/,$HMMs);
    my $HMMs=$HMMsCorName[0];
   }
   next if($TCDB_Protein ~~ @{$familyFPs{$HMMs}});
   push(@{$familyFPs{$HMMs}},$TCDB_Protein);
  } else {
    die "Something wrong with the script used previously. Protein $TCDB_Protein with the status \'Other\' is not in TPs (family: $TCDB_Family)\n";
  }
 }
 elsif ($Status eq 'Cross'){
  next;
 } else {
  die "Something wrong in line $_\n";
 }
}

close(TABLE);

open(OUTFILE,">",$outfile);

print OUTFILE "Family\tValue\tType_Of_Value\n";

#Computation of sensitivity and PPV
foreach my $family (keys %familyMembers){
 my $Sensitivity='';
 my $PPV='';
 if($familyTPs{$family}){
  my $TPcount = scalar(@{$familyTPs{$family}});
  if($familyFNs{$family}){
   my $FNcount = scalar(@{$familyFNs{$family}});
   $Sensitivity = $TPcount/($TPcount+$FNcount);
   print OUTFILE "$family\t$Sensitivity\tsensitivity\n";
  } else {
   $Sensitivity = $TPcount/($TPcount+0);
   print OUTFILE "$family\t$Sensitivity\tsensitivity\n";
  }
  if($familyFPs{$family}){
   my $FPcount = scalar(@{$familyFPs{$family}});
   $PPV = $TPcount/($TPcount+$FPcount);
   print OUTFILE "$family\t$PPV\tPPV\n";
  } else {
   $PPV = $TPcount/($TPcount+0);
   print OUTFILE "$family\t$PPV\tPPV\n";
  }
 } else {
  print OUTFILE "$family\tNO_TRUE_POSITIVES\tNO_TRUE_POSITIVES\n";
 }
}

close(OUTFILE);

sub usage{
    print STDERR "$0 version $version, Copyright (C) 2016 Renato Augusto Correa dos Santos\n";
    print STDERR "$0 comes with ABSOLUTELY NO WARRANTY; for details type `$0 -l'.\n";
    print STDERR "This is free software, and you are welcome to redistribute it under certain conditions;\n";
    print STDERR "type `$0 -l' for details.\n";
    print STDERR <<EOF;
NAME
    $0

USAGE
    $0 --outfile test.txt --mapfile mappingFamilyProteins.txt --intable hmmsearchSummary.txt 

OPTIONS
    --intable       -i    REQUIRED
    --mapfile       -m    REQUIRED
    --outfile       -o    REQUIRED
    --help          -h    This help
    --license       -l    License

EOF
}

sub license{
    print STDERR <<EOF;

 Copyright (C) 2016 Renato Augusto Correa dos Santos
 http://bce.bioetanol.org.br/
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
