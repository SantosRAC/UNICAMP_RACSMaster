#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use XML::Twig;

my $version='0.1';
my $infile='';
my $outfile='';
my $outtype='';
my $license='';
my $help='';
my $debug=0;

my @reactions=();
my @species=();

my %reactant2reactions;
my %reaction2products;

GetOptions(
    'license|l'    => \$license,
    'help|h|?'     => \$help,
    'debug|d:i'    => \$debug,
    'out_type:s'   => \$outtype,
    'outfile|o:s'    => \$outfile,
    'infile|i=s'   => \$infile
);

if($help){
	&usage();
    exit 1;
}
if($license){
	&license();
    exit 1;
}
if(!-s $infile){
    print STDERR "FATAL: you must provide a XML file as input\n";
    &usage();
    exit 0;	
}
if(!$outtype){
   print STDERR "FATAL: You must provide a type for the output file (edgelist or adjlist)\n";
   &usage();
   exit 0;
}
if(!$outfile){
   print STDERR "FATAL: You must provide a name for output file\n";
   &usage();
   exit 0;
}
if(-e $outfile){
   print STDERR "FATAL: Output file already exists\n";
   &usage();
   exit 0;
}

parseXML();

if ($outtype eq 'adjlist') {
 foreach my $key (keys %reaction2products) {
  open(OUTFILE,">>",$outfile);
  print OUTFILE "$key\t".join("\t",@{$reaction2products{$key}})."\n";
  close(OUTFILE);
 }
 
 foreach my $key (keys %reactant2reactions) {
  open(OUTFILE,">>",$outfile);
  print OUTFILE "$key\t".join("\t",@{$reactant2reactions{$key}})."\n";
  close(OUTFILE);
 }
}
elsif ($outtype eq 'edgelist') {
 foreach my $key (keys %reaction2products) {
  foreach my $value (@{$reaction2products{$key}}) {
   open(OUTFILE,">>",$outfile);
   print OUTFILE "$key\t$value\n";
   close(OUTFILE);
  }
 }

 foreach my $key (keys %reactant2reactions) {
  foreach my $value (@{$reactant2reactions{$key}}) {
   open(OUTFILE,">>",$outfile);
   print OUTFILE "$key\t$value\n";
   close(OUTFILE);
  }
 }
}
else {
 die "Unrecognized output type: \"$outtype\"\n";
}

sub parseXML{

 my $twig=XML::Twig->new();
 $twig->parsefile($infile) or die "cannot parse $infile";
 my $root = $twig->root;
 my @models = $root->children('model');

 foreach my $mod (@models){
  my @reactionLists = $mod->children('listOfReactions');
  foreach my $real (@reactionLists) {
   my @AllReactions = $real->children('reaction');
   foreach my $rea (@AllReactions) {
    my $id = '';
    $id = $rea->att('id');
    if ($id ~~ @reactions) {
     die "Duplicated reaction: $id\nCheck your input file\n";
    } else {
     push(@reactions,$id);
    }
    my @reactantList = $rea->children('listOfReactants');
    foreach my $reactl (@reactantList) {
     my @AllReactants = $reactl->children('speciesReference');
     my $numberReactants=0;
     foreach my $reaspe (@AllReactants) {
      $numberReactants++;
      if ($id ~~ @{$reactant2reactions{$reaspe->att('species')}}) {
       die "Reaction $id already exists in reactant .".$reaspe->att('species')."\n";
      } else {
       push(@{$reactant2reactions{$reaspe->att('species')}},$id);
      }
     }
    }
    my @productList = $rea->children('listOfProducts');
    foreach my $productl (@productList) {
     my @AllProducts = $productl->children('speciesReference');
     my $numberProducts=0;
     foreach my $productspe (@AllProducts) {
      $numberProducts++;
      if ($productspe->att('species') ~~ @{$reaction2products{$id}}) {
       die "Metabolite ".$productspe->att('species')." already exists in reaction $id\n";
      } else {
       push(@{$reaction2products{$id}},$productspe->att('species'));
      }
     }
    }
   }
  }
 }

}

##############################################
#Usage
##############################################
sub usage{
    print STDERR "$0 version $version, Copyright (C) 2016 Renato Augusto Correa dos Santos\n";
    print STDERR "$0 comes with ABSOLUTELY NO WARRANTY; for details type `$0 -l'.\n";
    print STDERR "This is free software, and you are welcome to redistribute it under certain conditions;\n";
    print STDERR "type `$0 -l' for details.\n";
    print STDERR <<EOF;
NAME
    $0   parses PRIAM XML input (SBML level 3 version 1, from PRIAM) and generates an edge or adjacency list
    
USAGE
    $0 -i infile.xml -out_type edgelist --outfile outfile.edgelist

OPTIONS
    --infile          -i    Input file (SBML level 3 version 1)     REQUIRED
    --outfile         -o    Output file name                        REQUIRED
    --out_type              Output file format                      REQUIRED
                            (allowed options: edgelist or adjlist)
    --debug           -d    debug (INT).
    --help            -h    This help.
    --license         -l    License.

EOF
}


##############################################
#License
##############################################
sub license{
    print STDERR <<EOF;

Copyright (C) 2016 Renato Augusto Correa dos Santos
e-mail: diriano\@gmail.com; renatoacsantos\@gmail.com

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
