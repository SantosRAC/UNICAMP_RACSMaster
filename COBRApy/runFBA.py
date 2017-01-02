import sys
import cobra
from cobra import Model, Reaction, Metabolite
import argparse

#Making sure you are running a version of python that works with this script.
if sys.version_info[0] != 2 or sys.version_info[1] < 7 or sys.version_info[2] < 8:
    print("This script requires Python version 2.7.8")
    sys.exit(1)

parser = argparse.ArgumentParser(description='Flux Balance Analysis', add_help=True)
parser.add_argument('-m','--model', dest='model', metavar='importedModel', type=str, help='Model (SBML3)', required=True)

# Get information from the argparse (arguments)
args = parser.parse_args()
impModel = cobra.io.read_sbml_model(args.model)

# Print number of reactions, metabolites, and genes in model
print("Number of Reactions: %d" % len(impModel.reactions))
print("Number of Metabolites: %d" % len(impModel.metabolites))
print("Number of Genes: %d" % len(impModel.genes))

# Print reactions and corresponding metabolites
for reaction in impModel.reactions:
 for metabolite in reaction.metabolites:
  print("%s\t%s" % (reaction, metabolite))
