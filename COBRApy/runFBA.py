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
parser.add_argument('-t','--model_type', default='sbml', dest='fileformat', metavar='formatFile', type=str, help='Model file format (sbml, mat)', required=False)
parser.add_argument('--optimize', dest='optimizeBoolean', action='store_true', help='Turn on optimization (GLPK)', required=False)

args = parser.parse_args()

# Get information from the argparse (arguments)
if args.fileformat == 'sbml':
 impModel = cobra.io.read_sbml_model(args.model)
elif args.fileformat == 'mat':
 impModel = cobra.io.load_matlab_model(args.model)
else:
 print("Something wrong...")
 exit(1)

# Print number of reactions, metabolites, and genes in model
print("#### Summary about the imported model ####")
print("Number of Reactions: %d" % len(impModel.reactions))
print("Number of Metabolites: %d" % len(impModel.metabolites))
print("Number of Genes: %d" % len(impModel.genes))

# Run Optimize for the current objective function
if args.optimizeBoolean:
 impModel.optimize()
 print("### Optimization (objective function)")
 #print(impModel.summary())
 print(impModel.objective)
 print(impModel.solution.status)
 print(impModel.solution.f)
 for metabolite in impModel.metabolites:
  print(metabolite.summary())

#Print details about reactions and their metabolites (reactants and products)
#for reaction in impModel.reactions:
# print("Reaction catalyzed: %s" % (reaction.reaction))
# print("Reaction identifier (id): %s" % (reaction))
# print("Reaction lower bound: %s" % (reaction.lower_bound))
# print("Reaction upper bound: %s" % (reaction.upper_bound))
# print("Reaction reversibility: %s" % (reaction.reversibility))
# print("Reaction Check Mass Balance: %s" % (reaction.check_mass_balance()))
# print("Reaction Gene Reaction Rule: %s" % (reaction.gene_reaction_rule))
# print("Reaction Associated Genes: %s" % (reaction.genes))
# for metabolite in reaction.metabolites:
#  print("Metabolite identifier (id): %s" % (metabolite))
#  print("Metabolite name: %s" % (metabolite.name))
#  print("Metabolite compartment: %s" % (metabolite.compartment))
#  print("Metabolite charge: %s" % (metabolite.charge))

# Print reactions in which each metabolite is present
#for metabolite in impModel.metabolites:
# print(metabolite.reactions)

# Print all genes in model
#for gene in impModel.genes:
# print(gene)
