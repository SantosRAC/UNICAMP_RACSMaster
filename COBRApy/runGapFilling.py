import sys
import cobra
from cobra import Model, Reaction, Metabolite
import argparse
import cobra.test

#Making sure you are running a version of python that works with this script.
if sys.version_info[0] != 2 or sys.version_info[1] < 7 or sys.version_info[2] < 8:
    print("This script requires Python version 2.7.8")
    sys.exit(1)

parser = argparse.ArgumentParser(description='Flux Balance Analysis', add_help=True)
parser.add_argument('-m1','--model1', dest='model1', metavar='importedModel1', type=str, help='Model (SBML3 or MAT)', required=True)
parser.add_argument('-m2','--model2', dest='model2', metavar='importedModel2', type=str, help='Model (SBML3 or MAT)', required=True)
parser.add_argument('-t1','--model1_type', default='sbml', dest='fileformat1', metavar='formatFile1', type=str, help='Model file format (sbml, mat)', required=False)
parser.add_argument('-t2','--model2_type', default='sbml', dest='fileformat2', metavar='formatFile2', type=str, help='Model file format (sbml, mat)', required=False)
parser.add_argument('--gapfill', dest='gapFillBoolean', action='store_true', help='', required=False)

args = parser.parse_args()

# Choose gurobi as solver (better for MILP)
solver = cobra.solvers.gurobi_solver

# Get information from the argparse (arguments) about type of first model
if args.fileformat1 == 'sbml':
 impModel1 = cobra.io.read_sbml_model(args.model1)
elif args.fileformat1 == 'mat':
 impModel1 = cobra.io.load_matlab_model(args.model1)
else:
 print("Something wrong...")
 exit(1)

# Get information from the argparse (arguments) about type of second model
if args.fileformat2 == 'sbml':
 impModel2 = cobra.io.read_sbml_model(args.model2)
elif args.fileformat2 == 'mat':
 impModel2 = cobra.io.load_matlab_model(args.model2)
else:
 print("Something wrong...")
 exit(1)

# Run Optimize for the current objective function
if args.gapFillBoolean:
 r = cobra.flux_analysis.growMatch(impModel1, impModel2)
 for e in r[0]:
  print(e.id)
