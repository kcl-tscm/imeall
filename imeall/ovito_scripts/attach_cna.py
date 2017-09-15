import argparse
from ovito.io import export_file, import_file
from ovito.modifiers import CommonNeighborAnalysisModifier

parser = argparse.ArgumentParser()
parser.add_argument('--input_file', '-i', help='Input grain boundary struct file to produce cumulative energy.')
args = parser.parse_args()

#initialize node and import file
node = import_file(args.input_file)
cna = CommonNeighborAnalysisModifier()
node.modifiers.append(cna)
node.compute()
export_file(node, "output.xyz", "xyz", columns = ["Particle Identifier", "Particle Type",
																									"Position.X", "Position.Y", "Position.Z",
                                                  "Potential Energy", "Structure Type"])
