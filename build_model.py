from modeller import *
import os
import argparse 
import warnings
from modeller.automodel import *

# from RamachanDraw import phi_psi, plot

warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments.') 
parser.add_argument("-ali", "--alignment", type=str, required=True)
parser.add_argument("-tmp","--template", nargs='+', default= ["all_tools"])
args = parser.parse_args()
env = Environ()

templates = args.template

seq_code = args.alignment

name1 = seq_code.replace("-multiple_templates","")
print(os. getcwd())
a = AutoModel(env,
              alnfile=f"{seq_code}.ali",
              knowns=tuple(templates),
              sequence=name1,
              assess_methods=(assess.DOPE,
                              # soap_protein_od.Scorer(),
                              assess.GA341))

a.starting_model = 1
a.ending_model = 3
a.make()
# def ramch(accession):

#     ax = plot(f"./{accession}.pdb")
   
    
#     plt.savefig(f"./{accession}_ramachandran.png")

#     pdb_file = f"./{accession}.pdb"

#     phi_psi_dict = phi_psi(pdb_file)
#     aa = phi_psi_dict.keys()
#     aa=list(aa)
#     ph_ps = phi_psi_dict.values()
#     ph_ps=list(ph_ps)

#     with open(f"./{accession}_phi_psi.tsv","w") as wh:
#         wh.write("Amino_acids\tphi_values\tpsi_values\n")
#         for l in range(len(aa)):
#             wh.write(f"{aa[l]}\t{ph_ps[l][0]}\t{ph_ps[l][1]}\n")

with open(f"{name1}_dope_scores.txt", "w") as f:
    f.write("Model_Name\tDOPE_Score\n")
    for mdl in a.outputs:
        name = mdl['name']
        dope_score = mdl['DOPE score']
        f.write(f"{name}\t{dope_score}\n")
        print(f"{name}: DOPE = {dope_score}")
        n = name.replace(".pdb","")
        # ramch(n)
