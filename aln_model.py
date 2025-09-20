from modeller import *
import os
import argparse 
import warnings
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments.') 
parser.add_argument("-ali", "--alignment", type=str, required=True)
parser.add_argument("-tmp","--template", nargs='+', default= ["all_tools"])
args = parser.parse_args()
env = Environ()
aln = Alignment(env)

templates = args.template

ali = args.alignment
aln.append(file=f'./{ali}.ali', align_codes=f'{ali}')
# Append all 4 templates
for template_code in templates:
    # env.io.atom_files_directory = ['./structure_predict_out']
    mdl = Model(env, file=template_code, model_segment=('FIRST:A','LAST:A'))
    aln.append_model(mdl, align_codes=template_code, atom_files=f'{template_code}.pdb')

# Append target sequence
# aln.append(file=f'./{ali}.ali', align_codes=f'{ali}')

# Align target to templates
aln.align2d()

# Write alignment outputs
aln.write(file=f'{ali}-multiple_templates.ali', alignment_format='PIR')
aln.write(file=f'{ali}-multiple_templates.pap', alignment_format='PAP')
