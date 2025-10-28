# WGS - Proteins - Blast - MW - G%A%S%T% - GA,GC,GT - GAGAGS - Structure

import os
import requests
from modeller import *
import argparse  
import warnings
import os
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments.')
parser.add_argument("-f", "--input", type=str, required=True, help="Input: Enter your WGS file name (It should be in FASTA format with '.fasta' or '.fna' extension)")
parser.add_argument("-e","--evalue", type=float, default=0.01, help="E-value: E-value for BLAST in structure prediction Default = 0.01")
parser.add_argument("-i","--ident", type=float, default=0.01, help="Identity: Minimum Pecentage of Identity (0-100) Default = 60")
parser.add_argument("-a","--alen", type=int, default=0.01, help="Alignment-Length: Minimum Alignment Length Default = 100")

args = parser.parse_args()
wgsfile = args.input
evalue = args.evalue
pident = args.ident
aln_lngth = args.alen
if ".fasta" in wgsfile or ".fna" in wgsfile:
    pass
else:
    print("Invalid input file name.")
    exit()
name = wgsfile.replace(".fasta","").replace(".fna","")

if not os.path.exists(f'./{name}_output'):
    os.mkdir(f'./{name}_output')
# os.system(f"augustus --species=tribolium  {wgsfile} > ./{name}_output/{name}.gff")
print("***Augustus Run Completed***")
with open(f"./{name}_output/{name}.gff","r") as r:
    data = r.readlines()
start = []
end = []
for i in range(len(data)):
    if "# protein sequence =" in data[i]:
        start.append(i)
    if "# Evidence" in data[i]:
        end.append(i)
if len(end)==0:
    for i in range(len(data)):
        if "# end gene" in data[i]:
            end.append(i)
wh = open(f"./{name}_output/{name}_protein.faa","w")
n = 0
for i in range(len(start)):
    n = n+1
    for j in range(start[i],end[i]):
        a = data[j].replace("# ","").replace("protein sequence = [",f">Protein_{n}\n").replace("]","")
        wh.write(a)
wh.close()
print("***GFF file is converted to FASTA file***")
os.system(f"diamond blastp --query ./{name}_output/{name}_protein.faa --db ./fibroin_db/fibroin_db.dmnd -o ./{name}_output/{name}_fibroin_blast.tsv --evalue {evalue} --quiet")
print("***BLASTP is DONE***")

with open(f"./{name}_output/{name}_fibroin_blast.tsv","r") as r:
    data = r.readlines()
Prot = []
ident = []
aln = []
for dat in data:
    a = dat.split("\t")
    Prot.append(a[0].replace("\n",""))
    ident.append(float(a[2]))
    aln.append(float(a[3]))

carote_ser_blast = []
for i in range(len(ident)):
    if ident[i] > pident and aln[i] > aln_lngth:
        
        carote_ser_blast.append(f"{Prot[i]}\n")
        
carote_ser_blast = set(carote_ser_blast)
carote_ser_blast = list(carote_ser_blast)

with open(f"./{name}_output/{name}_protein.faa","r") as r1:
    data1 = r1.read()

wh = open(f"./{name}_output/{name}_fibroin_blast.fasta","w")
sequence = data1.split(">")
for acc in carote_ser_blast:
    for seq in sequence:
        if acc in seq:
            wh.write(f">{seq}")
wh.close()
print("***BLAST filteration is DONE***")

wh_prpty = open(f"./{name}_output/{name}_fibroin_properties.tsv", 'w')
wh_fibroin = open(f"./{name}_output/{name}_fibroin_protein.fasta", 'w')
# WGS - Proteins - Blast - MW - G%A%S%T% - GA,GC,GT - GAGAGS - Structure
wh_prpty.write("Protein_name\tMolecular_weight\tAliphatic_index\tInstability_index\tIso-electric_point\tGRAVY\tG%\tA%\tS%\tGA\tGS\tGT\tGAGAGS\n")
try:
    with open(f"./{name}_output/{name}_fibroin_blast.fasta", 'r') as fh:
        data = fh.readlines()
except IOError:
    print("Unable to open the file. Try again.")
    exit()
if(">" not in data[0]):
        print("Invalid input\nMissing '>' from your fasta file")
        exit()
acc=[]
for i in data:
    if(">" in i):
        a = i.replace(">","").replace("\n","")
        acc.append(f"{a}")
pos=[]
for i in range(len(data)):
    if(">" in data[i]):
        pos.append(i)
pos.append(len(data))

seq=[]

for i in range(len(pos)-1):   
    for j in range(pos[i]+1,pos[i+1]):
        seq.append(data[j])
    seq.append("\t")

protein_list=[]
gr=[]
for i in seq:
    if i != "\t":
        gr.append(i)
    else:
        protein_list.append("".join(gr).replace(' ', '').replace('\n', ''))
        gr=[]
if gr:
    protein_list.append("".join(gr).replace(' ', '').replace('\n', ''))
mw = {
      'A': 71.08, 'R': 156.20, 'N': 114.11, 'D': 115.09, 'C': 103.14,
      'E': 129.12, 'Q': 128.41, 'G': 57.06, 'H': 137.15, 'I': 113.17,
      'L': 113.17, 'K': 128.1741, 'M': 131.21, 'F': 147.18, 'P': 97.12,
      'S': 87.08, 'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.14
    }
hyd = {
    'A': 1.8,  # Alanine
    'R': -4.5, # Arginine
    'N': -3.5, # Asparagine
    'D': -3.5, # Aspartic acid
    'C': 2.5,  # Cysteine
    'Q': -3.5, # Glutamine
    'E': -3.5, # Glutamic acid
    'G': -0.4, # Glycine
    'H': -3.2, # Histidine
    'I': 4.5,  # Isoleucine
    'L': 3.8,  # Leucine
    'K': -3.9, # Lysine
    'M': 1.9,  # Methionine
    'F': 2.8,  # Phenylalanine
    'P': -1.6, # Proline
    'S': -0.8, # Serine
    'T': -0.7, # Threonine
    'W': -0.9, # Tryptophan
    'Y': -1.3, # Tyrosine
    'V': 4.2   # Valine
}
for ac in range(len(acc)):
    plen=len(protein_list[ac])
    add=0
    #molecular weight
    for b in range(0,plen):
        base=protein_list[ac][b]
        try:
            add=add+mw[f"{base}"]
        except(KeyError):
            add=add
    mwprot=(add+18.01524)/1000

    #aliphatic index
    a = protein_list[ac].count("A")
    i = protein_list[ac].count("I")
    l = protein_list[ac].count("L")
    v = protein_list[ac].count("V")

    MolPA=a/plen*100
    MolPI=i/plen*100
    MolPL=l/plen*100
    MolPV=v/plen*100
    rvV=2.9
    rvIL=3.9
    AI=MolPA + (rvV*MolPV) + (rvIL*(MolPI+MolPL))
    
    #Instability index
    out =[]
    for i in range(0,plen-1):
        try:
            if(protein_list[ac][i]=="W"):
                value={ 'W': 1.0, 'C': 1.0, 'M': 24.68, 'H': 24.68, 'Y': 1.0, 'F': 1.0, 'Q': 1.0, 'N': 13.34, 'I': 1.0, 'R': 1.0, 'D': 1.0, 'P': 1.0, 'T': -14.03, 'K': 1.0, 'E': 1.0, 'V': -7.49, 'S': 1.0, 'G': -9.37, 'A': -14.03, 'L': 13.34 }
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="C"):
                value={ 'W': 24.68, 'C': 1.0, 'M': 33.6, 'H': 33.6, 'Y': 1.0, 'F': 1.0, 'Q': -6.54, 'N': 1.0, 'I': 1.0, 'R': 1.0, 'D': 20.26, 'P': 20.26, 'T': 33.6, 'K': 1.0, 'E': 1.0, 'V': -6.54, 'S': 1.0, 'G': 1.0, 'A': 1.0, 'L': 20.26 }
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="M"):
                value={'W':1.0,'C':1.0,'M':-1.88,'H':58.28,'Y':24.68,'F':1.0,'Q':-6.54,'N':1.0,'I':1.0,'R':-6.54,'D':1.0,'P':44.94,'T':-1.88,'K':1.0,'E':1.0,'V':1.0,'S':44.94,'G':1.0,'A':13.34,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="H"):
                value={'W':-1.88,'C':1.0,'M':1.0,'H':1.0,'Y':44.94,'F':-9.37,'Q':1.0,'N':24.68,'I':44.94,'R':1.0,'D':1.0,'P':-1.88,'T':-6.54,'K':24.68,'E':1.0,'V':1.0,'S':1.0,'G':-9.37,'A':1.0,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="Y"):
                value={'W':-9.37,'C':1.0,'M':44.94,'H':13.34,'Y':13.34,'F':1.0,'Q':1.0,'N':1.0,'I':1.0,'R':-15.91,'D':24.68,'P':13.34,'T':-7.49,'K':1.0,'E':-6.54,'V':1.0,'S':1.0,'G':-7.49,'A':24.68,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="F"):
                value={'W':1.0,'C':1.0,'M':1.0,'H':1.0,'Y':33.6,'F':1.0,'Q':1.0,'N':1.0,'I':1.0,'R':1.0,'D':13.34,'P':20.26,'T':1.0,'K':-14.03,'E':1.0,'V':1.0,'S':1.0,'G':1.0,'A':1.0,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="Q"):
                value={'W':1.0,'C':-6.54,'M':1.0,'H':1.0,'Y':-6.54,'F':-6.54,'Q':20.26,'N':1.0,'I':1.0,'R':1.0,'D':20.26,'P':20.26,'T':1.0,'K':1.0,'E':20.26,'V':-6.54,'S':44.94,'G':1.0,'A':1.0,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="N"):
                value={'W':-9.37,'C':-1.88,'M':1.0,'H':1.0,'Y':1.0,'F':-14.03,'Q':-6.54,'N':1.0,'I':44.94,'R':1.0,'D':1.0,'P':-1.88,'T':-7.49,'K':24.68,'E':1.0,'V':1.0,'S':1.0,'G':-14.03,'A':1.0,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="I"):
                value={'W':1.0,'C':1.0,'M':1.0,'H':13.34,'Y':1.0,'F':1.0,'Q':1.0,'N':1.0,'I':1.0,'R':1.0,'D':1.0,'P':-1.88,'T':1.0,'K':-7.49,'E':44.94,'V':-7.49,'S':1.0,'G':1.0,'A':1.0,'L':20.26}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="R"):
                value={'W':58.28,'C':1.0,'M':1.0,'H':20.26,'Y':-6.54,'F':1.0,'Q':20.26,'N':13.34,'I':1.0,'R':58.28,'D':1.0,'P':20.26,'T':1.0,'K':1.0,'E':1.0,'V':1.0,'S':44.94,'G':-7.49,'A':1.0,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="D"):
                value={'W':1.0,'C':1.0,'M':1.0,'H':1.0,'Y':1.0,'F':-6.54,'Q':1.0,'N':1.0,'I':1.0,'R':-6.54,'D':1.0,'P':1.0,'T':-14.03,'K':-7.49,'E':1.0,'V':1.0,'S':20.26,'G':1.0,'A':1.0,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="P"):
                value={'W':-1.88,'C':-6.54,'M':-6.54,'H':1.0,'Y':1.0,'F':20.26,'Q':20.26,'N':1.0,'I':1.0,'R':-6.54,'D':-6.54,'P':20.26,'T':1.0,'K':1.0,'E':18.38,'V':20.26,'S':20.26,'G':1.0,'A':20.26,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="T"):
                value={'W':-14.03,'C':1.0,'M':1.0,'H':1.0,'Y':1.0,'F':13.34,'Q':-6.54,'N':-14.03,'I':1.0,'R':1.0,'D':1.0,'P':1.0,'T':1.0,'K':1.0,'E':20.26,'V':1.0,'S':1.0,'G':-7.49,'A':1.0,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="K"):
                value={'W':1.0,'C':1.0,'M':33.6,'H':1.0,'Y':1.0,'F':1.0,'Q':24.68,'N':1.0,'I':-7.49,'R':33.6,'D':1.0,'P':-6.54,'T':1.0,'K':1.0,'E':1.0,'V':-7.49,'S':1.0,'G':-7.49,'A':1.0,'L':-7.49}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="E"):
                value={'W':-14.03,'C':44.94,'M':1.0,'H':-6.54,'Y':1.0,'F':1.0,'Q':20.26,'N':1.0,'I':20.26,'R':1.0,'D':20.26,'P':20.26,'T':1.0,'K':1.0,'E':33.6,'V':1.0,'S':20.26,'G':1.0,'A':1.0,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="V"):
                value={'W':1.0,'C':1.0,'M':1.0,'H':1.0,'Y':-6.54,'F':1.0,'Q':1.0,'N':1.0,'I':1.0,'R':1.0,'D':-14.03,'P':20.26,'T':-7.49,'K':-1.88,'E':1.0,'V':1.0,'S':1.0,'G':-7.49,'A':1.0,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="S"):
                value={'W':1.0,'C':33.6,'M':1.0,'H':1.0,'Y':1.0,'F':1.0,'Q':20.26,'N':1.0,'I':1.0,'R':20.26,'D':1.0,'P':44.94,'T':1.0,'K':1.0,'E':20.26,'V':1.0,'S':20.26,'G':1.0,'A':1.0,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="G"):
                value={'W':13.34,'C':1.0,'M':1.0,'H':1.0,'Y':-7.49,'F':1.0,'Q':1.0,'N':-7.49,'I':-7.49,'R':1.0,'D':1.0,'P':1.0,'T':-7.49,'K':-7.49,'E':-6.54,'V':1.0,'S':1.0,'G':13.34,'A':-7.49,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="A"):
                value={'W':1.0,'C':44.94,'M':1.0,'H':-7.49,'Y':1.0,'F':1.0,'Q':1.0,'N':1.0,'I':1.0, 'R':1.0,'D':-7.49,'P':20.26,'T':1.0,'K':1.0,'E':1.0,'V':1.0,'S':1.0,'G':1.0,'A':1.0,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
            elif(protein_list[ac][i]=="L"):
                value={'W':24.68,'C':1.0,'M':1.0,'H':1.0,'Y':1.0,'F':1.0,'Q':33.6,'N':1.0,'I':1.0,'R':20.26,'D':1.0,'P':20.26,'T':1.0,'K':-7.49,'E':1.0,'V':1.0,'S':1.0,'G':1.0,'A':1.0,'L':1.0}
                out.append(value[f"{protein_list[ac][i+1]}"])
        except(KeyError):
            out.append(0)
    add=sum(out)
    II=(10/plen)*add

    #Iso-electric Point
    d=e=c=y=h=k=r=u=0
    for aa in protein_list[ac]:
        if(aa == 'D'):
            d=d+1
        elif(aa == 'E'):
            e=e+1
        elif(aa == 'C'):
            c=c+1
        elif(aa == 'Y'):
            y=y+1
        elif(aa == 'H'):
            h=h+1
        elif(aa == 'K'):
            k=k+1
        elif(aa == 'R'):
            r=r+1
        else:
            u=u+1

    ph=0
    end=0
    eqd=eqe=eqc=eqy=eqh=eqk=eqr=eqcooh=eqnh=0.0
    while(end==0):
        eqd=-d/(1+(10**(3.65-ph)))
        eqe=-e/(1+(10**(4.25-ph)))
        eqc=-c/(1+(10**(8.18-ph)))
        eqy=-y/(1+(10**(10.7-ph)))
        eqh=h/(1+(10**(ph-6.0)))
        eqk=k/(1+(10**(ph-10.53)))
        eqr=r/(1+(10**(ph-12.48)))
        eqcooh=-1/(1+(10**(2.34-ph)))
        eqnh=1/(1+(10**(ph-9.69)))
        tpi=eqd+eqe+eqc+eqy+eqh+eqk+eqr+eqcooh+eqnh
        if(tpi>0):
            ph=ph+0.01
        else:
            end+=1

    iso_point = ph
    #gravy
    grv = 0
    for z in range(0,plen):
        base=protein_list[ac][z]
        try:
            grv=grv+hyd[f"{base}"]
        except(KeyError):
            grv=grv
    gravy=grv/plen
    # WGS - Proteins - Blast - MW - G%A%S%T% -  -  - Structure
    #G%A%S%T%
    G_ = protein_list[ac].count("G")*100/len(protein_list[ac])
    A_ = protein_list[ac].count("A")*100/len(protein_list[ac])
    S_ = protein_list[ac].count("S")*100/len(protein_list[ac])
    

    # GA,GS,GT
    GA = protein_list[ac].count("GA")
    GS = protein_list[ac].count("GS")
    GT = protein_list[ac].count("GT")

    #GAGAGS
    GAGAGS = protein_list[ac].count("GAGAGS")
    wh_prpty.write(f"{acc[ac]}\t{mwprot}\t{AI}\t{II}\t{iso_point}\t{gravy}\t{G_}\t{A_}\t{S_}\t{GA}\t{GS}\t{GT}\t{GAGAGS}\n")
    if GAGAGS > 1 and G_ > 35 and A_ > 25 and S_ > 8 and GA > 1 and GS > 1 and GT > 1:
        wh_fibroin.write(f">Hevay_Chain_{acc[ac]}\n{protein_list[ac]}\n")
    elif mwprot>20 and mwprot<30:
        wh_fibroin.write(f">Light_Chain_{acc[ac]}\n{protein_list[ac]}\n")

wh_fibroin.close()
wh_prpty.close()


os.system(f"cp ./aln_model.py ./{name}_output")
os.system(f"cp ./build_model.py ./{name}_output")
try:
    with open(f"./{name}_output/{name}_fibroin_protein.fasta", 'r') as fh:
        data = fh.readlines()
except IOError:
    print("Unable to open the file. Try again.")
    exit()
if(">" not in data[0]):
        print("Invalid input\nMissing '>' from your fasta file")
        exit()
accs=[]
for i in data:
    if(">" in i):
        a = i.replace(">","").replace("\n","").replace("|","_").replace(" ","_")
        if len(a)>30:
            accs.append(a[:30])
        else:
            accs.append(a)

pos=[]
for i in range(len(data)):
    if(">" in data[i]):
        pos.append(i)
pos.append(len(data))

seq=[]

for i in range(len(pos)-1):   
    for j in range(pos[i]+1,pos[i+1]):
        seq.append(data[j])
    seq.append("\t")

protein_list=[]
gr=[]
for i in seq:
    if i != "\t":
        gr.append(i)
    else:
        protein_list.append("".join(gr).replace(' ', '').replace('\n', ''))
        gr=[]
if gr:
    protein_list.append("".join(gr).replace(' ', '').replace('\n', ''))

i = 0

for i in range(len(accs)):

    with open(f"./{name}_output/{accs[i]}.fasta","w") as wh:
        wh.write(f">{accs[i]}\n{protein_list[i]}")
    wrt = open(f"./{name}_output/{accs[i]}.ali","w")
    wrt.write(f">P1;{accs[i]}\nsequence:{accs[i]}:::::::0.00: 0.00\n{protein_list[i]}*")
    wrt.close()
    
    os.system(f"diamond blastp --query ./{name}_output/{accs[i]}.fasta --db ./fibroin_db/fibroin_db.dmnd -o ./{name}_output/{accs[i]}_blast.tsv --evalue {evalue} --quiet")
    


for i in range(len(accs)):
    
    print(os.getcwd())
    print(accs[i])
    try:
        with open(f'./{name}_output/{accs[i]}_blast.tsv', 'r') as fh:
            data = fh.readlines()
    except IOError:
        print("Unable to open the tsv file. Try again.")
        continue
    ac = []
    ident = []
    for m in range(len(data)):
        a = data[m].split("\n")
        b = a[0].split("\t")
        ident.append(float(b[2]))
        try:
            c = b[1].split('|')
            ac.append(c[1])
        except(IndexError):
            c = b[1].split('_')
            ac.append(c[1])

    acc = []
    try:
        if max(ident) > pident:
            for cut in range(len(ident)):

                if ident[cut] > pident:
                    acc.append(ac[cut])
        else:
            err = open(f"./{name}_output/{accs[i]}_empty.txt","w")
            err.write(f"No templates avialable for {accs[i]} when percent of identity is {pident}")
            err.close()
            continue
    except(ValueError):
        err = open(f"./{name}_output/{accs[i]}_empty.txt","w")
        err.write(f"No templates avialable for {accs[i]} when percent of identity is {pident}")
        err.close()
        continue
    # alph=[]
    print(acc)
    uni =[]
    templates = []
    for id in range(len(acc)):
        alph = f"https://alphafold.ebi.ac.uk/files/AF-{acc[id]}-F1-model_v6.pdb"
        a=requests.get(alph)
        f2=a.text
        if "<Error>" not in f2:
            f=open(f"./{name}_output/{accs[i]}_{acc[id]}_template.pdb","w")
            f.write(f2)
            f.close()
            templates.append(f"{accs[i]}_{acc[id]}_template")
        
    if len(templates) == 0:
        err = open(f"./{name}_output/{accs[i]}_empty.txt","w")
        err.write(f"No templates avialable for {accs[i]} when percent of identity is {pident}")
        err.close()
        continue
        
  
    jn = " ".join(templates)
    os.chdir(f"./{name}_output")
    os.system(f"python aln_model.py -ali {accs[i]} -tmp {jn}")
    os.system(f"python build_model.py -ali {accs[i]}-multiple_templates -tmp {jn}")

    os.chdir("../")
