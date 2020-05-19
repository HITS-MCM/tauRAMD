import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import contacts,align,rms
import sys,os,glob


def Usage():
    print("Usage: python superimpose_swapp.py  MD_pdb  Struct_pdb [ligands_dir]")
    print("       ligands_dir - a directory with a set of ligand file in PDB")
    return


if len(sys.argv) >2:
    file_MD =sys.argv[1]
    file_structure = sys.argv[2]
else:
    Usage()
    sys.exit()

u = mda.Universe(file_structure)
uMD = mda.Universe(file_MD)
ul = []
fl = []
if  len(sys.argv) > 3:
    if sys.argv[3].find('pdb')< 0:
        d = glob.glob(sys.argv[3]+"/*pdb")
    else:
        d=[sys.argv[3]]
    for file_ligand in d:
        if os.path.exists(file_ligand):
            print("Ligand "+file_ligand+"  will be superimposed")
            ul.append(mda.Universe(file_ligand))
            fl.append(file_ligand)

ref_CA = uMD.select_atoms("name CA")
ref0 = ref_CA.positions - ref_CA.center_of_mass()
u_CA = u.select_atoms("name CA")
for ul_i in ul:  ul_i.atoms.translate(-u_CA.center_of_mass())
u.atoms.translate(-u_CA.center_of_mass())

u0 =  u_CA.positions - u_CA.center_of_mass()
R, rmsd = align.rotation_matrix(u0,ref0)  # compute rotation matrix

for i,ul_i in enumerate(ul): 
    ul_i.atoms.rotate(R)
    ul_i.atoms.translate(ref_CA.center_of_mass())
    ligand = ul_i.select_atoms("all")
    ligand = ul_i.select_atoms("all")
    ligand.write(fl[i][:-4]+"_superimposed.pdb")

u.atoms.rotate(R)
u.atoms.translate(ref_CA.center_of_mass()) # translate back to the old center of mass position
protein = u.select_atoms("all")
protein.write(file_structure[:-4]+"_superimposed.pdb")



