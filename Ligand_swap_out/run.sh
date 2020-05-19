#!/bin/sh

#---------------- Input data - 
Structure='Ligands1/RdRp-APVD-173.pdb '
Ligan_frcmod='Ligands1/A17_resp.frcmod'
Ligand_mol2='Ligands1/A17_resp.mol2'
Ligand_resi_name='A17'


# MD parameters
MD_prmtop='MD/input.prmtop'
MD_rst7='MD/equil.rst7'
MD_resi_name='A20'  #redidue name of the ligand to be exchanged

# generate MD.pdb from rest7
echo "============= generate equilibrated pdb from rest7 and prmtop ============="
cp $MD_prmtop MD.prmtop
cp $MD_rst7 MD.rst7
cpptraj < Scripts/cpptraj.in

echo "============= save BOX parameters ======================================="
grep --after-context=2 BOX  $MD_prmtop > MD_box.dat

# Superimpose 
echo "=============== superimpose ============================="
grep $Ligand_resi_name $Structure > Ligand.pdb
python  Scripts/superimpose_swapp.py  MD.pdb  $Structure Ligand.pdb

echo "================ change coordinates in the mol2 file ===================="
cp $Ligand_mol2 INH.mol2
python Scripts/change_coordinate.py  INH.mol2  Ligand_superimposed.pdb   Ligand_superimposed.mol2

#p generate ligand topology
echo "==================  generate ligand topology============================"
cp $Ligan_frcmod INH.frcmod
tleap -f Scripts/tleap_ligand_in


echo "====================== creat new coordinates and topology===================="
# Generation new PRMTOP and Coordinates
python Scripts/parmed_swapp.py MD.prmtop MD.rst7 $MD_resi_name  INH.prmtop INH.inpcrd  $Ligand_resi_name.prmtop $Ligand_resi_name.rst7


