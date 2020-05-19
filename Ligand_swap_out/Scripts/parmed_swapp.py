import parmed as pmd
import sys, os

def Usage():
    print("Usage:\n")
    print("python parmed_swapp old_parmtop old_rst7 old_ligand_name ligand_parmtop ligand.pdb new_prmtop new_rst7")
    return

try:
    old_prmtop = sys.argv[1]
    old_rst7 =  sys.argv[2]
    name_old_ligand =  sys.argv[3]

    ligand_prmtop = sys.argv[4]
    ligand_crd =  sys.argv[5]

    new_prmtop = sys.argv[6]
    new_rst7 = sys.argv[7]
except:
    Usage()

amber = pmd.load_file(old_prmtop, old_rst7) 
o = amber['!:'+name_old_ligand]
s = pmd.load_file(ligand_prmtop, ligand_crd)
n = o+s
if os.path.exists(new_prmtop): os.system("mv "+new_prmtop+" "+new_prmtop+"_saved")
if os.path.exists(new_rst7): os.system("mv "+new_rst7+" "+new_rst7+"_saved")
n.save(new_prmtop)
n.save(new_rst7)

f=open("MD_box.dat")
lines=f.readlines()
l = "sed \'/%FORMAT(5E16.8)/{n;s/.*/" +lines[-1][0:-1]+ "/}\' "+new_prmtop+"  > "+new_prmtop[:-7]+"_BOX.prmtop"
#print(l)
os.system(l)


