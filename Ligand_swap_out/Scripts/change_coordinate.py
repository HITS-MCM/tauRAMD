#!/usr/bin/env python

import os
import sys
import string
import re

def Usage():
    print ("Usage:  python ../change_charge.py  template.mol2 new_coordinates.pdb  new_name.mol2")

if(len(sys.argv)>2):
    try:
         file_input1=sys.argv[1]  # mol2 file with correct charges
         file_input2=sys.argv[2]  # pdb file to use as template and keep coorsinates
         file_input3=sys.argv[3]  # mol2 file for output
    except IndexError:
         Usage()
         sys.exit(1)       
else:
  Usage()
  sys.exit(1)

f1=open(file_input1,"r")
f2=open(file_input2,"r")
lines1 = f1.readlines()
lines2 = f2.readlines()
atom = []
coordinates = []
for l in lines2:
    s = re.split("\s+",l)
 #   print l,begin
    if (l.find("HETATM")>=0) or (l.find("ATOM")>=0):
      atom.append(s[2].strip())
      x = float(l[30:38].strip())   # atomic coordinatess = s.strip()
      y = float(l[38:46].strip())
      z = float(l[46:54].strip())
      coordinates.append([x,y,z])
f1.close()
f2.close()
#print charge

begin = 0
atom_number = 0
fout=open(file_input3,"w")
for l in lines1:
    s = re.split("\s+",l)
 #   print l,begin
    if (l.find("@<TRIPOS>ATOM")>=0):
      begin = 1
      fout.write(l)
    elif (l.find("@<TRIPOS>BOND")>=0):
      begin = 2
      fout.write(l)
    elif(begin == 1):
        if s[2].strip() == atom[atom_number]:
             str_coor = ("%9.4f %9.4f %9.4f" %(coordinates[atom_number][0],coordinates[atom_number][1],coordinates[atom_number][2]))
             l1=l[:17]+str_coor+l[46:]
             fout.write(l1)
             atom_number += 1
        else:
           print("Different atom names in mol2 and pdb files:",s[2],atom[atom_number])
    else:   
       fout.write(l)
fout.close()
