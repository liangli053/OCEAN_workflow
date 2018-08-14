#!/Users/liangli/anaconda3/envs/my_pymatgen/bin/python3
import sys
import os
import shutil
import numpy as np
from sys import argv, exit
from pymatgen.ext.matproj import MPRester as MP
from MPEntry import MPEntry

here = os.getcwd()
sys.path.append(here)
ocean_tmplate_loc = here

##-----Broadening parameters for K-edge-----##
K_broad = {'Li':0.5, 'O':0.55, 'K':0.74, 'Ca':0.81, 'Sc':0.86, 'Ti':0.94, 'V':1.01, 'Cr':1.08, 'Mn':1.16, 'Fe':1.25, 'Co':1.33, 'Cu':1.55}

##-----U_eff paramters for DFT+U-----##
U_eff = {'Ti':4.4, 'V':2.7, 'Cr':3.5, 'Mn':4.0, 'Fe':4.6, 'Co':5.0, 'Ni':5.1, 'Cu':4.0, 'Zn':7.5, 'Ga':3.9, 'Nb':2.1, 'Mo':2.4, 'Tc':2.7, 'Ru':3.0, 'Rh':3.3, 'Pd':3.6, 'Cd':2.1, 'In':1.9, 'Ta':2.0, 'W':2.2, 'Re':2.4, 'Os':2.6, 'Ir':2.8, 'Pt':3.0, 'La':7.5, 'Ce':6.3, 'Pr':5.5, 'Gd':6.0, 'Nd':6.2, 'Sm':6.4, 'Eu':5.4, 'Tm':6.0, 'Yb':6.3, 'Lu':3.8, 'Sc':2.9, 'Sn':3.5}

##-----Compounds to be calculated, taken from command line or keyboard-----##
if len(argv) == 1:
  from_keyboard = input("please specify the compounds/elements:\n")
  Mater2Do = from_keyboard.split()
else:
  Mater2Do = argv[1:]

ctr_atom = input("please specify the centered atom species:\n")

##-----Reformat the input, change lowercase atomic symbols to uppercase-----##
for i in range(len(Mater2Do)):
  curr_Mater = list(Mater2Do[i])
  for j in range(len(curr_Mater)):
    if (j == 0) & (curr_Mater[j].islower()):
      curr_Mater[j]=curr_Mater[j].upper()
    elif (curr_Mater[j].islower()) & (curr_Mater[j-1].isdigit()):
      curr_Mater[j]=curr_Mater[j].upper()
  Mater2Do[i] ="".join(curr_Mater)

##-----Determine te MaterProj ID for the materials on the convex hull-----##
MP_IDs = []
for curr_Mater in Mater2Do:
  Eng_ID_4curr_Mater = MP().get_data(curr_Mater,"vasp","e_above_hull")
  for i in Eng_ID_4curr_Mater:
    if i['e_above_hull'] == 0:
      MP_IDs.append(i['material_id'])
##-----Get data from Materials Project Entries-----##
all_compounds_info = []
for i in MP_IDs:
  curr_Mater = MPEntry(i, ctr_atom)
  all_compounds_info.append(curr_Mater)

##-----Write OCEAN input file-----##
here = os.getcwd()
for i in  all_compounds_info:
  os.chdir(here)
  new_folder = i.formula+'_'+i.MP_ID
  new_file = i.formula+'.in'
  os.mkdir(new_folder)
  os.chdir(here+'/'+new_folder)
  shutil.copy(ocean_tmplate_loc+'/'+"OCEAN.head", new_file)
  fout = open(new_file, "a")
  i.write_to_file(fout, ctr_atom)
  fout.write("# Spectral broadening in eV\n")
  fout.write("cnbse.broaden   "+str(K_broad[ctr_atom])+'\n')
  fout.write("#LDA+U for QE\nldau{\nlda_plus_u=.true. ")
  j = 1
  for j in range(0,i.ntypat):
      fout.write(",\nHubbard_U("+str(j+1)+")="+str(U_eff.get(i.elements[j], 0)))
  fout.write('\n}\n')
  fout.close()
