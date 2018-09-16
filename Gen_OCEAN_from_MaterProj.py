import sys
import os
import shutil
import numpy as np
from sys import argv, exit
from pymatgen.ext.matproj import MPRester as MP
from MPEntry import MPEntry

""" Main function to generate input files for X-ray core-level spectra computation.
    Compounds to compute are read from stdin.
    Interfaces with Materials Project (https://materialsproject.org/) to extract
    useful material informtion obtained from first-principles calculations.
    For polymorphic strcutrues, the ones with zero formation energy are selected.
    Heavy work is done by the helper function MPEntry.

    Arguments:
    K_broad : dict {str : float}
      K-edge broadening.

    U_eff : dict {str : float}
      Hubbard U.

    Mater2Do : list [str]
      Formulas of compounds to be computed (can be more than 1).
      Can be input interactively.

    MP_IDs : list [str]
      Materials project ID of compounds (formation eng == 0) to be computed.

    ctr_atom :

    all_compounds_info :
"""

here = os.getcwd()
sys.path.append(here)
ocean_tmplate_loc = here

# Broadening parameters for K-edge
K_broad = {'Li': 0.5, 'O': 0.55, 'K': 0.74, 'Ca': 0.81, 'Sc': 0.86, 'Ti': 0.94, 'V': 1.01, 'Cr': 1.08, 'Mn': 1.16, 'Fe': 1.25, 'Co': 1.33, 'Cu': 1.55}

# U_eff paramters for DFT+U
U_eff = {'Ti': 4.4, 'V': 2.7, 'Cr': 3.5, 'Mn': 4.0, 'Fe': 4.6, 'Co': 5.0, 'Ni': 5.1, 'Cu': 4.0, 'Zn': 7.5, 'Ga': 3.9, 'Nb': 2.1, 'Mo': 2.4, 'Tc': 2.7, 'Ru': 3.0, 'Rh': 3.3, 'Pd': 3.6, 'Cd': 2.1, 'In': 1.9, 'Ta': 2.0, 'W': 2.2, 'Re': 2.4, 'Os': 2.6, 'Ir': 2.8, 'Pt': 3.0, 'La': 7.5, 'Ce': 6.3, 'Pr': 5.5, 'Gd': 6.0, 'Nd': 6.2, 'Sm': 6.4, 'Eu': 5.4, 'Tm': 6.0, 'Yb': 6.3, 'Lu': 3.8, 'Sc': 2.9, 'Sn': 3.5}

# Compounds to be calculated, taken from command line or keyboard
if len(argv) == 1:
  from_keyboard = input("please specify the compounds/elements:\n")
  Mater2Do = from_keyboard.split()
else:
  Mater2Do = argv[1:]

ctr_atom = input("please specify the centered atom species:\n")

# Determine te MaterProj ID for the materials on the convex hull
MP_IDs = []
for formula in Mater2Do:
  Eng_ID = MP().get_data(formula, "vasp", "e_above_hull")
  for i in Eng_ID:
    # select materials with zero formation energy (residing at convex hull)
    if i['e_above_hull'] == 0:
      MP_IDs.append(i['material_id'])

# Get data from Materials Project database
all_compounds_info = []  # list of MPEntry class
for i in MP_IDs:
  curr_Mater = MPEntry(i, ctr_atom)
  all_compounds_info.append(curr_Mater)

# Write OCEAN input file
here = os.getcwd()
for i in all_compounds_info:
  os.chdir(here)
  new_folder = i.formula + '_' + i.ID
  new_file = i.formula + '.in'
  os.mkdir(new_folder)
  os.chdir(here + '/' + new_folder)
  shutil.copy(ocean_tmplate_loc + '/' + "OCEAN.head", new_file)
  fout = open(new_file, "a")
  i.write_to_file(fout, ctr_atom)
  fout.write("# Spectral broadening in eV\n")
  fout.write("cnbse.broaden   " + str(K_broad[ctr_atom]) + '\n')
  fout.write("#LDA+U for QE\nldau{\nlda_plus_u=.true. ")
  j = 1
  for j in range(0, i.ntypat):
    fout.write(",\nHubbard_U(" + str(j + 1) + ")=" + str(U_eff.get(i.elements[j], 0)))
  fout.write('\n}\n')
  fout.close()
