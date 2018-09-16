from pymatgen.ext.matproj import MPRester as MP
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as SpGA
import numpy as np

class Find_unique_sites:
  """ Return the atomic sites Wyckoff positions, to avoid

  """
  def __init__(self, MP_ID, ctr_atom):
    self.structure = MP().get_structures(MP_ID)[0]
    self.ctr_atom = ctr_atom
    self.ctr_atom_pos = self.structure.indices_from_symbol(ctr_atom)
    self.frac_coord = self.structure.frac_coords
    self.symm_finder = SpGA(self.structure)
    self.uniq_ctr_atom_coords = []
    self.uniq_ctr_atom_multiplicity = []
    self.uniq_ctr_atom_pos = []
    self.get_uniq_ctr_atom_coords()

# find the coordinates, multiplicity and postion of the unique absorbing atoms
  def get_uniq_ctr_atom_coords(self):
    symm_struct = self.symm_finder.get_symmetrized_structure()
    if not np.array_equal(self.frac_coord, symm_struct.frac_coords):
       print("The atomic arrangement changed after getting symmeitrized structure\n")
    equiv_site_seq = symm_struct.equivalent_sites
    for i in equiv_site_seq:
      if i[0].species_string == self.ctr_atom:
        self.uniq_ctr_atom_coords.append(i[0].frac_coords)
        self.uniq_ctr_atom_multiplicity.append(len(i))
    for i in self.uniq_ctr_atom_coords:
      for j in self.ctr_atom_pos:
        if np.array_equal(i, self.frac_coord[j,:]):
          self.uniq_ctr_atom_pos.append(j+1)
          break
