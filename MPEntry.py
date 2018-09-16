from pymatgen.ext.matproj import MPRester as MP
from pymatgen.core.periodic_table import Element
import numpy as np
from math import pi, ceil
from Find_unique_sites import Find_unique_sites


class MPEntry:

  def __init__(self, MP_ID, ctr_atom):
    self.ID = MP_ID
    self.formula = MP().get_data(self.ID, "vasp", "pretty_formula")[0]['pretty_formula']
    self.ctr_atom = ctr_atom
    self.structure = MP().get_structures(MP_ID)[0]
    self.ntypat = self.structure.ntypesp
    self.elements = self.structure.symbol_set
    self.natom = self.structure.num_sites
    self.band_gap = MP().query(MP_ID, ['band_gap'])[0]['band_gap']
    self.total_formula = MP().query(MP_ID, ['unit_cell_formula'])[0]['unit_cell_formula']
    self.typat = []
    self.pp_list = []
    self.occopt = 1 if self.band_gap > 0 else 3
    self.fband = 0.4 if self.occopt == 1 else 2
    self.rprim = self.structure.lattice.matrix
    self.xred = self.structure.frac_coords
    self.a = float(self.structure.lattice.a)
    self.b = float(self.structure.lattice.b)
    self.c = float(self.structure.lattice.c)
    self.ngkpt = []
    self.get_typat()
    self.get_pp_list()
    self.get_ngkpt()
    self.edge_info = Find_unique_sites(self.ID, self.ctr_atom)
    self.print_spacegroup_info()

  def print_spacegroup_info(self):
    SP_symbol = MP().query(self.ID, ['spacegroup'])[0]['spacegroup']['symbol']
    SP_number = MP().query(self.ID, ['spacegroup'])[0]['spacegroup']['number']
    computed_SP_number = self.edge_info.symm_finder.get_space_group_number()
    computed_SP_symbol = self.edge_info.symm_finder.get_space_group_symbol()
    print("spacegroup of {} is: {:3d}  {:7s}\n".format(self.formula, SP_number, SP_symbol))
    print("spacegroup of Materials Project computed cell is: {:3d}  {:7s}\n".format(computed_SP_number, computed_SP_symbol))
    if SP_number != computed_SP_number:
      print("WARNING: The spacegroup of {} found in Materials Project POSCAR may be inaccurate,\n \
       this can cause additional computational load\n".format(self.formula))

  def get_znucl(self):
    znucl = []
    for i in self.elements:
      tmp = Element(i)
      znucl.append(tmp.Z)
    return znucl

  def get_typat(self):
    num_of_element = []
    for i in range(len(self.elements)):
      num_of_element.append(int(self.total_formula[self.elements[i]]))
      self.typat.extend([i + 1] * num_of_element[i])

  def get_pp_list(self):
    """ Return pseudopotential file names
    """
    pp_list = []
    for i in self.elements:
      tmp = Element(i)
      pp_list.append(str('{:02d}'.format(tmp.Z)) + '-' + i.lower() + ".lda.fhi")
    return pp_list

  def get_ngkpt(self):
    if self.band_gap > 0:
      KSPACING = 0.3
    else:
      KSPACING = 0.15
    self.ngkpt.append(ceil(2 * pi / self.a / KSPACING))
    self.ngkpt.append(ceil(2 * pi / self.b / KSPACING))
    self.ngkpt.append(ceil(2 * pi / self.c / KSPACING))

# write_to_fiile as a separate module
  def write_to_file(self, fout, ctr_atom):
    znucl = self.get_znucl()
    pp_list = self.get_pp_list()
    absorb_ctr = Element(ctr_atom)
    fout.write("# number of element types\n"
               "ntypat  " + str(self.ntypat) + "\n"
               "znucl  { " + " ".join(str(e) for e in znucl) + ' }\n'
               "# number of total atoms\n"
               "natom  " + str(self.natom) + "\n"
               "pp_list  { " + "  ".join(pp_list) + ' }\n'
               "# occupation, 1 or 3(smearing)\n"
               "occopt  " + str(self.occopt) + "\n"
               "fband  " + str(self.fband) + "\n"
               "# opf control files\n"
               "opf.fill  { " + str(absorb_ctr.Z) + " " + ctr_atom + ".fill }\n"
               "opf.opts  { " + str(absorb_ctr.Z) + " " + ctr_atom + ".opts }\n"
               "acell  { 1.889725989  1.889725989  1.889725989 }\n"
               "rprim  {\n")
    for line in self.rprim:
      fout.write("    " + " ".join(str("%16.12f" % e) for e in line) + '\n')
    fout.write('}\n'
               "typat  { " + " ".join(str(e) for e in self.typat) + ' }\n'
               "xred  {\n")
    for line in self.xred:
      fout.write("    " + " ".join(str("%16.12f" % e) for e in line) + '\n')
    fout.write('}\n'
               "# Kpt mesh for ground state density calculation\n"
               "ngkpt { " + " ".join(str(e) for e in self.ngkpt) + ' }\n'
               "# Kpt mesh for final states\n"
               "nkpt  { }\n "
               "# Kpt mesh for screening calculation\n"
               "screen.nkpt  { 2 2 2 }\n"
               "# Total bands for final states\n"
               "nbands\n"
               "# Total bands for screening calculation\n"
               "screen.nbands\n"
               "# edge information # number of edges to calculate # atom number, n, l\n"
               "nedges  " + str(len(self.edge_info.uniq_ctr_atom_pos)) + '\n'
               "edges {\n")
    for i, j in zip(self.edge_info.uniq_ctr_atom_pos, self.edge_info.uniq_ctr_atom_multiplicity):
      fout.write("     {:2d}   1   0   # multiplicity: {:2d}\n".format(i, j))
    fout.write('}\n# Starting Magnetization Values for QE\nsmag{\n')
    for i in range(1, self.ntypat + 1):
      fout.write(
          "starting_magnetization(" + str(i) + ")=" + str(0) + "\n"
      )
    fout.write('}\n')
