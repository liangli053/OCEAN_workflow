from pymatgen.ext.matproj import MPRester as MP
from pymatgen.core.periodic_table import Element
import numpy as np
from math import pi
from Find_unique_sites import Find_unique_sites


class MPEntry:
  """ Use MPRester to interface with Materials Project database and extract all
      related material properties.
  """

  def __init__(self, MP_ID, ctr_atom):
    """
        Arguments:
        self.ID : str
          Materials Project ID of the compound.

        self.formula : str
          Chemical formula.

        self.ctr_atom : str
          Photon absorbing atom species.

        self.structure : MP class
          A materials project class containing all POSCAR information.

        self.elements : list [str]
          Element types in the simulation cell.

        self.edge_info : Find_unique_sites class
          Contains the coordinates of crystallographically distinct ctr_atom
          and corresponding multiplicity.
    """

    self.ID = MP_ID
    self.formula = MP().get_data(self.ID, "vasp", "pretty_formula")[0]['pretty_formula']
    self.ctr_atom = ctr_atom
    self.structure = MP().get_structures(self.ID)[0]
    self.elements = self.structure.symbol_set
    self.edge_info = Find_unique_sites(self.ID, self.ctr_atom)
    self.print_spacegroup_info()

  def print_spacegroup_info(self):
    """ Calculate the space group number using SpaceGroupAnalyzer and and compare with
        the data obtained by MP().query. Issue a warning if not consistent.

        Arguments:
        SP_symbol : str
          Space group symbol listed on Materials Project webpage.

        SP_number : str
          Space group number listed on Materials Project webpage.

        computed_SP_symbol : str
          Space group symbol obtained by SpaceGroupAnalyzer.

        computed_SP_number : str
          Space group number obtained by SpaceGroupAnalyzer.
    """
    SP_symbol = MP().query(self.ID, ['spacegroup'])[0]['spacegroup']['symbol']
    SP_number = MP().query(self.ID, ['spacegroup'])[0]['spacegroup']['number']
    computed_SP_number = self.edge_info.symm_finder.get_space_group_number()
    computed_SP_symbol = self.edge_info.symm_finder.get_space_group_symbol()
    print("spacegroup of {} is: {:3d}  {:7s}\n".format(self.formula, SP_number, SP_symbol))
    print("spacegroup of Materials Project computed cell is: {:3d}  {:7s}\n".format(computed_SP_number, computed_SP_symbol))
    if SP_number != computed_SP_number:
      print("WARNING: The spacegroup of {} found in Materials Project POSCAR may be inaccurate,\n \
       this can cause redundant calculation, but does not affet results\n".format(self.formula))

  def get_znucl(self):
    """ Return the znucl parameter - atomic numbers of atoms in the cell.

        Returns:
        znucl : list [int]
          Atomic numbers (Z) of atoms in the cell.
    """

    znucl = []
    for i in self.elements:
      tmp = Element(i)
      znucl.append(tmp.Z)
    return znucl

  def get_typat(self):
    """ Return typat parameter - list of atoms by the rank it has in znucl.
        e.g. in Cu4O2 - {1 1 1 1 2 2}

        Returns:
        typat : list [str]
    """

    num_of_element = []
    typat = []
    # non-reduced chemical formula according to total atom numbers
    total_formula = MP().query(self.ID, ['unit_cell_formula'])[0]['unit_cell_formula']
    for i in range(len(self.elements)):
      # number of atoms of specific species
      num_of_element.append(int(total_formula[self.elements[i]]))
      typat.extend([i + 1] * num_of_element[i])
    return typat

  def get_pp_list(self):
    """ Return pp_list parameter - pseudopotential file names.

        Returns:
          pp_list : list[str]
            File names of pseudo potential files.
    """

    pp_list = []
    for i in self.elements:
      tmp = Element(i)
      pp_list.append(str('{:02d}'.format(tmp.Z)) + '-' + i.lower() + ".lda.fhi")
    return pp_list

  def get_ngkpt(self):
    """ Return K-points for ground-state charge density calculation.
        Use a KSPACING parameter to define the lower bound of K-point spacing.
        Denser mesh if band_gap < 0 (metallic system).

        Arguments:
        band_gap : float
          Band gap value from material project.

        a, b, c: float
          Lattice parameters along a, b, c directions.

        Returns:
        ngkpt : list [int]
    """
    ngkpt = []
    a = float(self.structure.lattice.a)
    b = float(self.structure.lattice.b)
    c = float(self.structure.lattice.c)
    band_gap = MP().query(self.ID, ['band_gap'])[0]['band_gap']
    if band_gap > 0:
      KSPACING = 0.3
    else:
      KSPACING = 0.15
    ngkpt.append(int(2 * pi / a / KSPACING))
    ngkpt.append(int(2 * pi / b / KSPACING))
    ngkpt.append(int(2 * pi / c / KSPACING))
    return ngkpt
