from pymatgen.ext.matproj import MPRester as MP
from pymatgen.core.periodic_table import Element


def write_to_file(compound, fout, ctr_atom):
  """ Write the parametes into input file for X-ray spectra calculation.

      Arguments:
      compound : MPEntry class
        Defined in MPEntry.py file. A class with multiple methods to
        extract properties from database and to compute input parameters.

      fout : str
        Output file name.

      ctr_atom : str
        Photon absorbing atom species.

      znucl : list [int]
        Atomic numbers (Z) of atoms in the cell.

      pp_list : list[str]
        File names of pseudo-potential files.

      natom : int
        Number of total atoms in the cell.

      ntypat  : int
        Number of different types of atoms.

      typat  : list [int]
        List of atom indices by the rank it has in znucl.

      band_gap  : float
        DFT  band gap from materails project.

      occopt  : int
        Occupation used in DFT calculations. Determined by band gap value.

      fband  : float
        Ab-init parameter. Not used in Quantum Espresso.

      rprim  : array (float) (3 x 3)
        Primative vectors defining cell shape.

      xred : array (float) (natom x 3)
        Fractional coordinates of all atoms in the cell.

      ngkpt  : list[int]
        K-point mesh used in ground-state DFT calculation.
  """
  znucl = compound.get_znucl()
  pp_list = compound.get_pp_list()
  natom = compound.structure.num_sites
  ntypat = compound.structure.ntypesp
  typat = compound.get_typat()
  band_gap = MP().query(compound.ID, ['band_gap'])[0]['band_gap']
  occopt = 1 if band_gap > 0 else 3
  fband = 0.4 if occopt == 1 else 2
  rprim = compound.structure.lattice.matrix
  xred = compound.structure.frac_coords
  ngkpt = compound.get_ngkpt()

  fout.write("# number of element types\n"
             "ntypat  " + str(ntypat) + "\n"
             "znucl  { " + " ".join(str(e) for e in znucl) + ' }\n'
             "# number of total atoms\n"
             "natom  " + str(natom) + "\n"
             "pp_list  { " + "  ".join(pp_list) + ' }\n'
             "# occupation, 1 or 3(smearing)\n"
             "occopt  " + str(occopt) + "\n"
             "fband  " + str(fband) + "\n"
             "# opf control files\n"
             "opf.fill  { " + str(Element(ctr_atom).Z) + " " + ctr_atom + ".fill }\n"
             "opf.opts  { " + str(Element(ctr_atom).Z) + " " + ctr_atom + ".opts }\n"
             "acell  { 1.889725989  1.889725989  1.889725989 }\n"
             "rprim  {\n")
  for line in rprim:
    fout.write("    " + " ".join(str("%16.12f" % e) for e in line) + '\n')
  fout.write('}\n'
             "typat  { " + " ".join(str(e) for e in typat) + ' }\n'
             "xred  {\n")
  for line in xred:
    fout.write("    " + " ".join(str("%16.12f" % e) for e in line) + '\n')
  fout.write('}\n'
             "# Kpt mesh for ground state density calculation\n"
             "ngkpt { " + " ".join(str(e) for e in ngkpt) + ' }\n'
             "# Kpt mesh for final states\n"
             "nkpt  { }\n "
             "# Kpt mesh for screening calculation\n"
             "screen.nkpt  { 2 2 2 }\n"
             "# Total bands for final states\n"
             "nbands\n"
             "# Total bands for screening calculation\n"
             "screen.nbands\n"
             "# edge information # number of edges to calculate # atom number, n, l\n"
             "nedges  " + str(len(compound.edge_info.uniq_ctr_atom_pos)) + '\n'
             "edges {\n")
  for i, j in zip(compound.edge_info.uniq_ctr_atom_pos, compound.edge_info.uniq_ctr_atom_multiplicity):
    fout.write("     {:2d}   1   0   # multiplicity: {:2d}\n".format(i, j))
  fout.write('}\n# Starting Magnetization Values for QE\nsmag{\n')
  for i in range(1, ntypat + 1):
    fout.write(
        "starting_magnetization(" + str(i) + ")=" + str(0) + "\n"
    )
  fout.write('}\n')
