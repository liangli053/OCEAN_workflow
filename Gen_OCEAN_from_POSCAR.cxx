#include <algorithm>
#include <iomanip>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "GotoLine.h"
using namespace std;

int main() {
  bool POSCAR_exists = false;
  bool CONTCAR_exists = false;
  string input_file;
  string line;
  string word;
  int number, tot_atm_num;
  const int line_lattice_st = 3;
  const int line_atm_species = 6;
  const int line_atm_nums = 7;
  const int line_coords_st = 10;
  vector<string> lattice;
  vector<string> atm_species;
  vector<int> atm_nums;
  vector<int> typat;
  vector<vector<double>> atm_coords;
// read from POSCAR, if CONTCAR exists, read CONTCAR
  if (ifstream("POSCAR")) { POSCAR_exists = 1; input_file = "POSCAR";}
  if (ifstream("CONTCAR")) { CONTCAR_exists = 1; input_file = "CONTCAR";}
  if (POSCAR_exists == false && CONTCAR_exists == false) {cout << "no POSCAR or CONTCAR are found" << endl;} 
  fstream read_from(input_file.c_str());
// read lattice parameters
  GotoLine(read_from, line_lattice_st);
  for (int i=1; i<=3; ++i) {
    getline(read_from, line);
    lattice.push_back(line);
  }; 
// read atomic species form Line 6
  GotoLine(read_from, line_atm_species);
  getline(read_from, line);
  istringstream ss(line);
  while (ss >> word) { atm_species.push_back(word); }
// read atom numbers fro line 7  
  GotoLine(read_from, line_atm_nums);
  getline(read_from, line);
  istringstream ss1(line);
  while (ss1 >> number) { atm_nums.push_back(number); }
  tot_atm_num = accumulate(atm_nums.begin(), atm_nums.end(), 0);
// read atom coordinates
  GotoLine(read_from, line_coords_st);
  for (int i=1; i<= tot_atm_num; ++i) {
    getline(read_from, line);
    istringstream ss2(line);
    vector<double> coords_row;
    double x;
    for (int j=1; j<=3; ++j) {
      ss2 >> x;
      coords_row.push_back(x);
    };
    atm_coords.push_back(coords_row);
  };

  read_from.close();
  for (vector<int>::size_type i=0; i<atm_nums.size(); ++i) {
    for (int j=1; j <= atm_nums[i]; ++j)
      typat.push_back(i+1);
  };

// wirte to "ocean.in"
  cout << "Don't forget to adjust pp_list and edges" << endl;
  ofstream write_to("ocean.in");
  write_to << "control 0" << endl 
           << "para_prefix { mpiexec -n 160 }" << endl 
           << "# which stages of the calculation to run (one or more of): paw, dft, prep, screen, bse, all" << endl 
           << "stages { all }" << endl 
           << "dft { abinit }" << endl 
           << "# Ntypes of atoms" << endl 
           << "ntypat " << atm_species.size() << endl 
           << "# Z num for types" << endl 
           << "znucl { 29 8 }" << endl 
           << "# pseudo location" << endl 
           << "ppdir { '/lcrc/project/PhotoCu2O/OCEAN/bulk' }" << endl 
           << "# pseudopotentials"  << endl
           << "pp_list{  29-Cu.LDA.fhi" << endl 
           << "          08-O.LDA.fhi }" << endl 
           << "# Kinetic Energy cutoff (in Ry for QE)" << endl
           << "ecut 70" << endl
           << "# SCF Energy tolerance" << endl
           << "toldfe 1.1d-8" << endl
           << "# SCF wftol" << endl
           << "tolwfr 1.1d-15" << endl
           << "# SCF iterations" << endl
           << "nstep 250" << endl
           << "# Static dielectric const (taken from experiments)" << endl
           << "diemac 6.46" << endl
           << "# Kpt mesh for ground state density calculation" << endl
           << "ngkpt { 3 2 1  }" << endl
           << "# Kpt mesh for final states" << endl
           << "nkpt { 3 2 1 }" << endl
           << "# occupation, 1 or 3(smearing)" << endl
           << "occopt 3" << endl
           << "# SCF mixing" << endl
           << "#mixing { 0.2 }" << endl
           << "# Kpt mesh for screening calculation" << endl
           << "paw.nkpt { 3 2 1 }" << endl
           << "# paw control files" << endl
           << "paw.fill{ 29 Cu.fill }" << endl
           << "paw.opts{ 29 Cu.opts }" << endl
           << "# radius for paw reconstruciton, can be multiple values,only one is used as screening radius as in cnbse.ras" << endl
           << "paw.shells{ 4.0 }" << endl
           << "cnbse.rad{ 4.0 }" << endl
           << "# xmesh" << endl
           << "#cnbse.xmesh { 10 10 10 }" << endl
           << "# spectral broadening in eV" << endl
           << "cnbse.broaden 1.55" << endl
           << "# The code will figure out a good plot range" << endl
           << "#cnbse.spect_range{ 1500  -20  50 }" << endl
           << "#Scaling Factor" << endl
           << "scfac 0.8" << endl << endl

           << "# Mag. of latt. vec. in Bohr" << endl
           << "acell { 1.889725989  1.889725989  1.889725989 }" << endl
           << "# Cart comps. of latt. vecs." << endl
           << "rprim {"  << endl;
// write lattice parameters
  for (auto i = lattice.begin(); i!=lattice.end(); ++i) {
    write_to << *i << endl;
  };
  write_to << "}" << endl 
           << "# N atoms in unit cell" << endl 
           << "natom " << tot_atm_num << endl 
           << "# Type of each atom" << endl 
           << "typat { ";
// write typat
  for (auto i = typat.begin(); i!=typat.end(); ++i) {
    write_to << *i << " ";
  };
  write_to <<"}" << endl 
           <<"# Relative positions of atoms" << endl 
           <<"xred {" << endl;
// write atom coordinates
  for (auto i = atm_coords.begin(); i!=atm_coords.end(); ++i) {
    write_to << fixed;
    write_to <<"   ";
    for (auto j = (*i).begin(); j!=(*i).end(); ++j) {
       write_to << setw(18) << setprecision(16) << *j << "   ";
    };
    write_to << endl;
  };
  write_to << "}" << endl << endl 
           << "fband 2" << endl 
           << "# Total bands for final states"  << endl
           << "nbands 900" << endl
           << "# Total bands for screening wave functions" << endl
           << "paw.nbands 1600" << endl 
           << "# edge information # number of edges to calculate # atom number, n quantum number, l quantum number" << endl 
           << "nedges 1" << endl 
           << "edges { 1 1 0 }" << endl;
}
