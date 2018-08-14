These are part of the codes that provide an automated workflow to generate input files for X-ray spectra calculation,<br/>
run calculation using OCEAN code (http://monalisa.phys.washington.edu/OCEAN/index.html), and perform post-processing on the output spectra. <br/>

The input files for OCEAN can either be genrated using Materials Projects database (https://github.com/materialsproject), or VASP output files. <br/>

To generate input from Materials Projects, an API key is required.<br/> 
   To use, just specify the formula of compounds, or Materials Project ID. This is useful for performing high-throughput calculations.
   e.g.:  Gen_OCEAN_from_MaterProj.py FeO Li5FeO4 Mn7FeCl3O10 <br/>
          please specify the centered atom species:<br/>
          Fe<br/>
          spacegroup of FeO is: 139  I4/mmm<br/>
          spacegroup of Materials Project computed cell is: 139  I4/mmm<br/><br/>

          spacegroup of Li5FeO4 is:  61  Pbca<br/>
          spacegroup of Materials Project computed cell is:  61  Pbca<br/><br/>   

          spacegroup of Mn7FeCl3O10 is: 225  Fm-3m<br/>  
          spacegroup of Materials Project computed cell is: 225  Fm-3m<br/><br/>
   The code will automatically figure out the intrinsic core-hole lifetime broadening and Hubbard+U parameters if there are transition metal elements.<br/>
   K-point mesh for Brouillon zone integration and smearing parameters (metal or insulator) will be determined based on Materials Project entries. <br/>
   The space group and multiplicity of the absorbing species will also be checked, to avoid redundant BSE calculations on atoms with identical Wyckoff positions.<br/>

To generate input from VASP strucutres (POSCAR or CONTCAR). <br/>
  compile Gen_OCEAN_from_POSCAR.cxx: g++ -std=c++11 -o Gen_OCEAN_from_POSCAR.exe Gen_OCEAN_from_POSCAR.cxx <br/>
  execute in the VASP output directory. CONTCAR file will be used if it exists.<br/>

To obtain the averge of perhaps hundreds of simulated spectra, use Average_OCEAN_spectra.sh and Calc_Ocean_Ave.cxx.<br/>
   To use, Average_OCEAN_spectra.sh  $folder $absorbing_atom $starting_index $ending_index $n,l quantum numbers<br/>
   e.g.:   Average_OCEAN_spectra.sh ./Fe2O3 Fe 1 6 1s<br/>
