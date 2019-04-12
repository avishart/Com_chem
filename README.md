# Com_chem
Short scripts for computational chemistry calculations. Mainly designed for the softwares Gaussian09 and Gaussian16, http://gaussian.com/.


-boltzmann.py: Calculate the Boltzman distribution by giving the script a list of output files.  

-bond_scan_file.py: Use the script in a folder of output files (.out or .log) calculated by Gaussian. The bond scan is between two atoms defined in the script. See generate_bond_scan.py.

-Change_basis_func.py: Change the basis set and functionals for a computational chemistry calculation in Gaussian by giving the script an input file (.com). New input files are generated in the same folder.

-Convergence_geo_opt.py: Plot the convergence of geometry optimization calculation in Gaussian. Give the script the output file.

-generate_bond_scan.py: The script generate input files for a compressed and stretched bond. Use a input file (.com) from Gaussian and seperate the structure in two fragments. Then define the bond stretched by giving the two atoms in the script.

-Output_to_input.py: Make a new input file by using the geometry from an output file. Give the script the output file (.out or .log) from Gaussian and the new filename of the input file.

