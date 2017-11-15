%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DISCONTINUOUS GALERKIN FINITE ELEMENT METHOD OVER POLYEHEDRAL MESHES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Project of Numerical Analisys for Partial Differential Equations - Prof. Perotto

Author: Andrea Vescovini
Supervisor Prof. Paola Antonietti

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION OF THE CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

In this folder:
- "main3D.m" is used to set the input data and start the solver.
- "print_result.m" is used to perform a convergence analysis and prints the
  results on a file "results*.txt".
- "plotmesh.m" should plot the mesh but is not very effective.
- "hexmesh_gen.m" is a routine that generates a hexahedral mesh starting from
  a tetrahedral one.
- the other files contain functions used by the solver.

In the folder "../meshes":
- "cube_str*t.mesh" contains a structured tetrahedral mesh of * elements.
- "cube_str*h.mesh" contains a structured hexahedral mesh of */6 elements.
- "cube_str*ht.mesh" contains a tetrahedral/hexahedral mesh of 7*/12 elements.
- "cube_str*p.mesh" contains a polyhedral mesh created from "cube_str*t.mesh".
- "metis_script" is a bash script to launch the software METIS for the
  generation of the polyhedral meshes.
- "structured_cube_gen.edp" is a FreeFem++ source file used to generate the
  tetrahedral meshes.
- "tethexmesh_gen.m" is a routine that generates a mixed tetrahedral/hexahedral
  mesh starting from a hexahedral one.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
