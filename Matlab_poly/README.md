## Discontinuous Galerkin finite element method over polyhedral meshes

Project of Numerical Analisys for Partial Differential Equations - Prof. Perotto

**Author :** Andrea Vescovini

**Supervisor :** Prof. Paola Antonietti

### Descrption of the code

- `main3D.m` is used to set the input data and start the solver.
- `print_result.m` is used to perform a convergence analysis and prints the
  results on a file `results*.txt`.
- `export_solution.m` is used by main3D and exports the solution and the mesh
  in a file output.vtk, in order to be visualized with a visualization software
  like Paraview.
- `hexmesh_gen.m` is a routine that generates a hexahedral mesh starting from
  a tetrahedral one.
- the other files contain functions used by the solver.
