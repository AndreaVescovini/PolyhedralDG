## Meshes
- `cube_strXt.mesh` contains a structured tetrahedral mesh made of X elements.
- `cube_strXh.mesh` contains a structured hexahedral mesh of X/6 elements.
- `cube_strXht.mesh` contains a tetrahedral/hexahedral mesh of 7X/12 elements.
- `cube_strXp.mesh` contains a polyhedral mesh created from `cube_strXt.mesh`.
- `structured_cube_gen.edp` is a FreeFem++ source file used to generate the
  tetrahedral meshes, it requires some functions contained in `cube.idp`.
- `tethexmesh_gen.m` is a Matlab routine that generates a mixed tetrahedral/hexahedral
  mesh starting from a hexahedral one.
- `metis_script.sh` is a bash script that can be used to launch the software METIS
  for the generation of the polyhedral meshes.

  METIS (http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) is used in order
  to produce a polyhedral mesh starting from a tetrahedral one. In can work
  partitioning the dual graph of the mesh (i.e. each element becomes a node of the
  graph) with a mutilevel k-way algorithm.
