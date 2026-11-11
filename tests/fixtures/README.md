# Fixtures

| file | what | source and licence |
|---|---|---|
| `gmsh_square.{geo,msh}`, `gmsh_cube.{geo,msh}` | the unit square (98 nodes, 162 triangles, 4 side tags, 4 point tags) and cube (341 nodes, 1125 tetrahedra, 6 face tags), gmsh 4.12.2, MSH 2.2 | proteustoolkit/proteus `test/mesh_tests/importer_meshes/` (C7.1 importer tests), MIT (`LICENSE.proteus`) |
| `square2.geo`, `cube.geo` | the two small cases below: the unit square split into two regions, the unit cube | written for this repository |
| `square2_o{1,2}.{msh2,msh40,msh41}.msh`, `cube_o{1,2}.{msh2,msh40,msh41}.msh` | those cases at element order 1 and 2 in MSH 2.2, 4.0 and 4.1 ASCII | made by `generate.py` with gmsh 4.15.2 |

Regenerate the second group with `python generate.py` (the `gmsh` executable on `PATH`).

`gmsh_{square,cube}_reference.npz`: seven mesh arrays (`nodeArray`, `elementNodesArray`, `elementMaterialTypes`, `nodeMaterialTypes`, `elementBoundaryNodesArray`, `exteriorElementBoundariesArray`, `elementBoundaryMaterialTypes`) built from `gmsh_square.msh` / `gmsh_cube.msh` by proteus `bf846a47`'s `msh2simplex`, `tetgen -Vfeen` (3D) and the triangle/tetgen readers. Pinned so the check needs no proteus.
