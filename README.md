# ymf-import-gmsh

Convert gmsh `.msh` meshes (MSH 2.2 and 4.x) to ymf meshes.

A `ymf-import-<code>` converter of proteus's ymf mesh campaign (C7.1): ymf is the only mesh representation
proteus reads, and each external format gets its own converter repository. This repository depends on ymf;
ymf and proteus never depend on it. It reads the file format and does not link gmsh.

```
ymf-import-gmsh MESH.msh [-o PREFIX]      # writes PREFIX.ymf + PREFIX.h5
```

```python
from ymf_import_gmsh import convert
convert("mesh")          # mesh.msh -> mesh.ymf + mesh.h5
```

## What it reads

ASCII MSH 2.2, 4.0 and 4.1: triangles and tetrahedra of order 1 and 2 (`Tri_6`, `Tet_10`), with the
segments and points that carry physical groups. Binary files and other element types (quadrilaterals,
hexahedra, prisms, ...) are refused by name.

| gmsh | ymf mesh |
|---|---|
| nodes, in tag order | `nodes` (always x y z), zero-based |
| highest-dimension elements | `elements`; `Triangle`/`Tri_6` or `Tetrahedron`/`Tet_10` |
| a cell's first physical group | `elementMaterialTypes` (0: none) |
| tagged segments (2D) or triangles (3D) | `facets` + `facetTags` (an element in no physical group is not a facet) |
| physical points, else tagged boundary faces (3D), else tagged segments | `nodeMaterialTypes`, as proteus's removed `msh2simplex` did |
| second-order edge nodes | reordered to ymf's edge order (gmsh's `Tet_10` swaps the last two) |
| file name and SHA-256 | `provenance` |

## Acceptance

`tests/fixtures/gmsh_{square,cube}_reference.npz` are the seven mesh arrays proteus's former in-tree path
(`msh2simplex`, tetgen `-Vfeen` in 3D, the triangle/tetgen readers) built from `gmsh_square.msh` and
`gmsh_cube.msh`, generated at proteustoolkit/proteus `bf846a47`, the parent of the commit that removed
`msh2simplex` (C7.1 round 5). The converted nodes, cells and material tags are compared to them bit-equal in
every run, with no proteus; where proteus is installed, `TriangularMesh`/`TetrahedralMesh.generateFromYmfFile`
is compared on all seven arrays, boundary topology and tags included.

Fixtures and how they were made: `tests/fixtures/README.md`.
