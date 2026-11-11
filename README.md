# ymf-import-gmsh

Convert gmsh `.msh` meshes (MSH 2.2 and 4.x) to ymf meshes.

A `ymf-import-<code>` converter of proteus's ymf mesh campaign (C7.1): ymf is the only mesh representation
proteus reads, and each external format gets its own converter repository. This repository depends on ymf;
ymf and proteus never depend on it. It reads the file format and does not link gmsh.
