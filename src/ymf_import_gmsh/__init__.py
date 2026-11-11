"""Convert gmsh ``.msh`` meshes (ASCII MSH 2.2 and 4.0/4.1) to ymf meshes.

The ymf mesh written is the one proteus reads (``proteus.MeshTools.readMeshYmf``): a
format 2 archive whose one instant holds one uniform space, ``Spatial_Domain``, with
``prefix.ymf`` describing it and ``prefix.h5`` holding ``nodes``, ``elements``,
``nodeMaterialTypes``, ``elementMaterialTypes`` and, when the file has tagged boundary
elements, ``facets`` and ``facetTags``. The mesh carries ymf's ``provenance``.

What is read, and how it maps:

* Cells are the highest-dimension elements: triangles (``Triangle``, or ``Tri_6`` for
  second order) or tetrahedra (``Tetrahedron``, ``Tet_10``). Other element types are
  refused by name rather than dropped.
* Numbers are zero-based and nodes are in the order of their gmsh tags, so a file
  numbered 1..N keeps its order. Cell order is the file's.
* A cell's tag (``elementMaterialTypes``) is its first physical group, 0 if it has none.
* Facets are the tagged boundary elements: segments (``Polyline``/``Edge_3``) under
  triangles, triangles (``Triangle``/``Tri_6``) under tetrahedra, tagged by their first
  physical group. An element in no physical group is not a facet (tag 0 is untagged).
* A node's tag (``nodeMaterialTypes``) follows proteus's in-tree reader
  (``msh2simplex``): a physical point's tag, else (3D) the first tagged boundary triangle
  that has the node, else the first tagged segment that has it, else 0.
* Second-order cells are reordered to ymf's edge-node order (VTK's; gmsh's ``Tet_10``
  numbers the last two edge nodes the other way round).
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

__version__ = "0.1.0.dev0"

__all__ = ["GmshMesh", "read_msh", "convert", "SPACE", "__version__"]

#: The space a ymf mesh file holds its mesh in, as in a proteus archive.
SPACE = "Spatial_Domain"

#: gmsh element type -> (ymf topology type, order, nodes, dimension)
ELEMENT_TYPES = {
    15: ("Polyvertex", 1, 1, 0),
    1: ("Polyline", 1, 2, 1),
    8: ("Edge_3", 2, 3, 1),
    2: ("Triangle", 1, 3, 2),
    9: ("Tri_6", 2, 6, 2),
    4: ("Tetrahedron", 1, 4, 3),
    11: ("Tet_10", 2, 10, 3),
}

#: the edges (vertex pairs) of gmsh's second-order cells, in gmsh's edge-node order
GMSH_EDGES = {
    "Tri_6": ((0, 1), (1, 2), (2, 0)),
    "Tet_10": ((0, 1), (1, 2), (2, 0), (0, 3), (2, 3), (1, 3)),
}


@dataclass
class GmshMesh:
    """A gmsh mesh as read: zero-based numbers, nodes in tag order."""

    topology: str                     #: "Triangle", "Tri_6", "Tetrahedron" or "Tet_10"
    nodes: List[List[float]]
    cells: List[List[int]]            #: node numbers per cell, in gmsh's order
    cell_tags: List[int]
    facet_topology: str = ""          #: "Polyline", "Edge_3", "Triangle", "Tri_6"
    facets: List[List[int]] = field(default_factory=list)
    facet_tags: List[int] = field(default_factory=list)
    node_tags: List[int] = field(default_factory=list)
    files: List[str] = field(default_factory=list)


def _tokens(path: str) -> Tuple[Dict[str, List[str]], List[str]]:
    """The file's ``$Section`` ... ``$EndSection`` blocks as line lists, and the format line's words."""
    sections: Dict[str, List[str]] = {}
    name = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("$End"):
                name = None
            elif line.startswith("$"):
                name = line[1:]
                sections[name] = []
            elif name is not None and line:
                sections[name].append(line)
    if "MeshFormat" not in sections:
        raise ValueError("%s: no $MeshFormat section; not a gmsh .msh file" % (path,))
    return sections, sections["MeshFormat"][0].split()


def _physical_of_entities(lines: List[str], legacy: bool) -> Dict[Tuple[int, int], int]:
    """``{(dim, entity tag): first physical tag}`` from a 4.x ``$Entities`` section (0: none)."""
    counts = [int(w) for w in lines[0].split()[:4]]
    result = {}
    pos = 1
    for dim, n in enumerate(counts):
        for _ in range(n):
            words = lines[pos].split()
            pos += 1
            skip = 4 if dim == 0 and not legacy else 7   # tag, then x y z (4.1 point) or a bounding box
            n_phys = int(words[skip])
            result[(dim, int(words[0]))] = int(words[skip + 1]) if n_phys else 0
    return result


def read_msh(filebase: str) -> GmshMesh:
    """Read ``filebase.msh`` (ASCII MSH 2.2, 4.0 or 4.1)."""
    path = filebase + ".msh"
    sections, fmt = _tokens(path)
    version = fmt[0]
    if len(fmt) > 1 and fmt[1] != "0":
        raise NotImplementedError("%s: binary MSH is not converted (write the mesh with "
                                  "Mesh.Binary = 0)" % (path,))
    major = version.split(".")[0]
    if version == "2.2" or major == "2":
        nodes, elements = _read_v2(path, sections)
    elif major == "4":
        nodes, elements = _read_v4(path, sections, version)
    else:
        raise NotImplementedError("%s: MSH version %s (supported: 2.2, 4.0, 4.1)" % (path, version))
    return _assemble(path, nodes, elements)


def _read_v2(path, sections):
    node_lines = sections["Nodes"]
    nodes = {}
    for line in node_lines[1:1 + int(node_lines[0])]:
        w = line.split()
        nodes[int(w[0])] = [float(x) for x in w[1:4]]
    elements = []     # (type, physical tag, node tags)
    element_lines = sections["Elements"]
    for line in element_lines[1:1 + int(element_lines[0])]:
        w = line.split()
        etype, n_tags = int(w[1]), int(w[2])
        n = ELEMENT_TYPES.get(etype, (None, None, 0, None))[2]
        elements.append((etype, int(w[3]) if n_tags else 0, [int(x) for x in w[3 + n_tags:]]))
    return nodes, elements


def _read_v4(path, sections, version):
    legacy = version in ("4", "4.0")      # gmsh writes the 4.0 format as "4 0 8"
    physical = _physical_of_entities(sections["Entities"], legacy)
    node_lines = sections["Nodes"]
    n_blocks = int(node_lines[0].split()[0])
    nodes = {}
    pos = 1
    for _ in range(n_blocks):
        head = [int(x) for x in node_lines[pos].split()]
        pos += 1
        n = head[3]                  # entityDim/entityTag swap places in 4.0; the count is last in both
        if legacy:
            for line in node_lines[pos:pos + n]:
                w = line.split()
                nodes[int(w[0])] = [float(x) for x in w[1:4]]
            pos += n
        else:
            tags = [int(node_lines[pos + i]) for i in range(n)]
            for i, t in enumerate(tags):
                nodes[t] = [float(x) for x in node_lines[pos + n + i].split()[:3]]
            pos += 2 * n
    element_lines = sections["Elements"]
    n_blocks = int(element_lines[0].split()[0])
    elements = []
    pos = 1
    for _ in range(n_blocks):
        head = [int(x) for x in element_lines[pos].split()]
        pos += 1
        dim, entity = (head[1], head[0]) if legacy else (head[0], head[1])
        etype, n = head[2], head[3]
        tag = physical.get((dim, entity), 0)
        for line in element_lines[pos:pos + n]:
            w = line.split()
            elements.append((etype, tag, [int(x) for x in w[1:]]))
        pos += n
    return nodes, elements


def _assemble(path, nodes, elements) -> GmshMesh:
    for etype, _, _ in elements:
        if etype not in ELEMENT_TYPES:
            raise NotImplementedError(
                "%s: gmsh element type %d is not converted (supported: triangles and tetrahedra of "
                "order 1 and 2, with segments and points); a mesh missing elements is worse than "
                "no mesh" % (path, etype))
    by_dim = {}
    for etype, tag, conn in elements:
        by_dim.setdefault(ELEMENT_TYPES[etype][3], []).append((etype, tag, conn))
    cell_dim = max(d for d in by_dim if d >= 2) if any(d >= 2 for d in by_dim) else None
    if cell_dim is None:
        raise ValueError("%s: no triangles or tetrahedra" % (path,))
    kinds = {e[0] for e in by_dim[cell_dim]}
    if len(kinds) != 1:
        raise NotImplementedError("%s: first- and second-order cells in one mesh" % (path,))
    cell_type = kinds.pop()
    topology = ELEMENT_TYPES[cell_type][0]
    order = ELEMENT_TYPES[cell_type][1]
    facet_dim = cell_dim - 1
    order_facets = {e[0] for e in by_dim.get(facet_dim, [])}
    if order_facets - {t for t, v in ELEMENT_TYPES.items() if v[1] == order and v[3] == facet_dim}:
        raise NotImplementedError("%s: boundary elements of another order than the cells" % (path,))

    tags = sorted(nodes)
    index = {t: i for i, t in enumerate(tags)}
    coordinates = [nodes[t] for t in tags]

    def numbers(conn):
        try:
            return [index[t] for t in conn]
        except KeyError as error:
            raise ValueError("%s: element node %s is not a node of the file" % (path, error))

    seen = set()
    cells, cell_tags = [], []
    for etype, tag, conn in by_dim[cell_dim]:
        row = numbers(conn)
        key = frozenset(row)
        if key in seen:                  # gmsh 2.2 writes a cell once per physical group it is in
            continue
        seen.add(key)
        cells.append(row)
        cell_tags.append(tag)
    if len(cells) != len(by_dim[cell_dim]):
        warnings.warn("%s: %d cells are in more than one physical group; the first group's tag is kept"
                      % (path, len(by_dim[cell_dim]) - len(cells)))

    n_vertices = {2: 3, 3: 4}[cell_dim]
    facets, facet_tags, seen = [], [], set()
    for etype, tag, conn in by_dim.get(facet_dim, []):
        if not tag:
            continue
        row = numbers(conn)
        key = frozenset(row[:n_vertices - 1])
        if key in seen:
            continue
        seen.add(key)
        facets.append(row)
        facet_tags.append(tag)
    facet_topology = ""
    if facets:
        facet_topology = ELEMENT_TYPES[next(iter(order_facets))][0]

    node_tags = [0] * len(tags)
    for etype, tag, conn in by_dim.get(0, []):
        node_tags[numbers(conn)[0]] = tag
    # a node takes the first tagged boundary face that has it (3D), then the first tagged segment
    first = [by_dim.get(2, [])] if cell_dim == 3 else []
    for group in first + [by_dim.get(1, [])]:
        for etype, tag, conn in group:
            if tag:
                for n in numbers(conn)[:2 if ELEMENT_TYPES[etype][3] == 1 else 3]:
                    if node_tags[n] == 0:
                        node_tags[n] = tag
    return GmshMesh(topology, coordinates, cells, cell_tags, facet_topology, facets, facet_tags,
                    node_tags, files=[path])


def convert(filebase: str, prefix: str | None = None) -> str:
    """Write the ymf mesh ``prefix.ymf`` + ``prefix.h5`` for ``filebase.msh``; return ``prefix``.

    ``prefix`` is ``filebase`` by default.
    """
    import h5py
    import numpy as np
    from ymf import archive
    from ymf.mesh import file_provenance, reorder_cells

    mesh = read_msh(filebase)
    prefix = filebase if prefix is None else prefix
    cells = np.asarray(mesh.cells, dtype="i")
    facets = np.asarray(mesh.facets, dtype="i")
    if mesh.topology in GMSH_EDGES:
        cells = reorder_cells(cells, mesh.topology, GMSH_EDGES[mesh.topology]).astype("i")
    arrays = [
        ("nodes", np.asarray(mesh.nodes, dtype="d").reshape(-1, 3)),
        ("elements", np.ascontiguousarray(cells)),
        ("nodeMaterialTypes", np.asarray(mesh.node_tags, dtype="i")),
        ("elementMaterialTypes", np.asarray(mesh.cell_tags, dtype="i")),
    ]
    if len(mesh.facets):
        arrays += [("facets", np.ascontiguousarray(facets)),
                   ("facetTags", np.asarray(mesh.facet_tags, dtype="i"))]
    h5_name = os.path.basename(prefix) + ".h5"
    with h5py.File(prefix + ".h5", "w") as f:
        for name, array in arrays:
            f.create_dataset(name, data=array)
    items = {name: archive.data_item_for(array, "%s:/%s" % (h5_name, name)) for name, array in arrays}
    facet_list = None
    if len(mesh.facets):
        facet_list = archive.facets(
            archive.topology(mesh.facet_topology, len(mesh.facets), items["facets"],
                             nodes_per_element=facets.shape[1]),
            items["facetTags"])
    grid = archive.uniform_grid(
        archive.topology(mesh.topology, len(mesh.cells), items["elements"]),
        archive.geometry(items["nodes"]),
        [archive.attribute("nodeMaterialTypes", items["nodeMaterialTypes"], center="Node"),
         archive.attribute("elementMaterialTypes", items["elementMaterialTypes"], center="Cell")],
        facets=facet_list,
        provenance=file_provenance("ymf-import-gmsh", mesh.files, version=__version__))
    domain = archive.new_approximation()
    archive.add_instant(domain, 0.0, {SPACE: grid})
    archive.write_ymf(domain, prefix + ".ymf")
    return prefix
