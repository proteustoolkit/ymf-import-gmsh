"""gmsh meshes to ymf meshes: checked against the files, against each other across MSH formats,
and, where proteus is installed, against proteus's in-tree path (``msh2simplex`` and the
triangle/tetgen readers): the acceptance test of C7.1 round 4."""

import os
import shutil
import subprocess

import numpy as np
import pytest

from ymf import archive
from ymf.mesh import check_facets, edge_midpoint_offsets

from ymf_import_gmsh import GMSH_EDGES, SPACE, __version__, convert, read_msh
from ymf_import_gmsh.cli import main

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
FORMATS = ("msh2", "msh40", "msh41")


def _copy(tmp_path, name):
    """Copy ``name.msh`` (a fixture) into tmp_path; return its file base."""
    shutil.copy(os.path.join(FIXTURES, name + ".msh"), tmp_path)
    return str(tmp_path / name)


def _read(prefix):
    """The space of a converted mesh, validated, its arrays by name, and its facet checks."""
    domain, _ = archive.read_ymf(prefix + ".ymf")
    space = archive.validate_approximation(domain)["times"][0]["spaces"][SPACE]
    base = os.path.dirname(prefix)
    arrays = {"nodes": archive.load_array(space["geometry"]["data_item"], base),
              "elements": archive.load_array(space["topology"]["data_item"], base)}
    for a in space["attributes"]:
        arrays[a["name"]] = archive.load_array(a["data_item"], base)
    if "facets" in space:
        arrays["facets"] = archive.load_array(space["facets"]["topology"]["data_item"], base)
        arrays["facetTags"] = archive.load_array(space["facets"]["tags"]["data_item"], base)
    return space, arrays


def _measure(arrays, dim):
    """Area (2D) or volume (3D) of the cells' vertices (the first dim+1 nodes of each)."""
    x = arrays["nodes"][arrays["elements"][:, :dim + 1]]
    a = x[:, 1:] - x[:, :1]
    return np.abs(np.linalg.det(a[:, :, :dim])).sum() / (2 if dim == 2 else 6)


# --------------------------------------------------------------------------
# the converter alone
# --------------------------------------------------------------------------


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("order", (1, 2))
def test_two_regions_square(tmp_path, fmt, order):
    """Cell tags 1/2 per region, side tags 1-4, a physical point; all three MSH formats."""
    prefix = convert(_copy(tmp_path, "square2_o%d.%s" % (order, fmt)))
    space, arrays = _read(prefix)
    assert space["topology"]["type"] == ("Triangle" if order == 1 else "Tri_6")
    assert space["facets"]["topology"]["type"] == ("Polyline" if order == 1 else "Edge_3")
    assert arrays["elements"].shape[1] == (3 if order == 1 else 6)
    assert _measure(arrays, 2) == pytest.approx(1.0)
    x = arrays["nodes"][arrays["elements"][:, :3]].mean(axis=1)
    assert (arrays["elementMaterialTypes"] == np.where(x[:, 0] < 0.5, 1, 2)).all()
    assert set(check_facets(space, str(tmp_path))) == {1, 2, 3, 4}
    # the side a facet is on: bottom y=0 -> 1 (two lines of tag 1 in the .geo), right x=1 -> 2, top -> 3, left -> 4
    f = arrays["nodes"][arrays["facets"][:, :2]].mean(axis=1)
    side = np.where(f[:, 1] == 0, 1, np.where(f[:, 0] == 1, 2, np.where(f[:, 1] == 1, 3, 4)))
    assert (arrays["facetTags"] == side).all()
    assert arrays["nodeMaterialTypes"][0] == 9                       # Physical Point(9): the origin
    assert (arrays["nodes"][0] == 0).all()
    assert arrays["nodes"][:, 2].any() == False                      # noqa: E712


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("order", (1, 2))
def test_cube(tmp_path, fmt, order):
    prefix = convert(_copy(tmp_path, "cube_o%d.%s" % (order, fmt)))
    space, arrays = _read(prefix)
    assert space["topology"]["type"] == ("Tetrahedron" if order == 1 else "Tet_10")
    assert space["facets"]["topology"]["type"] == ("Triangle" if order == 1 else "Tri_6")
    assert _measure(arrays, 3) == pytest.approx(1.0)
    assert set(arrays["elementMaterialTypes"].tolist()) == {7}
    assert set(check_facets(space, str(tmp_path))) == {1, 2, 3, 4, 5, 6}
    # face tags: 1 z=0, 2 z=1, 3 y=0, 4 x=1, 5 y=1, 6 x=0 (the extrusion's side order)
    f = arrays["nodes"][arrays["facets"][:, :3]].mean(axis=1)
    on = {1: f[:, 2] == 0, 2: f[:, 2] == 1}
    for tag in (1, 2):
        assert (arrays["facetTags"][on[tag]] == tag).all()
    walls = {tag: arrays["facets"][arrays["facetTags"] == tag] for tag in range(1, 7)}
    assert all(len(w) for w in walls.values())


@pytest.mark.parametrize("case, n_cell", [("square2", 3), ("cube", 4)])
@pytest.mark.parametrize("order", (1, 2))
def test_the_formats_give_one_mesh(tmp_path, case, n_cell, order):
    """MSH 2.2, 4.0 and 4.1 of one mesh: the same nodes, and the same cells and tags
    (as sets, since the formats order elements differently)."""
    results = {}
    for fmt in FORMATS:
        prefix = convert(_copy(tmp_path, "%s_o%d.%s" % (case, order, fmt)))
        results[fmt] = _read(prefix)[1]
    ref = results["msh2"]
    for fmt in FORMATS[1:]:
        got = results[fmt]
        np.testing.assert_array_equal(got["nodes"], ref["nodes"], err_msg=fmt)
        np.testing.assert_array_equal(got["nodeMaterialTypes"], ref["nodeMaterialTypes"], err_msg=fmt)
        for cells, tags in (("elements", "elementMaterialTypes"), ("facets", "facetTags")):
            a = sorted(zip(map(tuple, got[cells].tolist()), got[tags].tolist()))
            b = sorted(zip(map(tuple, ref[cells].tolist()), ref[tags].tolist()))
            assert a == b, (fmt, cells)


@pytest.mark.parametrize("name, topology", [("square2_o2.msh41", "Tri_6"), ("cube_o2.msh2", "Tet_10"),
                                            ("cube_o2.msh40", "Tet_10")])
def test_second_order_nodes_are_in_ymf_edge_order(tmp_path, name, topology):
    """Straight edges: each edge node is the midpoint of the edge ymf's order says. gmsh's Tet_10
    puts edge (2,3) before (1,3); without the reorder this offset is of order one."""
    prefix = convert(_copy(tmp_path, name))
    arrays = _read(prefix)[1]
    assert edge_midpoint_offsets(topology, arrays["nodes"], arrays["elements"]).max() < 1e-10
    raw = np.asarray(read_msh(os.path.join(str(tmp_path), name)).cells)
    if topology == "Tet_10":
        assert edge_midpoint_offsets(topology, arrays["nodes"], raw).max() > 0.1
        assert GMSH_EDGES["Tet_10"][4] == (2, 3)


def test_the_stored_proteus_fixtures_convert(tmp_path):
    """gmsh_square and gmsh_cube, as proteus's importer tests store them (4 point tags, side tags)."""
    s_space, s = _read(convert(_copy(tmp_path, "gmsh_square")))
    c_space, c = _read(convert(_copy(tmp_path, "gmsh_cube")))
    assert s["nodes"].shape == (98, 3) and s["elements"].shape == (162, 3) and len(s["facets"]) == 32
    assert c["nodes"].shape == (341, 3) and c["elements"].shape == (1125, 4) and len(c["facets"]) == 540
    assert set(s["elementMaterialTypes"].tolist()) == {1}
    assert set(check_facets(s_space, str(tmp_path))) == {1, 2, 3, 4}
    assert set(check_facets(c_space, str(tmp_path))) == {1, 2, 3, 4, 5, 6}
    assert s["nodeMaterialTypes"][:4].tolist() == [1, 2, 3, 4]       # the physical points


def test_provenance_names_the_file_it_read(tmp_path):
    import hashlib
    prefix = convert(_copy(tmp_path, "gmsh_square"))
    space, _ = _read(prefix)
    assert space["provenance"] == {
        "tool": "ymf-import-gmsh", "version": __version__,
        "files": [{"name": "gmsh_square.msh",
                   "sha256": hashlib.sha256(open(prefix + ".msh", "rb").read()).hexdigest()}]}


MSH_SEGMENTS_INTO_A_NODE = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
$EndNodes
$Elements
6
1 1 2 7 1 1 2
2 1 2 7 1 3 2
3 1 2 8 2 3 4
4 1 2 8 2 1 4
5 2 2 1 1 1 2 3
6 2 2 1 1 1 3 4
$EndElements
"""


def test_a_node_that_only_ends_segments_takes_their_tag(tmp_path):
    """Nodes 2 and 4 only ever end a segment (the second node must be read too)."""
    (tmp_path / "into.msh").write_text(MSH_SEGMENTS_INTO_A_NODE)
    _, arrays = _read(convert(str(tmp_path / "into")))
    assert arrays["nodeMaterialTypes"].tolist() == [7, 7, 7, 8]
    assert arrays["facetTags"].tolist() == [7, 7, 8, 8]


def test_unsupported_input_is_refused_by_name(tmp_path):
    quad = MSH_SEGMENTS_INTO_A_NODE.replace("6\n1 1 2 7", "5\n1 1 2 7").replace(
        "5 2 2 1 1 1 2 3\n6 2 2 1 1 1 3 4", "5 3 2 1 1 1 2 3 4")
    (tmp_path / "quad.msh").write_text(quad)
    with pytest.raises(NotImplementedError, match="element type 3"):
        convert(str(tmp_path / "quad"))
    binary = MSH_SEGMENTS_INTO_A_NODE.replace("2.2 0 8", "2.2 1 8")
    (tmp_path / "bin.msh").write_text(binary)
    with pytest.raises(NotImplementedError, match="binary"):
        convert(str(tmp_path / "bin"))
    (tmp_path / "old.msh").write_text(MSH_SEGMENTS_INTO_A_NODE.replace("2.2 0 8", "1.0 0 8"))
    with pytest.raises(NotImplementedError, match="version 1.0"):
        convert(str(tmp_path / "old"))


def test_the_command_line(tmp_path, capsys):
    base = _copy(tmp_path, "gmsh_square")
    assert main([base + ".msh", "-o", str(tmp_path / "out")]) == 0
    assert capsys.readouterr().out.strip() == str(tmp_path / "out") + ".ymf"
    assert os.path.exists(str(tmp_path / "out.h5"))
    assert main([str(tmp_path / "missing")]) == 1


# --------------------------------------------------------------------------
# acceptance: the same mesh as proteus's in-tree path (skipped without proteus)
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def proteus_meshtools():
    MeshTools = pytest.importorskip("proteus.MeshTools")
    from proteus import Comm
    Comm.init()
    return MeshTools


def _in_tree(MeshTools, filebase, nd):
    """The in-tree path of a use_gmsh domain: msh2simplex, tetgen (3D), the triangle/tetgen reader."""
    MeshTools.msh2simplex(filebase, nd=nd)
    if nd == 2:
        mesh = MeshTools.TriangularMesh()
        mesh.generateFromTriangleFiles(filebase, 1)
        return mesh
    tetgen = shutil.which("tetgen")
    if tetgen is None:
        pytest.skip("no tetgen")
    subprocess.run([tetgen, "-Vfeen", filebase + ".ele"], check=True, capture_output=True)
    for ext in ("ele", "node", "face"):
        os.replace("%s.1.%s" % (filebase, ext), "%s.%s" % (filebase, ext))
    mesh = MeshTools.TetrahedralMesh()
    mesh.generateFromTetgenFiles(filebase, 1)
    return mesh


ACCEPTANCE = [("gmsh_square", 2), ("gmsh_cube", 3), ("square2_o1.msh2", 2), ("cube_o1.msh2", 3)]
MESH_ARRAYS = ("nodeArray", "elementNodesArray", "elementMaterialTypes", "nodeMaterialTypes",
               "elementBoundaryNodesArray", "exteriorElementBoundariesArray",
               "elementBoundaryMaterialTypes")


@pytest.mark.parametrize("name, nd", ACCEPTANCE)
def test_proteus_reads_the_mesh_msh2simplex_reads(tmp_path, monkeypatch, proteus_meshtools, name, nd):
    """The round 4 acceptance: through the ymf mesh, the mesh msh2simplex and the triangle/tetgen
    reader build (nodes, cells, material and boundary tags), bit-equal."""
    monkeypatch.chdir(tmp_path)
    filebase = _copy(tmp_path, name)
    prefix = convert(filebase, str(tmp_path / "converted"))
    via_ymf = (proteus_meshtools.TriangularMesh if nd == 2 else proteus_meshtools.TetrahedralMesh)()
    via_ymf.generateFromYmfFile(prefix)
    in_tree = _in_tree(proteus_meshtools, filebase, nd)
    for key in MESH_ARRAYS:
        np.testing.assert_array_equal(getattr(via_ymf, key), getattr(in_tree, key), err_msg=key)
