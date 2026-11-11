"""``ymf-import-gmsh MESH [-o PREFIX]``: convert a gmsh .msh mesh to a ymf mesh."""

import argparse
import sys

from . import __version__, convert


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="ymf-import-gmsh",
        description="Convert MESH.msh (ASCII MSH 2.2 or 4.x) to PREFIX.ymf and PREFIX.h5.")
    parser.add_argument("mesh", help="the mesh file, with or without its .msh suffix")
    parser.add_argument("-o", "--output", help="output prefix (default: the mesh file's)")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args(argv)
    filebase = args.mesh[:-4] if args.mesh.endswith(".msh") else args.mesh
    try:
        prefix = convert(filebase, args.output)
    except (OSError, ValueError, NotImplementedError) as error:
        print("ymf-import-gmsh: %s" % (error,), file=sys.stderr)
        return 1
    print("%s.ymf" % (prefix,))
    return 0


if __name__ == "__main__":
    sys.exit(main())
