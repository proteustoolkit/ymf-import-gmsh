"""Regenerate the small gmsh fixtures: ``python generate.py`` (needs the ``gmsh`` executable on PATH).

Writes, beside this script, ``<case>_o<order>.<format>.msh`` for each case, element order (1, 2)
and MSH format (msh2 = 2.2, msh40, msh41 = 4.0 and 4.1 ASCII). The cases are in the code below:
the unit square split into two regions (a cell tag per region, a tag per side), and the unit
cube (a tag per face). ``gmsh_square`` and ``gmsh_cube`` (copied from proteus's importer tests,
gmsh 4.12.2, msh2) are not made here.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

CASES = {
    "square2": ("-2", """
h = 0.5;
Point(1) = {0,0,0,h}; Point(2) = {0.5,0,0,h}; Point(3) = {1,0,0,h};
Point(4) = {1,1,0,h}; Point(5) = {0.5,1,0,h}; Point(6) = {0,1,0,h};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,5}; Line(5) = {5,6};
Line(6) = {6,1}; Line(7) = {2,5};
Line Loop(1) = {1,7,5,6}; Plane Surface(1) = {1};
Line Loop(2) = {2,3,4,-7}; Plane Surface(2) = {2};
Physical Line(1) = {1,2}; Physical Line(2) = {3}; Physical Line(3) = {4,5}; Physical Line(4) = {6};
Physical Surface(1) = {1}; Physical Surface(2) = {2};
Physical Point(9) = {1};
"""),
    "cube": ("-3", """
h = 1.0;
Point(1) = {0,0,0,h}; Point(2) = {1,0,0,h}; Point(3) = {1,1,0,h}; Point(4) = {0,1,0,h};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Line Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};
out[] = Extrude {0,0,1} { Surface{1}; };
Physical Surface(1) = {1}; Physical Surface(2) = {out[0]};
Physical Surface(3) = {out[2]}; Physical Surface(4) = {out[3]};
Physical Surface(5) = {out[4]}; Physical Surface(6) = {out[5]};
Physical Volume(7) = {out[1]};
"""),
}


def main():
    for case, (dim, geo) in CASES.items():
        with open(os.path.join(HERE, case + ".geo"), "w") as f:
            f.write(geo.lstrip())
        for order in (1, 2):
            for fmt in ("msh2", "msh40", "msh41"):
                out = os.path.join(HERE, "%s_o%d.%s.msh" % (case, order, fmt))
                subprocess.run(["gmsh", os.path.join(HERE, case + ".geo"), dim, "-order", str(order),
                                "-format", fmt, "-nopopup", "-v", "0", "-o", out], check=True)
                print(out, os.path.getsize(out))


if __name__ == "__main__":
    sys.exit(main())
