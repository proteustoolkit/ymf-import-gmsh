h = 1.0;
Point(1) = {0,0,0,h}; Point(2) = {1,0,0,h}; Point(3) = {1,1,0,h}; Point(4) = {0,1,0,h};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Line Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};
out[] = Extrude {0,0,1} { Surface{1}; };
Physical Surface(1) = {1}; Physical Surface(2) = {out[0]};
Physical Surface(3) = {out[2]}; Physical Surface(4) = {out[3]};
Physical Surface(5) = {out[4]}; Physical Surface(6) = {out[5]};
Physical Volume(7) = {out[1]};
