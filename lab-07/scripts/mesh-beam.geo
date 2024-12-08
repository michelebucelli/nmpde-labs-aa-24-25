L = 1.0;  // Square side.

N = 20;
strN = Sprintf("%.0f", N);
h = 1.0 / N;  // Mesh size.

Point(1) = {0, 0, 0, h};
Point(2) = {0, L, 0, h};

Line(1)   = {2, 1};
Extrude {0, 0, L} { Line{1}; }
Extrude {10, 0, 0} { Surface{5}; }

// Set tags to the boundaries. //////////////////////////////////////////////////

Physical Surface(0) = {5};
Physical Surface(1) = {27};
Physical Surface(2) = {18};
Physical Surface(3) = {26};
Physical Surface(4) = {14};
Physical Surface(5) = {22};

Physical Volume(10) = {1};

// Generate mesh ////////////////////////////////////////////////////////////////

Mesh 3;
Save StrCat("../mesh/mesh-beam-", strN, ".msh");