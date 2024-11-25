L = 1.0;  // Square side.

N = 40;
strN = Sprintf("%.0f", N);
h = 1.0 / N;  // Mesh size.

Point(1) = {0, 0, 0, h};
Point(2) = {0, L, 0, h};

Line(1)   = {2, 1};
Extrude {1, 0, 0} { Line{1}; }

// Set tags to the boundaries. //////////////////////////////////////////////////

Physical Curve(0) = {1};
Physical Curve(1) = {2};
Physical Curve(2) = {4};
Physical Curve(3) = {3};

Physical Surface(10) = {5};

// Generate mesh ////////////////////////////////////////////////////////////////

Mesh 2;
Save StrCat("../mesh/mesh-square-", strN, ".msh");