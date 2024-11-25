L = 1.0;  // Square side.

N = 20;
strN = Sprintf("%.0f", N);
h = 1.0 / N;  // Mesh size.

Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};

Line(1)   = {1, 2};
Extrude {0, L, 0} { Line{1}; }
Extrude { {0, 0, 1}, {2 * L, L, 0}, -Pi / 2 } { Line{2}; }
Extrude { {0, 0, 1}, {2 * L, L, 0}, -Pi / 2 } { Line{6}; }
Extrude {0, -L, 0} { Line{10}; }

// Set tags to the boundaries. //////////////////////////////////////////////////

Physical Curve(0) = {1};
Physical Curve(1) = {14};
Physical Curve(2) = {4, 8, 12, 16};
Physical Curve(3) = {3, 7, 11, 15};

Physical Surface(10) = {5, 9, 13, 17};

// Generate mesh ////////////////////////////////////////////////////////////////

Mesh 2;
Save StrCat("../mesh/mesh-u-", strN, ".msh");