a = 1;
b = 1;
c = 1;
d = 5;
e = 2;

N = 3;
strN = Sprintf("%.0f", N);
h = 1.0 / N;  // Mesh size.

Point(1) = {0, 0, c, h};
Point(2) = {0, 0, a + c, h};
Line(1) = {2, 1};

Extrude {b, 0, 0} { Line{1}; }

Point(5) = {b, 0, 0, h};
Line(5) = {4, 5};

Extrude {d, 0, 0} { Line{2, 5}; }

Extrude {0, e, 0} { Surface{5, 9, 13}; }

// Set tags to the boundaries. //////////////////////////////////////////////////

Physical Surface(0) = {22};
Physical Surface(1) = {5, 9, 13, 26, 34, 35, 56, 57, 66, 70, 79};
Physical Surface(2) = {52, 74};

Physical Volume(10) = {1, 2, 3};

// Generate mesh ////////////////////////////////////////////////////////////////

Mesh 3;
Save StrCat("../mesh/mesh-step-", strN, ".msh");
