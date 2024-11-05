L = 1.0;  // Square side.
R = 0.75; // Circular hole radius.

h = 0.1;  // Mesh size.

// Points. //////////////////////////////////////////////////////////////////////

Point(1) = {0, 0, 0, h};
Point(2) = {0, L, 0, h};
Point(3) = {L, L, 0, h};
Point(4) = {L, R, 0, 0.25 * h};
Point(5) = {L - R, 0, 0, 0.25 * h};
Point(6) = {L, 0, 0, h};

// Lines ////////////////////////////////////////////////////////////////////////

Line(1)   = {1, 2};
Line(2)   = {2, 3};
Line(3)   = {3, 4};
Circle(4) = {4, 6, 5};
Line(5)   = {5, 1};

Line Loop(1) = {1, 2, 3, 4, 5};

// Surfaces /////////////////////////////////////////////////////////////////////

Plane Surface(1) = {1};

// Set tags to the boundaries. //////////////////////////////////////////////////

Physical Line(0) = {1};
Physical Line(1) = {3};
Physical Line(2) = {5};
Physical Line(3) = {2};
Physical Line(4) = {4};
Physical Surface(1) = {1};

// Generate mesh ////////////////////////////////////////////////////////////////

Mesh 2;
Save "mesh.msh";