L = 1.0; // Channel width.
R = 2.0; // Channel central radius.

h = 0.1; // Mesh size.

Point(1) = {0,  L/2, 0, h};
Point(2) = {0, -L/2, 0, h};

Line(1)  = {1, 2};

Extrude { {0, 0, 1}, {0, -R, 0}, -Pi / 2 } { Line{1}; }

Physical Line(0) = {1};
Physical Line(1) = {2};
Physical Line(2) = {3};
Physical Line(3) = {4};
Physical Surface(10) = {5};

Mesh 2;
Save "mesh-pipe.msh";