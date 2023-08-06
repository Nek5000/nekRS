// geometrical parameter
r_ratio = 0.7;
L = 2*Pi;
d = 1.0;

// mesh parameter
n_axial = 36;
n_azimuthal = 40;
n_radial = 12;
radial_growth_rate = 1.9;

//
r_i = d*r_ratio/(1-r_ratio);
r_o = d/(1-r_ratio);

// 
Point(1) = {r_o,0,0};
Point(2) = {0,r_o,0};
Point(3) = {-r_o,0,0};
Point(4) = {0,-r_o,0};

Point(5) =  {(r_i+r_o)/2,0,0};
Point(6) = {0,(r_i+r_o)/2,0};
Point(7) = {-(r_i+r_o)/2,0,0};
Point(8) = {0,-(r_i+r_o)/2,0};

// 
Point(9) = {r_i,0,0};
Point(10) = {0,r_i,0};
Point(11) = {-r_i,0,0};
Point(12) = {0,-r_i,0};

// org point
Point(21) = {0,0,0};

//
Circle(1) = {1,21,2};
Circle(2) = {2,21,3};
Circle(3) = {3,21,4};
Circle(4) = {4,21,1};
//
Circle(5) = {5,21,6};
Circle(6) = {6,21,7};
Circle(7) = {7,21,8};
Circle(8) = {8,21,5};
//
Circle(9) = {9,21,10};
Circle(10) = {10,21,11};
Circle(11) = {11,21,12};
Circle(12) = {12,21,9};
//

Line(21) = {1,5};
Line(22) = {2,6};
Line(23) = {3,7};
Line(24) = {4,8};
//
Line(25) = {5,9};
Line(26) = {6,10};
Line(27) = {7,11};
Line(28) = {8,12};

//
Line Loop(1) = {1,22,-5,-21};
Line Loop(2) = {2,23,-6,-22};
Line Loop(3) = {3,24,-7,-23};
Line Loop(4) = {4,21,-8,-24};

//
Line Loop(5) = {5,26,-9,-25};
Line Loop(6) = {6,27,-10,-26};
Line Loop(7) = {7,28,-11,-27};
Line Loop(8) = {8,25,-12,-28};
//

For iface In {1:8:1}
	Plane Surface(iface) = {iface};
	Transfinite Surface {iface};
	Recombine Surface {iface};
EndFor

For icurve In {1:12:1}
	Transfinite Curve{icurve} = (n_azimuthal/4)+1;
EndFor

For icurve In {21:24:1}
	Transfinite Curve{icurve} = n_radial/2+1 Using Progression radial_growth_rate;
EndFor

For icurve In {25:28:1}
	Transfinite Curve{icurve} = n_radial/2+1 Using Progression 1/radial_growth_rate;
EndFor


For iface In {1:8:1}
Extrude {0,0,L}{
	Surface {iface};
	Layers{n_axial};
	Recombine;
}
EndFor

Physical Surface("innerWall", 1) = {155,177,199,133};
Physical Surface("outerWall", 2) = {59,37,81,103};
Physical Surface("inflow",3) = {1,2,3,4,5,6,7,8};
Physical Surface("outflow",4) = {72,160,50,138,182,94,204,116};
Physical Volume("fluid",5) = {1,2,3,4,5,6,7,8};

Mesh 1;
Mesh 2;
Mesh 3;

SetOrder 2;

Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;
Mesh.SaveAll = 0;
Mesh.Binary = 1;

Save "tcf.msh";
