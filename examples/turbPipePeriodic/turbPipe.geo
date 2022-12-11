/*
   ***** gmsh Script for generating 2D/3D mesh for straight pipe. *****
   ** Saleh Rezaeiravesh, salehr@kth.se
   >>> To generate mesh run: gmsh turbPipe.geo -3 -order 2   
*/

//constants
PI=3.14159265359;

meshDim=3; 

// GRID SETTINGS ///////////////////////////////////////

// Geometrical parameters
// Note: r<RB<R
R=0.5;   //Pipe radius
r=0.35;
RB=0.48917;   
lambda=0.6;   //=R_{arc}/R
Lz=6;   //length in z-dir

// Grid Paramaters
Nc=8;   // no. of nodes (=#elem+1) in azimuthal direction
NB=2;   // no. of elemtns adjacent to the wall
NM=5;   // no. of nodes (=#elem+1) between the near wall layer and central square part

compressRatio_B=0.85;  //ratio of grid compression toward the wall (<1)
compressRatio_M=0.87;  //compression ratio in the middle layer
Nz=35;  //no of elements in z-dir (streamwise)

///////////////////////////////////////////////////

th=PI/4.;  //theta
dx=r*Cos(th);
dy=r*Sin(th);
dxB=RB*Cos(th);
dyB=RB*Sin(th);
Dx=R*Cos(th);
Dy=R*Sin(th);

//***** define points coordinates
//auxiliary points (only help define the geometry)
Point(1) = {0, 0, 0, 1.0};
Point(2) = {lambda*R, 0, 0, 1.0};
Point(3) = {0, -lambda*R, 0, 1.0};
Point(4) = {-lambda*R, 0, 0, 1.0};
Point(5) = {0 , lambda*R, 0, 1.0};

//blocks vertices
Point(6)={dx, dy, 0.0, 1.0};
Point(7)={dx, -dy, 0.0, 1.0};
Point(8)={-dx,-dy, 0.0, 1.0};
Point(9)={-dx, dy, 0.0, 1.0};
Point(10)={dxB, dyB, 0.0, 1.0};
Point(11)={dxB, -dyB, 0.0, 1.0};
Point(12)={-dxB,-dyB, 0.0, 1.0};
Point(13)={-dxB, dyB, 0.0, 1.0};
Point(14)={Dx, Dy, 0.0, 1.0};
Point(15)={Dx, -Dy, 0.0, 1.0};
Point(16)={-Dx,-Dy, 0.0, 1.0};
Point(17)={-Dx, Dy, 0.0, 1.0};

//***** define lines and curves
Circle(1)={9, 3, 6};   //Circle()={startNode, circleCenter, endNode}
Circle(2)={6, 4, 7};   
Circle(3)={7, 5, 8};   
Circle(4)={8, 2, 9};   
Circle(5)={13, 1, 10};
Circle(6)={10, 1, 11};
Circle(7)={11, 1, 12};
Circle(8)={12, 1, 13};
Circle(9)={17, 1, 14};
Circle(10)={14, 1, 15};
Circle(11)={15, 1, 16};
Circle(12)={16, 1, 17};
Line(13)={6, 10};
Line(14)={7, 11};
Line(15)={8, 12};
Line(16)={9, 13};
Line(17)={10, 14};
Line(18)={11, 15};
Line(19)={12, 16};
Line(20)={13, 17};

//***** assign number of mesh on the created lines/arcs
Transfinite Line {1, 2, 3, 4} = Nc Using Bump 1.0;   
Transfinite Line {5, 6, 7, 8} = Nc;   
Transfinite Line {9, 10, 11, 12} = Nc;   
Transfinite Line {13, 14, 15, 16} = NM Using Progression compressRatio_M;  
Transfinite Line {17, 18, 19, 20} = NB Using Progression compressRatio_B;   //over the radial lines near the wall. Note: "For example Progression 2 meaning that each line element in the series will be twice as long as the preceding one)".

//***** create surfaces
// Note: use a negative sign if a line is swept in the opposite direction of the original definition
Line Loop(1)={1, 2, 3 ,   4};   Plane Surface(1)={1}; //central part of the mesh
Line Loop(2)={5, -13, -1, 16};  Plane Surface(2)={2}; 
Line Loop(3)={13, 6, -14 , -2}; Plane Surface(3)={3}; 
Line Loop(4)={-3, 14, 7, -15};  Plane Surface(4)={4}; 
Line Loop(5)={-16, -4, 15, 8};  Plane Surface(5)={5}; 
Line Loop(6)={9, -17, -5, 20};  Plane Surface(6)={6}; 
Line Loop(7)={17, 10, -18, -6}; Plane Surface(7)={7}; 
Line Loop(8)={-7, 18, 11, -19}; Plane Surface(8)={8}; 
Line Loop(9)={-8, 19, 12, -20}; Plane Surface(9)={9}; 

If (meshDim==2)
  //Line Loop (50)={9,10,11,12};
   Physical Line("wall")={9, 10, 11, 12};
// Physical Line("wall")={50};
   Physical Surface(1)={1:9};
EndIf

Recombine Surface "*";
Transfinite Surface "*";

If (meshDim==3)

   //make a 3d mesh by extrusion in z-dir
   mesh3D[]=Extrude {0,0,Lz} 
   {
       Surface{1:9};
       Layers{Nz}; 
       Recombine; 
   };

   //Physical Surfaces & Volume (Note: gmsh only generates mesh for the physical entities)
   // 1: inlet
   // 2: outlet
   // 3: wall
   Physical Surface("inlet") = {6, 7, 8, 9, 5, 2, 3, 4, 1};
   Physical Surface("outlet") = {152, 174, 196, 218, 64, 86, 108, 130, 42};
   Physical Surface("wall") = {139, 213, 191, 165};
   Physical Volume("flowDomain") = {6, 2, 1, 4, 8, 3, 7, 5, 9};

   Recombine Volume "*";

EndIf

Coherence;

Mesh 1;
Mesh 2;
Mesh 3;

SetOrder 2;

Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;
Mesh.SaveAll = 0;
Mesh.Binary = 0;

Save "turbPipe.msh";
