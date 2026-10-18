/*
The domain: a zero thickness flat plate 
in  the middle. A block above it and below it.The 
interface becomes  the plate using OpenFOAMs baffles
utility. 6 radial connectors , 8 blocks in total. Four
connectors to the 45 degree points on the circle on
each quadrant. Two from the leading and trailing edge.
*/  

// ----------------------------------//
//------- Domain Art ----------------//
// ----------------------------------//

/*          _________________
    /                                \
 /        " "                  " "    \
|            " "__ __ _ _ _ " "         |
|___ _ __ ___ __|_________|__ _ ___ __ _|
|               |_ _ _ _ _|             |
 \           " "            " "        / 
    \                               /
          " "_________________" "
*/

//-----------------------------------//
//-------Geometric Parameters--------//
//-----------------------------------//

lc = 0.1 ;           //characteristic length
a = 0.5 ;            //half chord
R = 40*a ;           //farfield radius

//-----------------------------------//
//-----------Mesh Parameters---------//
//-----------------------------------//

Re = 1000;           // flow reynolds number
b = a/2 ;            // upper and lower block height
r_side = 1.02 ;       // inflation layer for the BL
h = 0.05/Sqrt(Re) ;  // first cell height 
r_block = 1.02  ;     // within block progression ratio
r_far = 1.03  ;       // farfield progression ratio
Nside = 101;          // number of points on block sides
Nplate = 101;         // number of points on the plate 
Nspokes = 200;
beta = 0.075;            // tip resolution

origin = newp;
Point(origin) = {0, 0, 0, lc};

//-----------The block---------------//

left = newp; Point(left) = {-a, 0 , 0, lc};
right = newp; Point(right) = {a,0, 0, lc};
right_up = newp; Point(right_up) = {a, b, 0, lc};
left_up = newp; Point(left_up) = {-a, b, 0, lc};
left_down = newp; Point(left_down) = {-a, -b, 0, lc};
right_down = newp; Point(right_down) = {a, -b, 0, lc};


//-----------Far field points--------//

left_far = newp;
Point(left_far) = {-R, 0, 0, lc};

right_far = newp;
Point(right_far) = {R, 0, 0, lc};

right_up_far = newp;
Point(right_up_far) = {R/Sqrt(2), R/Sqrt(2), 0, lc};

left_up_far = newp;
Point(left_up_far) = {-R/Sqrt(2), R/Sqrt(2), 0, lc};

left_down_far = newp;
Point(left_down_far) = {-R/Sqrt(2), -R/Sqrt(2), 0, lc};

right_down_far = newp;
Point(right_down_far) = {R/Sqrt(2), -R/Sqrt(2), 0, lc};

//---- lines and arcs -------//

block_upper = newl;
Line(block_upper) = {left_up, right_up};

block_lower = newl;
Line(block_lower) = {left_down, right_down};

block_interface = newl;
Line(block_interface)={left, right};

block_up_right = newl;
Line(block_up_right) = {right, right_up}; 

block_up_left = newl;
Line(block_up_left) = {left, left_up};

block_down_right = newl;
Line(block_down_right) = {right, right_down}; 

block_down_left = newl;
Line(block_down_left) = {left, left_down};

//---- six circular arcs. clockwise-------//
arc_up = newl;
Circle(arc_up) = {left_up_far, origin, right_up_far};

arc_right_up = newl;
Circle(arc_right_up) = {right_up_far, origin,right_far};

arc_right_down = newl;
Circle(arc_right_down) = {right_far, origin,right_down_far };

arc_down = newl;
Circle(arc_down) = {right_down_far, origin,left_down_far};

arc_left_down = newl;
Circle(arc_left_down) = {left_down_far, origin, left_far};

arc_left_up = newl;
Circle(arc_left_up) = {left_far, origin, left_up_far};

//----------radial spokes--------------//

quadrant1_radial = newl;
Line(quadrant1_radial) = {right_up, right_up_far};

quadrant2_radial = newl;
Line(quadrant2_radial) = {left_up, left_up_far};

quadrant3_radial = newl; 
Line(quadrant3_radial) = {left_down, left_down_far};

quadrant4_radial = newl;
Line(quadrant4_radial)={right_down, right_down_far};

left_radial = newl;
Line(left_radial) = {left, left_far};

right_radial = newl; 
Line(right_radial) = {right, right_far};

upper_loop = newll;
Curve Loop(upper_loop) = {
    quadrant2_radial,
    arc_up,
    -quadrant1_radial,
    -block_upper
};

upper_surface = news;
Plane Surface(upper_surface) = {upper_loop};

upper_block_loop =  newll;
Curve Loop(upper_block_loop) = {
    block_upper,
    -block_up_right,
    -block_interface,
    block_up_left
};
upper_block_surface = news;
Plane Surface(upper_block_surface) = {upper_block_loop};

lower_block_loop = newll;
Curve Loop(lower_block_loop) = {
    block_lower,          // LD -> RD
    -block_down_right,    // RD -> R
    -block_interface,     // R -> L
    block_down_left       // L -> LD
};

lower_block_surface = news;
Plane Surface(lower_block_surface) = {lower_block_loop};


//---------- Lower outer block ----------//

lower_loop = newll;
Curve Loop(lower_loop) = {
    block_lower,
    quadrant4_radial,
    arc_down,
    -quadrant3_radial
};

lower_surface = news;
Plane Surface(lower_surface) = {lower_loop};


//---------- Upper-right block ----------//

upper_right_loop = newll;
Curve Loop(upper_right_loop) = {
    quadrant1_radial,
    arc_right_up,
    -right_radial,
    block_up_right
};

upper_right_surface = news;
Plane Surface(upper_right_surface) = {upper_right_loop};


//---------- Lower-right block ----------//

lower_right_loop = newll;
Curve Loop(lower_right_loop) = {
    quadrant4_radial,
    -arc_right_down,
    -right_radial,
    block_down_right
};

lower_right_surface = news;
Plane Surface(lower_right_surface) = {lower_right_loop};


//---------- Upper-left block ----------//

upper_left_loop = newll;
Curve Loop(upper_left_loop) = {
    quadrant2_radial,
    -arc_left_up,
    -left_radial,
    block_up_left
};

upper_left_surface = news;
Plane Surface(upper_left_surface) = {upper_left_loop};


//---------- Lower-left block ----------//

lower_left_loop = newll;
Curve Loop(lower_left_loop) = {
    quadrant3_radial,
    arc_left_down,
    -left_radial,
    block_down_left
};

lower_left_surface = news;
Plane Surface(lower_left_surface) = {lower_left_loop};


Transfinite Curve{
    arc_up,
    block_upper,
    block_interface,
    block_lower,
    arc_down
} = Nplate Using Bump beta;

Transfinite Curve{
    block_up_left,
    block_up_right,
    block_down_left,
    block_down_right,
    -arc_right_up,
    arc_right_down,
    arc_left_up,
    -arc_left_down
} = Nside Using Progression r_side;


Transfinite Curve{
    left_radial,
    right_radial,
    quadrant1_radial,
    quadrant2_radial,
    quadrant3_radial,
    quadrant4_radial
}=Nspokes Using Progression r_far;

Transfinite Surface{
    lower_block_surface,
    upper_block_surface,
    upper_surface,
    lower_surface,
    upper_right_surface,
    upper_left_surface, 
    lower_right_surface,
    lower_left_surface
};


Recombine Surface {
    upper_surface,
    lower_surface,
    upper_block_surface,
    lower_block_surface,
    upper_right_surface,
    lower_right_surface,
    upper_left_surface,
    lower_left_surface
}; //to generate quad elements

Physical Curve("farfield") = {
    arc_up,
    arc_right_up,
    arc_right_down,
    arc_down,
    arc_left_down,
    arc_left_up
};

Physical Curve("plate_interface") = {
    block_interface
};

Physical Surface("fluid") = {
    upper_surface,
    lower_surface,
    upper_block_surface,
    lower_block_surface,
    upper_right_surface,
    lower_right_surface,
    upper_left_surface,
    lower_left_surface
};

Mesh 2;