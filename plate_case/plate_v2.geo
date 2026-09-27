// ===========================================================================
//  Eldredge & Wang (AIAA 2010-4281) plate, Re = 1000
//  Structured O-grid around the plate, Delaunay beyond.
//  Pitch axis at the leading edge (their Fig. 3 condition).
// ===========================================================================
SetFactory("OpenCASCADE");

c    = 1.0;         // chord
th   = 0.023;       // thickness / chord
r    = th*c/2;      // radius of the semicircular edges = 0.0115
xp   = 0.0;         // PITCH AXIS
d    = 0.5*c;       // thickness of the structured block (see note at the end)
R    = 30.0*c;      // far field
span = 1.0;

// Radial: 39 nodes at ratio 1.12 over 0.5 gives a first cell of 8.2e-4,
// i.e. the same wall spacing as before, and the first 20 layers reproduce
// the old boundary-layer field exactly.
nRad = 39;  gRad = 1.12;
// Circumferential: inner and outer sides of each sector must agree, which
// is automatic here because the outer curve is a parallel offset of the inner.
nArc = 25;  nFlank = 181;

// ---- plate outline: a stadium ---------------------------------------------
pLEc=newp; Point(pLEc)={r,0,0};      pTEc=newp; Point(pTEc)={c-r,0,0};   // arc centres
pLE =newp; Point(pLE) ={0,0,0};      pTE =newp; Point(pTE) ={c,0,0};
pTL =newp; Point(pTL) ={r,r,0};      pTR =newp; Point(pTR) ={c-r,r,0};
pBL =newp; Point(pBL) ={r,-r,0};     pBR =newp; Point(pBR) ={c-r,-r,0};

cLEu=newl; Circle(cLEu)={pLE,pLEc,pTL};   cTop=newl; Line(cTop)={pTL,pTR};
cTEu=newl; Circle(cTEu)={pTR,pTEc,pTE};   cTEl=newl; Circle(cTEl)={pTE,pTEc,pBR};
cBot=newl; Line(cBot)={pBR,pBL};          cLEl=newl; Circle(cLEl)={pBL,pLEc,pLE};

// ---- block boundary: a PARALLEL OFFSET of the plate ------------------------
// Same two arc centres, radius r+d, same flank length.  This is the whole
// trick: offsetting a stadium gives a stadium, so corresponding inner and
// outer points share a normal and the radial mesh lines come out orthogonal.
// A circle here instead pushes the average non-orthogonality from 4 to 22.
ro = r + d;
qW =newp; Point(qW) ={r-ro,   0,  0};   qE =newp; Point(qE) ={c-r+ro, 0,  0};
qTL=newp; Point(qTL)={r,      ro, 0};   qTR=newp; Point(qTR)={c-r,    ro, 0};
qBL=newp; Point(qBL)={r,     -ro, 0};   qBR=newp; Point(qBR)={c-r,   -ro, 0};

oLEu=newl; Circle(oLEu)={qW,pLEc,qTL};    oTop=newl; Line(oTop)={qTL,qTR};
oTEu=newl; Circle(oTEu)={qTR,pTEc,qE};    oTEl=newl; Circle(oTEl)={qE,pTEc,qBR};
oBot=newl; Line(oBot)={qBR,qBL};          oLEl=newl; Circle(oLEl)={qBL,pLEc,qW};

// ---- the two cuts that make the ring four-sided ----------------------------
// Drawn OUTWARD from the plate: a Progression runs from a curve's first point,
// so this puts the fine spacing at the wall and not at the block boundary.
spLE=newl; Line(spLE)={pLE,qW};    spTE=newl; Line(spTE)={pTE,qE};

// ---- the two sectors -------------------------------------------------------
llUp=newll; Curve Loop(llUp)={cLEu,cTop,cTEu, spTE, -oTEu,-oTop,-oLEu, -spLE};
sUp =news;  Plane Surface(sUp)={llUp};
llLo=newll; Curve Loop(llLo)={cTEl,cBot,cLEl, spLE, -oLEl,-oBot,-oTEl, -spTE};
sLo =news;  Plane Surface(sLo)={llLo};

// ---- far field: circle about the PITCH AXIS, block as a hole ---------------
// Concentric with the rotation axis so the domain maps onto itself as the mesh
// turns.  The block boundary is internal and is free to be a stadium.
cFar =newl;  Circle(cFar)={xp,0,0,R};
llFar=newll; Curve Loop(llFar)={cFar};
llBlk=newll; Curve Loop(llBlk)={oLEu,oTop,oTEu,oTEl,oBot,oLEl};
sFar =news;  Plane Surface(sFar)={llFar,llBlk};

// ---- structured meshing of the block ---------------------------------------
Transfinite Curve {cLEu,cTEu,cTEl,cLEl, oLEu,oTEu,oTEl,oLEl} = nArc;
Transfinite Curve {cTop,cBot, oTop,oBot} = nFlank Using Bump 0.2;
Transfinite Curve {spLE,spTE} = nRad Using Progression gRad;
Transfinite Surface {sUp} = {pLE,pTE,qE,qW};
Transfinite Surface {sLo} = {pTE,pLE,qW,qE};
Recombine Surface {sUp,sLo};

// ---- sizing outside the block ----------------------------------------------
// Field 2 starts at the block's outer spacing so there is no jump at the
// interface.  Field 3 holds the wake resolved; this is what was missing before
// and why the vortex looked faceted in the animation.
Field[1]=Distance;  Field[1].CurvesList={oLEu,oTop,oTEu,oTEl,oBot,oLEl};
Field[2]=Threshold; Field[2].InField=1;
Field[2].SizeMin=0.014; Field[2].SizeMax=1.5; Field[2].DistMin=0; Field[2].DistMax=22;
Field[3]=Ball;
Field[3].XCenter=0; Field[3].YCenter=0; Field[3].ZCenter=0;
Field[3].Radius=3.0; Field[3].VIn=0.03; Field[3].VOut=1.5; Field[3].Thickness=1.0;
Field[4]=Min;  Field[4].FieldsList={2,3};
Background Field=4;

Mesh.MeshSizeExtendFromBoundary=0;   // or the wall spacing leaks into the far field
Mesh.MeshSizeFromPoints=0;
Mesh.MeshSizeFromCurvature=0;
Mesh.Algorithm=5;        // Delaunay, for the unstructured part
Mesh.RecombineAll=0;     // only the sectors are recombined, above
Mesh.Smoothing=10;

// ---- 2D -> one cell in z ---------------------------------------------------
// Extruded in ONE call so the sectors and the far field share their interface
// faces instead of duplicating them.  Index layout per surface is
// [top, volume, laterals in curve-loop order]; sUp and sLo have 8 laterals
// each, sFar has 7 (the far circle, then the 6 block curves).
out[] = Extrude {0,0,span} { Surface{sUp,sLo,sFar}; Layers{1}; Recombine; };

Physical Volume ("internal")     = {out[1], out[11], out[21]};
Physical Surface("frontAndBack") = {sUp, sLo, sFar, out[0], out[10], out[20]};
Physical Surface("farfield")     = {out[22]};
Physical Surface("plate")        = {out[2],out[3],out[4], out[12],out[13],out[14]};
