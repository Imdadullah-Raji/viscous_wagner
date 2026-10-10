# Sample case: plate, Re = 500, 4 deg

Dictionaries only, from `re500/steady_aoa4`. Built by `run_sweep.py --re 500` (which calls
`static_airfoil/make_steady.sh`), with 6000 SIMPLE iterations. The mesh is
`plate_mesher.py <out.msh> --re 500` (63,360 hexes, first cell 0.05/sqrt(Re)).
