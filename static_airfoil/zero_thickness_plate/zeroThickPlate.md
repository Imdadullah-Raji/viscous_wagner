# Zero Thickness Plate

To test the effect of finite nose radius on a3 and to set up the step response simulations and compare against classical wagner result a zero thickness flat plate is desired.

## Where the Code Stands

I wrote an initial gmsh geofile. It is currently
- unparameterized, h is defined with Re but has no effect. Also other numerical parameters too are set by hand within the script.
- Also the tip resolution question.I have placed a bump and tried to match the first cell layer on the radial direction. But one problem stands, *what is the length scale of the resolution? the nose has zero radius* 
- The code is still 2D, for OpenFOAM we need to extrude.I also intend to test it out in Nektar++ later, to be logged, not to be acted upon now. (Nektar++ is in a docker container on this machine)



---

## Python mesher (added by Claude, 2026-10-18)

`ztp_mesher.py` ports `zeroThickPlate.geo` to the gmsh Python API. The `.geo` stays unchanged as the
original. The block topology is the same; the plate is the internal interface.

```bash
/home/raji/miniconda3/envs/research/bin/python ztp_mesher.py out.msh --re 500 [--tip 5.6e-4 --nplate 240 --nspokes 140]
BAFFLES=$PWD/createBafflesDict ../make_steady.sh <case> 500 <AoA> out.msh 6000    # from static_airfoil/
```

**What changed from the `.geo`:**
- Chord on [0, 1], far field at 30c centred at 0.25, as on the other meshes.
- h1 = 0.05/√Re sets the first cell normal to the plate.
- The tip cell is an explicit length (`--tip`). Each plate edge is split at mid-chord, and each half
  grows away from its tip.
- The first spoke cell matches the cell it continues.
- Extruded one cell in z. Physical groups: `internal`, `farfield`, `frontAndBack`, and `plate` (internal).

**Conversion:**
- `gmshToFoam` turns the internal `plate` faces into a faceZone.
- `createBafflesDict` (here) splits it into the walls `plate_upper` and `plate_lower`, both in group
  `airfoil`. So `make_case.sh`'s `airfoil` boundary conditions, forces, Cp and wallShearStress cover both
  sides.
- `make_case.sh` runs `createBaffles` only when `BAFFLES=` is set. The default path is byte-identical to
  before (checked).

**Mesh check (corrected 2026-10-18).** All meshes pass `checkMesh`, all hexes. An earlier version of
this line gave 64° and 1.13 for T3; those were T1's numbers. The real values:

| mesh | max non-orthogonality | max skewness | max aspect ratio |
|---|---:|---:|---:|
| T1 | 64.3° | 1.13 | 9.2 |
| T3 (the sweep) | 73.6° | 2.10 | 36.9 |
| T4 | 76.2° | 2.39 | 73.8 |

**Where the non-orthogonality is.** It sits in the four side blocks, the wedges upstream and downstream
of the plate, and peaks along the diagonal spokes.
- On T4, about 25 % of faces exceed 60°.
- Close to the tips: 34° within 0.005c, 54° within 0.02c, 63° within 0.05c.
- The blocks above and below the plate are close to orthogonal.
- The cause is the block topology: each wedge maps a short vertical block edge onto a 45° far-field
  arc.
- Centring the far field on mid-chord (76.1°) or making the far-field spacing uniform (77.0°) does
  not help.

**Tip resolution: still a convergence parameter.** See `A3_RE500.md` for the 6° tip study (T1–T4).
- The order-of-magnitude viscous scales at the edge are c/Re at zero incidence and c/(Re·α)² at
  incidence α. At Re 500 they span about 2e-3 down to 4e-4 c over 0–6°.
- The sweep uses T3 (tip 5.6e-4c).

## Elliptic-coordinate mesher (added by Claude, 2026-10-19)

`ell_mesher.py` is the mesher to use. It builds confocal ellipses and hyperbolas with foci at the
tips, the conformal map z = c/2 + (c/2)·cosh(μ + iν). It is orthogonal everywhere (max
non-orthogonality about 5°) and clusters at the tips by itself.

- Usage:
  ```bash
  python3 ell_mesher.py out.msh --re 500 --nnu 800 --nmu 210
  BAFFLES=$PWD/createBafflesDict ../make_steady.sh <case> 500 <AoA> out.msh 6000
  ```
- The plate is still an internal surface, so the baffle workflow is unchanged.
- Only the four cells touching the tips are skewed: each has a 180° corner where the map doubles
  angles. Their skewness falls with nν, and nν = 800 passes `checkMesh`.
- The H-type `ztp_mesher.py` above is superseded. Its tip distortion produced a spurious
  non-cubic lift curve (see `A3_RE500.md`).
