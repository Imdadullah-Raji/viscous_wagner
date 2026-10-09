# Re = 500, angle of attack = 30 degrees: flow animation

`re500_aoa30_q_vorticity.mp4` is a 1920 × 1080 H.264 video, 25 fps, 450 frames, 18 s. The PNG is a representative still for a slide or video poster.

This case is two-dimensional (frontAndBack is an empty boundary). The visualization therefore shows filled regions where Q c²/U∞² > 0.10, outlined at that threshold and colored by signed spanwise vorticity ωz c/U∞. It is not a three-dimensional isosurface. Red is positive (counterclockwise); blue is negative (clockwise). The color scale is fixed at ±8; extended colorbar ends indicate saturation.

Coordinates and velocity components are rotated by +30 degrees into the freestream-aligned frame using flowkit. The displayed crop is −0.75 ≤ xs/c ≤ 5.75 and −1.9 ≤ ys/c ≤ 1.9. The airfoil outline uses the ordered wall face centers. The case specifies c = 1, U∞ = 1 and ν = 0.002.

The final 450 consecutive saved velocity snapshots are used. Exact times are recorded in metadata.json. Playback advances one saved snapshot per frame (approximately 2.43 convective time units per playback second); no temporal interpolation or artificial loop closure is applied.

Q and vorticity are derived for visualization from saved U, rather than read from solver-generated diagnostic fields. Native cell-center velocity gradients use inverse-distance-squared weighted least squares over 12 neighbors. Neighbor connections crossing the airfoil are excluded. Q = −½ tr[(∇U)²]; ωz = ∂v/∂x − ∂u/∂y. Derived scalars are linearly interpolated onto a 0.0125c display grid. Solid-crossing interpolation triangles and a one-cell margin are masked. The immediate wall region is consequently excluded; this is intended for vortex visualization, not wall-gradient measurement. No spatial smoothing is applied.

Recreate with:

```bash
conda run -n flowkit python render_flow.py --out .
```

Use `--preview` to render only the still; `--case PATH` selects the input case. The script needs flowkit, fluidfoam, NumPy, SciPy, xarray, Matplotlib and ffmpeg.

Validation: the gradient/Q implementation was checked against analytic solid-body rotation (Q = +4, ωz = +4) and pure strain (Q = −4, ωz = 0). Video encoding, dimensions, duration and frame count were checked with ffprobe, with a full decode check using ffmpeg.
