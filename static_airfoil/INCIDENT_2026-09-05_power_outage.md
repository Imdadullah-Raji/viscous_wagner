# Incident: power outages during the Re=500 sweep

Two events, 2026-09-05 (AoA 22) and 2026-09-07 (AoA 16). Sections 1-9 are the
first; section 10 is the second and the two verification lessons it added.

**Date:** 2026-09-05, 22:36 -- 23:34 (down 56 min, recovery 2 min)
**Impact:** ~35 min of compute re-run. No data lost, no data corrupted.
**Latent bug found:** phase-2 resume would have silently recorded the wrong
number of orbits. Fixed.

---

## 1. Timeline

| time | event |
|---|---|
| 21:15:58 | `run_even.sh` launched: AoA 22 (warm-started from AoA 20), then AoA 24 |
| 21:16:18 | phase 0 starts |
| 21:26 | phase 0 done, `Co_max=0.025918`, dt set to 5.4016e-4 |
| 22:01 | chunk to t=20.5 completes, writes its restart point |
| ~22:36 | **power lost** mid-chunk, ~t=40.18 of the run to t=40.5 |
| 23:32 | power returns, machine boots |
| 23:34:11 | resumed from t=20.5 |

`last -x reboot` reports the **boot**, not the loss, so it showed 23:32 and looked
like the outage happened an hour after the run stopped. The interrupted chunk log's
mtime (22:36) is what actually dates the event. Worth remembering: to date a crash,
read the mtime of whatever the process was writing, not the boot record.

## 2. What survived, and why

Everything except the chunk in flight. That is a property of the staged design, not
luck: **phase 1 writes exactly one field snapshot per chunk, at the chunk end, purely
as a restart point.** The worst case for any interruption is therefore re-running the
chunk in progress -- at most `CHUNK=20` time units, here 35 minutes.

Recovery cost had the design been "no writes at all until convergence": the entire
transient from t=0.5, and the warm start seeded into time 0 as well. For AoA 30 that
would have been about 8 hours.

## 3. Integrity verification (before resuming, not after)

A killed process leaves a clean filesystem. A power cut can leave a half-written file,
so the restart point was checked rather than assumed:

| check | result |
|---|---|
| t=20.5 present on all 6 processors | yes -- `U p phi U_0 phi_0 Cp wallShearStress uniform/` on each |
| per-file byte sizes across ranks | proportional to per-rank cell counts, no short file |
| `uniform/time` present per rank | yes (without it the restart has no time state) |
| interrupted t=40.5 chunk | no `End` in log, no field directory -- nothing half-written |
| `re500_aoa20` 1026 snapshots | 0 problems |
| `re500_aoa30` 1289 snapshots | 0 problems |
| zero-byte files anywhere | only `re500_aoa30.OpenFOAM`, a ParaView marker meant to be empty |
| truncated `forceCoeffs` .dat | parses -- `limit_cycle.py` reads the series to t=40.18 |

The force file was left mid-line by the outage. `load_series()` survives it because it
`float()`s each field and skips rows that are short or unparseable; it does not assume
the last line is complete. That was already there for the restart-seam problem and
turned out to cover this too.

## 4. "If it was phase 1, why was there any write at all?"

Fair question, and the naming does not help -- the driver prints *"PHASE 1: transient,
**field-free**, chunks of 20"* immediately before writing a field.

The decision was that **the POD dataset** starts after the transient. That is enforced,
and provably so. What phase 1 writes is a restart point, which is a different thing
serving a different purpose:

| | phase 1 chunk write | phase 2 snapshot |
|---|---|---|
| how often | 1 per **37 025** steps (AoA 22) | 1 per **223** steps |
| purpose | restart point | POD sample |
| reconstructed? | **no** -- stays in `processor*/` | yes |
| in the POD set? | **no** | yes |
| cost | 6.2 MB per write, ~31 MB per case | ~13 MB per write, 8-16 GB per case |

The mechanism is `writeInterval == endTime`: a chunk's write interval is set equal to
its end time, so exactly one write can occur, on the last step. That is what "field-free"
means in the log -- no field *series*, one restart snapshot. Phase 0 is the same: its
"(no fields)" line also ends with a single write at t=0.5.

Exclusion from the POD set is structural, not a filter someone has to remember to
apply. `reconstructPar` runs only for `t >= t_conv`, so the earlier chunk writes are
never reconstructed and cannot enter the dataset. Verified on the finished AoA 30 case:

```
phase-1 chunk writes, decomposed only: 20.5 40.5 60.5 80.5 100.5
reconstructed below t_conv=120.5:      0
reconstructed at or above t_conv:      1289
```

So the answer is: yes, phase 1 writes -- one snapshot per chunk, 166x rarer than phase
2, 0.2% of the case's storage, and structurally barred from the POD set. It bought a
35-minute recovery instead of an 8-hour one.

(One related caveat, unchanged by this incident: `t_conv` itself **is** a chunk write,
and it sits off the phase-2 sampling grid -- 0.0898 before the first sample instead of
0.0972. It gets reconstructed because the range is inclusive. Drop it from the POD set.)

## 5. The latent bug the outage exposed

Phase-1 chunk boundaries and phase-2 snapshots are both just time directories.
**Nothing on disk distinguishes them.** The resume logic read the latest time
directory and assumed it was a chunk boundary:

```bash
LATEST=$(ls -d processor0/[0-9]* | sed 's|.*/||' | sort -g | tail -1)
if [ -n "$LATEST" ] && [ "$LATEST" > "$T0" ]; then
  T=$LATEST                      # <-- treated as a phase-1 chunk boundary
  ...                            #     re-enters the phase-1 chunk loop
```

Had the outage landed 90 minutes later, during phase 2, the resume would have:

1. re-entered the phase-1 loop from a phase-2 snapshot time;
2. called `limit_cycle.py`, which would immediately return `CONVERGED=yes` -- the flow
   *is* converged;
3. set `T_CONV` to the **interruption time** rather than the real one;
4. run phase 2 for `NCYCLES` orbits from there, appending a second, overlapping block
   of snapshots on a different grid offset.

The result would have been a POD set with the wrong orbit count and a seam in its
sampling -- and, critically, **no error anywhere**. Same failure signature as the
period-doubling bug: a plausible number, silently wrong.

## 6. The fix

`run_staged.sh` now stamps a marker as phase 2 begins:

```bash
printf 'T_CONV=%s\nDT=%s\nWI=%s\nT2=%s\nORDER=%s\nPERIOD=%s\nSAMPLE_PERIOD=%s\n' \
       "$T_CONV" "$DT" "$WI" "$T2" "$ORDER" "$PERIOD" "$SAMPLE_PERIOD" > .phase2
```

and the resume branch reads it before anything else:

```bash
if [ -f .phase2 ] && [ -n "$LATEST" ]; then
  . ./.phase2
  if [ "$LATEST" > "$T_CONV" ]; then
    # finish phase 2 on the ORIGINAL dt, writeInterval and endTime
    patch_cd latestTime $DT $T2 runTime $WI
    run log.phase2_resume
    sawtooth_gate "$LATEST" "phase 2 resume"
    reconstructPar -newTimes -time "${T_CONV}:"
    exit 0
  fi
fi
```

`dt` and `writeInterval` are deliberately re-used rather than re-derived, for the same
reason `extend_phase2.sh` holds them fixed: snapshots either side of the interruption
must stay on one uniform grid. A re-derived interval shifts by a step or two and puts
a seam in what is meant to be a uniformly weighted ensemble.

## 7. Operational changes

- **`setsid nohup ... &`**, not plain `nohup`. The job needs its own session or it dies
  with whatever launched it.
- **`run_even.sh` is now restart-safe.** It builds and warm-starts an angle only if
  that case has no `processor0` yet. This matters more than it sounds: `make_case.sh`
  opens with `rm -rf` on the case directory, so blindly re-running the original chain
  script after an interruption would have destroyed both the progress *and* the warm
  start -- turning a 35-minute loss into a full re-run.

## 8. Residual risk

**AoA 22's own phase 2 will not write the `.phase2` marker.** The patch was installed
via temp-file + `mv` so the running process's inode stayed untouched -- which means
that process is still executing the pre-patch copy. The marker will be written by hand
when AoA 22 reaches phase 2. AoA 24 gets it automatically, since `run_even.sh` invokes
the script fresh for each angle.

Editing a running bash script in place is itself a hazard worth naming: bash reads the
file incrementally by file offset, so an in-place rewrite (anything that truncates and
rewrites the same inode, `python open(p,'w')` included) can make a running script
resume mid-token in the wrong place. Write a new file and `mv` it over.

## 9. Checklist for the next one

1. `pgrep -c foamRun` -- confirm nothing is still running.
2. Date the event from the mtime of the interrupted log, not from `last`.
3. Verify the restart point on **every** processor: field set, byte sizes proportional
   to per-rank cell counts, `uniform/time` present.
4. Sweep finished cases for zero-byte and short files before trusting them.
5. Check whether `.phase2` exists -- it decides which resume path is correct.
6. Relaunch the restart-safe chain script with `setsid nohup`. Never re-run
   `make_case.sh` on a case that has progress.


## 10. Second occurrence, 2026-09-07 (AoA 16)

**Down 12:03:55, back 12:45:55, resumed 12:49:20.** ~43 min of compute lost, no data
lost. AoA 16 was mid-chunk t=20.5 -> 40.5 of phase 1; AoA 18 had finished cleanly at
11:10 and was untouched.

The §9 checklist worked. What it did not cover was the corruption mode, which was
different from the first outage's, and two verification steps that were themselves
wrong.

### The corruption: NUL padding, not a torn line

On 2026-09-05 the force file was left **mid-line** -- an ordinary short write, and
`load_series()` already tolerated it. This time
`postProcessing/forceCoeffs/20.5/forceCoeffs.dat` ended in **80 NUL bytes** after the
last good sample at t=22.220256: 247547 bytes on disk, 247467 of them real. That is the
journalling signature of a power cut -- the file's *length* was committed while the data
blocks never reached the platter -- and it is worth recognising on sight, because a NUL
run does not look like truncation in an editor and `tail` prints it as blank lines.

`load_series()` would in fact have survived it (`float()` on a NUL field raises
`ValueError` and the row is skipped), but the directory was quarantined to
`.interrupted/` anyway rather than trusted, on the same reasoning as the diverged-chunk
quarantine in §15 of `FINDINGS.md`: the restart re-runs t=20.5 -> 40.5 and appends to
that same directory, so leaving damaged bytes in it buys nothing and risks a hybrid
line.

### Verification lesson 1: `grep -P '\x00'` found nothing, and was wrong

`grep` on this box is **ugrep**, and `grep -qP '\x00'` matched none of the 47 force
files including the one that was visibly full of NULs. The scan reported a clean bill of
health for a corrupt file. Byte-level integrity checks go through Python
(`b'\x00' in open(f,'rb').read()`), not through `grep`, on this machine.

### Verification lesson 2: the fields are BINARY, so NUL bytes there mean nothing

Re-run in Python, the same scan flagged **36932 files** as containing NUL bytes --
including `constant/polyMesh/points`, `faces` and `owner`. None of it was corruption.
`controlDict` sets `writeFormat binary`, so NUL bytes are the normal contents of a field
file, and scanning them for NULs is meaningless by construction. Only the ASCII
`forceCoeffs` .dat and the other `postProcessing` text output can be checked that way.

**The valid truncation test for a binary OpenFOAM field is its footer.** A complete file
ends with the `// *****` banner; a write cut off by a power loss does not. Applied to
AoA 16:

| write | files with a valid footer | without |
|---|---|---|
| t=20.5 (the restart point) | 42 | 0 |
| t=0.5 | 42 | 0 |
| t=0 | 24 | 0 |

across all six processors -- so the restart point was provably whole, not merely
plausible.

### Checklist additions

7. Check the ASCII `postProcessing` files for a **NUL tail** as well as a short last
   line, and do it in Python -- ugrep's `-P '\x00'` silently matches nothing here.
8. Do **not** scan binary field files for NUL bytes; `writeFormat` is `binary` and the
   result is meaningless. Check the `// *****` footer instead, on every rank.
9. Quarantine the interrupted chunk's `postProcessing/forceCoeffs/<restart>/` to
   `.interrupted/` before relaunching. The restart appends to that same directory.
