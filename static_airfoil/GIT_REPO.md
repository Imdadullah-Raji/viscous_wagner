# The git repository: how it is organised and how to keep it up to date

2026-10-10. This workspace (`/home/raji/Research/Thesis`) is **not** under version control. Code and notes are
copied by hand into a separate, clean repository:

- **Local:** `/home/raji/Research/Repositories/viscous-wagner`
- **Remote:** `git@github.com:Imdadullah-Raji/viscous_wagner.git` (public), branch `main`
- **Copied:** scripts, notes and the figures the notes use. **Not copied:** fields, force histories,
  meshes, logs, PDFs. Those stay here.

The history on GitHub is fixed. **Add new commits on top; never rewrite or force-push.**

## Layout: the repo mirrors the workspace

Each file sits at the **same relative path** as in `/home/raji/Research/Thesis`. Relative links in the notes
and `../plate_case` paths in the scripts therefore work in both places.

```
viscous-wagner/
├── README.md                   programme, main finding, machines, how the history was built, authorship
├── HYPOTHESIS_PROVENANCE.md    where each idea came from (user / ChatGPT / Claude)
├── .gitignore                  *.msh, processor*/, postProcessing/, polyMesh/, logs, PDFs
├── plan/simulation_plan.md
├── static_airfoil/             = Thesis/static_airfoil, code and notes only
│   ├── *.sh                    drivers (make_case, run_staged, make_steady, run_steady, ...)
│   ├── *.md                    CLAUDE.md, FINDINGS.md, STEADY_*.md, ...
│   ├── meshing/  analysis/  lsb/  instructions/  bizarre_experiments/
│   ├── analysis/figs/  lsb/figs/        figures (PNG only, no csv/txt outputs)
│   ├── reports/*/*.tex                  LaTeX sources only
│   ├── NACA0012/<sweep>/  NACA0006/<sweep>/   only the driver scripts and notes of each sweep
│   └── sample_cases/           dictionaries of a steady and an unsteady case (exists only in the repo)
└── plate_case/                 = Thesis/plate_case, code and notes only
    ├── plate_mesher, run_sweep.py/.sh, run_batch.sh, analysis/, mesh_conv_re500/, unsteady_re500/
    ├── re<RE>/RESULTS.md, re500/analysis/, re1000/analysis/     notes and scripts inside case folders
    └── sample_cases/           (exists only in the repo)
```

Five files differ from the workspace on purpose. All others are byte-identical copies.
- `README.md`, `HYPOTHESIS_PROVENANCE.md` and the two `sample_cases/README.md` exist only in the repo.
- `plate_case/run_batch.sh` was rebuilt from a session transcript, because cfd-test was unreachable.

## What goes in

| put it in | leave it out |
|---|---|
| new `.py` / `.sh` drivers and analysis scripts | case folders, `processor*/`, time directories, `postProcessing/` |
| notes (`*.md`) and LaTeX sources | `.msh` meshes (they regenerate from `mesher.py` / `plate_mesher.py`) |
| PNG figures that a note embeds or names | `.csv` / `.txt` script outputs, logs, PDFs |
| the `instructions/*.md` briefs | `*.bak`, `*.orig`, `.pre_*` copies, `__pycache__/` |
| driver scripts that live in a sweep folder (e.g. `NACA0006/steady_re_scaling/build.sh`) | `.vscode/` |

**Watch for scripts and notes inside case folders.** The first import missed `plate_case/re500/analysis/`
because it sat among the cases. The check below catches that.

## Routine: bring the repo up to date

```bash
T=/home/raji/Research/Thesis
R=/home/raji/Research/Repositories/viscous-wagner

# 1. What is new or changed? (code and notes; case folders skipped)
cd $T
find static_airfoil plate_case \
     \( -name 'processor*' -o -name 'steady_aoa*' -o -name 'aoa[0-9]*' -o -name 're500_aoa*' \
        -o -name 'ship_re500_odd_aoa' -o -name 'claude_md' -o -name '__pycache__' \) -prune -o \
     -type f \( -name '*.py' -o -name '*.sh' -o -name '*.md' -o -name '*.tex' -o -name '*.geo' \) -print |
while read f; do
  if [ ! -e "$R/$f" ]; then echo "NEW      $f"
  elif ! cmp -s "$f" "$R/$f"; then echo "CHANGED  $f"; fi
done

# 2. Copy what belongs (same relative path), for example:
cp --parents static_airfoil/analysis/new_script.py "$R/"
cp --parents static_airfoil/analysis/figs/new_dir/plot.png "$R/"

# 3. Commit (message: see below), then push
cd $R && git add -A && git status --short      # check nothing unwanted is staged
git commit && git push
```

The check prints nothing when the repo is up to date. As of 2026-10-10 it reports two lines:
- `NEW` for this file, `GIT_REPO.md`;
- `CHANGED` for `instructions/oct10_gitVersionControl.md`, which you edited after it was committed.

Neither has been copied yet. Copy them if you want them in the repo.

To find figures a note embeds that are missing from the repo:

```bash
cd $R
git grep -h -o -I -E '[A-Za-z0-9_./-]+\.png' -- '*.md' '*.tex' | sort -u    # then look for each in $T
```

## Commit messages

The convention so far:

1. **Subject: the research decision or finding, not the software change.** For example, "Keep only odd
   powers: the symmetric section makes C_L odd in alpha", not "add fit_odd_cubic.py". A real bug fix
   can still say so.
2. **Body:** what was decided or found, with the key numbers.
3. **An `Authorship:` block that names every file it adds or changes**, grouped by:
   - `Claude-generated (Claude Opus 5.5)`: Claude Code / Claude agent sessions
   - `Codex-generated`: OpenAI Codex (GPT 6 Astra)
   - `Written by the user`
   - `Figures (output of the scripts)`
4. **No hypothesis/idea sources in commits.** Those go in `HYPOTHESIS_PROVENANCE.md`; add a row when a new
   idea enters.
5. When Claude Code makes the commit, it appends its `Co-Authored-By:` / `Claude-Session:` trailer lines.

Template:

```
NACA0006: a3 plateau survives to 8 deg at Re 2000

<two or three lines: what was run, the result, the numbers>

Authorship:
- Claude-generated (Claude Opus 5.5): static_airfoil/analysis/foo.py, static_airfoil/FOO.md
- Written by the user: static_airfoil/instructions/oct12_foo.md
- Figures (output of the scripts): static_airfoil/analysis/figs/foo/a.png
```

**Dates.** Commits so far carry the date the research step happened, not the date they were committed.
For catching up on older work, set both dates when committing:
`GIT_AUTHOR_DATE="2026-10-12T15:00:00+0600" GIT_COMMITTER_DATE="$GIT_AUTHOR_DATE" git commit`.
For current work, the default date is right.

## Quirks of the existing history

- **Commits 1–47 were reconstructed on 2026-10-10.** Dates come from file times, run-log headers and
  session logs.
  - Where an older version survived (`.bak`, `.pre_*`), it was committed first.
  - Otherwise a file enters with its final content, so some early commits contain later edits.
- **The 2026-10-09 folder move is replayed in one commit**, "Group cases by body". Before it, files use
  the old flat paths (`steady_re500/`, `re500_aoa20/`, `steady_naca0006_re500/`).
- **Dates fall out of order after commit 47.** The two plate commits added later carry their Oct 6 and
  Oct 9 dates. Commit 47 has a nominal 12:00 timestamp, later than the commits after it.
- **Older commits say only `Codex-generated`.** The model (GPT 6 Astra) is named only from the README
  onward. That was a deliberate choice: no rewrite.
- **Known stale drivers:** `NACA0012/steady_re{100,1000,2000}/sweep.sh` still `cd ..`. That was right
  before the folder move; they now need `../..`. They are copied unchanged, both here and in the repo.
