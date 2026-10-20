# Repo copy: host names, logins and paths replaced by $HOST_B, $REMOTE_B_DIR, $CFD_USER (set them to run).
"""progress.py -- progress bars for the zero-thickness-plate queues, local and remote, read-only.

The runs are headless background jobs; this only READS their logs, so it can be started and stopped
(Ctrl-C) at any time without touching them.  Per queue: one bar for the chain (cases done, with the
running case counted fractionally) and one for the running case's iterations, an s/iteration measured
from the solver's own ExecutionTime, and the ETA of the queue.

  python progress.py              one snapshot (works through `!` in Claude Code)
  python progress.py --follow     live, redrawn every 30 s (remote machines every 60 s); Ctrl-C to stop

Edit QUEUES when a queue is added.  Remote machines are read over ssh with one call per refresh.
"""
import os
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

try:
    from tqdm import tqdm
    def bar(n, total, prefix):
        return tqdm.format_meter(n, total, 0, ncols=110, prefix=prefix,
                                 bar_format="{l_bar}{bar}| {n:.2f}/{total:g}")
except ImportError:                                                  # system python has no tqdm
    def bar(n, total, prefix):
        w = 40; k = int(w * n / total) if total else 0
        return f"{prefix}{100 * n / total:3.0f}%|{'█' * k}{' ' * (w - k)}| {n:.2f}/{total:g}"

HERE = Path(__file__).resolve().parent
CFD_USER = os.environ.get("CFD_USER", "set-CFD_USER")       # cfd-test login (OS Login, sudo)
CFD = ["ssh", "-o", "BatchMode=yes", "-l", CFD_USER, "cfd-test", "sudo", "bash", "-s"]
# label, ssh command (None = this machine), case folder, angles in run order, first and last iteration
QUEUES = [
    ("this machine  fold Re500, ascending", None, str(HERE / "fold/re500_asc"),
     [7, 8, 9, 9.5, 10, 10.5, 11, 12], 0, 10000),
    ("cfd-test      Re700 tight extension", CFD, "/home/raji/cases/ztp_tight/re700",
     [4, 4.25, 4.5], 9000, 12000),
    ("cfd-test      fold Re500, descending", CFD, "/home/raji/cases/ztp_fold/re500_desc",
     [12, 11, 10.5, 10, 9.5, 9], 0, 10000),
]

# For each case: last iteration, End present, s/iteration over the last ~50 iterations, log age (s).
PROBE = r'''
for C in "$@"; do L=$C/log.steady
  if [ -f "$L" ]; then
    tail -c 300000 "$L" | grep -a -E '^(Time = |ExecutionTime|End)' | awk -v c="$C" -v age=$(( $(date +%s) - $(stat -c %Y "$L") )) '
      /^Time = /{t=$3+0; ts[++n]=t} /^ExecutionTime/{e[n]=$3+0; m=n; if(!f) f=n} /^End/{end=1}
      END{k=(m-50>f?m-50:f); r=(m>k && ts[m]>ts[k])?(e[m]-e[k])/(ts[m]-ts[k]):0;
          printf "%s %s %d %.4f %d\n", c, t+0, end+0, r, age}'
  else echo "$C - 0 0 0"; fi
done
'''


def probe(ssh, dirs):
    script = "set -- " + " ".join(dirs) + "\n" + PROBE
    cmd = ssh if ssh else ["bash", "-s"]
    try:
        out = subprocess.run(cmd, input=script, capture_output=True, text=True, timeout=120).stdout
    except subprocess.TimeoutExpired:
        return None
    rows = {}
    for line in out.splitlines():
        f = line.split()
        if len(f) == 5 and f[0] in dirs:
            rows[f[0]] = (None if f[1] == "-" else float(f[1]), f[2] == "1", float(f[3]), int(f[4]))
    return rows


def render(cache):
    out = [datetime.now().strftime("%Y-%m-%d %H:%M:%S")]
    for label, ssh, folder, angles, it0, it1 in QUEUES:
        rows = cache.get(label)
        out.append(f"\n{label}")
        if rows is None:
            out.append("  (machine not reachable)"); continue
        done, cur, rate = 0, None, 0.0
        for a in angles:
            t, end, r, age = rows.get(f"{folder}/steady_aoa{a:g}", (None, False, 0, 0))
            if t is not None and end and t >= it1:
                done += 1
            elif t is not None and not end and t >= it0 and cur is None:
                cur, rate = (a, t, age), r
        frac = done + ((cur[1] - it0) / (it1 - it0) if cur else 0)
        out.append("  " + bar(frac, len(angles), f"cases {done}/{len(angles)} "))
        if cur:
            a, t, age = cur
            out.append("  " + bar(t - it0, it1 - it0, f"{a:>5g} deg   ").replace(".00/", "/"))
            left = (len(angles) - frac) * (it1 - it0) * rate
            eta = (datetime.now() + timedelta(seconds=left)).strftime("%a %H:%M") if rate else "?"
            stale = f"   (log silent {age // 60} min!)" if age > 600 else ""
            out.append(f"  {rate:.3f} s/it   queue left {timedelta(seconds=int(left))}   ETA {eta}{stale}")
        elif done == len(angles):
            out.append("  complete")
        else:
            out.append("  waiting (nothing running in this queue)")
    return "\n".join(out)


def main():
    follow = "--follow" in sys.argv
    cache, polled = {}, {}
    while True:
        now = time.time()
        for label, ssh, folder, angles, *_ in QUEUES:
            if label not in polled or now - polled[label] >= (60 if ssh else 30) or not follow:
                cache[label] = probe(ssh, [f"{folder}/steady_aoa{a:g}" for a in angles])
                polled[label] = now
        text = render(cache)
        if not follow:
            print(text); return
        print("\033[H\033[J" + text + "\n\n(Ctrl-C to stop; the runs are not affected)", flush=True)
        time.sleep(30)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        pass
