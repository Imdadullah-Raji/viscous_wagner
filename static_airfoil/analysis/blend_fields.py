#!/usr/bin/env python3
"""out = A + eps*(B - A)  on two OpenFOAM ASCII fields sharing a mesh.

WHY THIS EXISTS.  Near the bifurcation point the saturated amplitude collapses,
so the standard seed -- map the base flow from half a degree up -- stops being a
small perturbation.  At 11.5 deg it would be dCl ~ 0.0053 against an expected
saturated amplitude of 0.0018-0.0030: 2-3x TOO BIG.  Starting above the attractor
makes the perturbation DECAY toward it, which reads as sigma < 0 and would
wrongly call the angle stable -- the exact sign error the run is meant to settle.

Blending with eps ~ 0.1 gives a perturbation an order of magnitude below the
attractor, so what is measured is growth away from the fixed point, not
relaxation onto the limit cycle.

Rewrites every `nonuniform List<scalar|vector>` block (internalField and any
boundaryField `value` entries), leaving uniform entries, BC types and headers
untouched.  The two files must come from the same mesh -- block lengths are
checked, not assumed.
"""
import re
import sys

BLOCK = re.compile(
    r"(nonuniform\s+List<(scalar|vector)>\s*\n\s*(\d+)\s*\n\()(.*?)(\n\)\s*;)",
    re.S)


def numbers(body, kind):
    if kind == "scalar":
        return [[float(x)] for x in body.split()]
    return [[float(v) for v in m.split()]
            for m in re.findall(r"\(([^)]*)\)", body)]


def render(vals, kind):
    if kind == "scalar":
        return "\n" + "\n".join(f"{v[0]:.10g}" for v in vals)
    return "\n" + "\n".join("(" + " ".join(f"{c:.10g}" for c in v) + ")" for v in vals)


def main():
    if len(sys.argv) != 5:
        sys.exit("usage: blend_fields.py <fileA> <fileB> <eps> <out>")
    fa, fb, eps, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
    sa, sb = open(fa).read(), open(fb).read()
    ba, bb = list(BLOCK.finditer(sa)), list(BLOCK.finditer(sb))
    if len(ba) != len(bb):
        sys.exit(f"block count differs: {len(ba)} vs {len(bb)} -- different mesh?")

    pieces, last, nblend, nval = [], 0, 0, 0
    for ma, mb in zip(ba, bb):
        if ma.group(3) != mb.group(3) or ma.group(2) != mb.group(2):
            sys.exit(f"block mismatch: {ma.group(2)}[{ma.group(3)}] "
                     f"vs {mb.group(2)}[{mb.group(3)}]")
        kind = ma.group(2)
        va, vb = numbers(ma.group(4), kind), numbers(mb.group(4), kind)
        if len(va) != len(vb) or len(va) != int(ma.group(3)):
            sys.exit(f"length mismatch in a {kind} block")
        blended = [[a + eps * (b - a) for a, b in zip(ra, rb)]
                   for ra, rb in zip(va, vb)]
        pieces.append(sa[last:ma.start()])
        pieces.append(ma.group(1) + render(blended, kind) + ma.group(5))
        last = ma.end()
        nblend += 1
        nval += len(va)
    pieces.append(sa[last:])
    open(out, "w").write("".join(pieces))
    print(f"  {out}: blended {nblend} block(s), {nval} values, eps={eps}")


if __name__ == "__main__":
    main()
