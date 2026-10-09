#!/bin/bash
#
# Expand one-line OpenFOAM FoamFile headers so fluidfoam can parse them.
#
#     FoamFile { version 2.0; format ascii; class volScalarField; object p; }
# ->  FoamFile
#     {
#         version 2.0;
#         ...
#     }
#
# fluidfoam's header parser (readof.py:_parse_session) is a line-by-line state
# machine: it needs a line that is exactly "FoamFile", then one that is exactly
# "{", then one "key value;" per line. Given the compact form it finds no
# header at all, decides the file is raw data, and dies with a confusing
# "IndexError: index out of range" further down.
#
# Usage:  ./edit.sh [-n] [case-glob]      -n = dry run
#         ./edit.sh -n                    show what would change
#         ./edit.sh                       fix NACA0012/re500_aoa* (default)
#         ./edit.sh 'NACA0012/bi_point_search/*'   some other set

set -uo pipefail

DRY=0
[ "${1-}" = "-n" ] && { DRY=1; shift; }
GLOB="${1-NACA0012/re500_aoa*}"

# Top-level 0/ only. -maxdepth 3 keeps processor*/0/* out (those are depth 4
# and are written by decomposePar in the standard banner format anyway).
mapfile -t FILES < <(find . -maxdepth 3 -type f -path "./$GLOB/0/*" \
                          \( -name p -o -name U \) | sort)

if [ ${#FILES[@]} -eq 0 ]; then
    echo "No 0/p or 0/U under ./$GLOB — nothing to do."
    exit 1
fi

changed=0
skipped=0
for f in "${FILES[@]}"; do
    # Idempotent: only the compact one-line form needs touching.
    if ! grep -qE '^FoamFile[[:space:]]*\{.*\}[[:space:]]*$' "$f"; then
        skipped=$((skipped + 1))
        continue
    fi
    if [ "$DRY" -eq 1 ]; then
        echo "  would fix  $f"
    else
        cp -p "$f" "$f.bak"
        sed -i -E '/^FoamFile[[:space:]]*\{.*\}[[:space:]]*$/{
            s/^FoamFile[ \t]*\{[ \t]*/FoamFile\n{\n    /
            s/[ \t]*\}[ \t]*$/\n}/
            s/;[ \t]*/;\n    /g
            s/[ \t]+\n/\n/g
            s/\n[ \t]*\n\}/\n}/
        }' "$f"
        echo "  fixed      $f   (backup: $(basename "$f").bak)"
    fi
    changed=$((changed + 1))
done

echo
if [ "$DRY" -eq 1 ]; then
    echo "dry run: $changed file(s) would change, $skipped already fine"
else
    echo "$changed file(s) changed, $skipped already fine"
fi
