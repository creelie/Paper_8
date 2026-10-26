#!/usr/bin/env bash
# run_all.sh: every check in this repository. Outputs go to build/.
#   C:      checks/c/regtour.c lists the regular tournaments on d = 3, 5, 7, 9 vertices with
#           N+(0) = {1..m} and must reproduce checks/data/classes*.txt (one line per
#           characteristic polynomial of S); Table 1 of the paper is read off these files.
#   Python: wormald_hessian.py computes H from Wormald's definition in exact arithmetic and
#           compares it with Theorem 1.1; hessian_matrix_bounds.py compares the whole Hessian
#           matrix with E^T F E and checks Corollary 3.3; spectral_bounds.py checks Theorem 1.5
#           and Remark 1.6.
#   GP:     wormald.gp repeats the comparison with Wormald's definition symbolically and checks
#           the closed form for R_d against the eigenvalue product to 80 digits.
# Needs a C compiler, Python 3 with numpy, and PARI/GP.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
B="$ROOT/build"
PY="$ROOT/checks/python"
mkdir -p "$B"
cd "$ROOT"
CLASSES="checks/data/classes3.txt checks/data/classes5.txt checks/data/classes7.txt checks/data/classes9.txt"

echo "== C: regular tournaments on at most 9 vertices"
cc -O2 -Wall -Wextra -Werror -o "$B/regtour" checks/c/regtour.c
for d in 3 5 7 9; do
  "$B/regtour" "$d" > "$B/classes$d.txt"
  cmp "$B/classes$d.txt" "checks/data/classes$d.txt"
  head -1 "$B/classes$d.txt"
done

echo "== Python: H from Wormald's definition against Theorem 1.1"
python3 -I "$PY/wormald_hessian.py" $CLASSES | tee "$B/wormald_hessian.txt"

echo "== Python: the full Hessian matrix and Corollary 3.3"
python3 "$PY/hessian_matrix_bounds.py" | tee "$B/hessian_matrix_bounds.txt"

echo "== Python: Theorem 1.5 and Remark 1.6"
python3 -I "$PY/spectral_bounds.py" $CLASSES | tee "$B/spectral_bounds.txt"

echo "== PARI/GP: independent check"
gp -q checks/gp/wormald.gp | tee "$B/wormald_gp.txt"
grep -q '^failures 0$' "$B/wormald_gp.txt"
echo "all checks passed"
