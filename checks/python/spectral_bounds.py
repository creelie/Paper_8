"""Checks of the spectral consequences of Theorem 1.1 (Theorem 1.5 of the paper), in floating point
with a safety margin, plus Wormald's 9-vertex circulant values.

For a regular tournament with skew matrix S = A - A^T let mu_1 >= ... >= mu_m (m = (d-1)/2) be the
moduli of the eigenvalues of S in conjugate pairs, s_k = mu_k^2.  Theorem 1.1 gives
    rho = sqrt(d/H) = prod_k (d^2 + s_k)/(d^2 - s_k).
Checked here, for every regular tournament class on d <= 9 vertices (files from regtour.c) and for
random regular tournaments on 11, 13, 15 and 21 vertices:
  * sum_k s_k = d(d-1)/2 and mu_1 <= cot(pi/2d);
  * Theorem 1.5: rho >= ((d+1)/(d-1))^m, with equality exactly for the doubly regular classes
    (S S^T = dI - J), and e^{1/d} rho > e;
  * Theorem 1.5: rho <= U_d, the maximum of prod (d^2+s_k)/(d^2-s_k) under the two constraints;
  * the cyclic tournament R_d has the largest rho among the classes on d <= 9 vertices.
It also recomputes the seven 9-vertex "circulant" values in Wormald's paper and shows that six of the
connection sets D contain both r and 9 - r, so they do not define tournaments.

    usage: python3 -I spectral_bounds.py CLASSES_FILE [CLASSES_FILE ...]
"""
import math
import random
import sys

import numpy as np


def circ(d, D):
    return [[1 if (j - i) % d in D else 0 for j in range(d)] for i in range(d)]


def pairs(A):
    d = len(A)
    S = np.array([[A[i][j] - A[j][i] for j in range(d)] for i in range(d)], float)
    ev = sorted(np.linalg.eigvalsh(S @ S.T), reverse=True)
    return [max(ev[2 * k], 0.0) for k in range((d - 1) // 2)], S


def rho_of(s, d):
    return math.prod((d * d + x) / (d * d - x) for x in s)


def upper_bound(d):
    m, sigma, smax = (d - 1) // 2, d * (d - 1) / 2, 1 / math.tan(math.pi / (2 * d)) ** 2
    t = min(int(sigma // smax), m)
    r = sigma - t * smax
    phi = lambda x: math.log((d * d + x) / (d * d - x))
    return math.exp(t * phi(smax) + (phi(r) if t < m else 0.0))


def random_regular(d, steps, rng):
    A = circ(d, set(range(1, (d + 1) // 2)))
    for _ in range(steps):
        i, j, k = rng.sample(range(d), 3)
        if A[i][j] and A[j][k] and A[k][i]:
            A[i][j] = A[j][k] = A[k][i] = 0
            A[j][i] = A[k][j] = A[i][k] = 1
    return A


def classes(path):
    for line in open(path):
        t = line.split()
        if t and t[0] == "poly":
            d = (len(t) - 5) // 2
            yield [[int(ch) for ch in r] for r in t[d + 5:]]


def main():
    bad = 0
    groups = {}
    for path in sys.argv[1:]:
        for A in classes(path):
            groups.setdefault(len(A), []).append(("class", A))
    rng = random.Random(7)
    for d, cnt in ((11, 40), (13, 20), (15, 10), (21, 4)):
        groups.setdefault(d, []).extend(("random", random_regular(d, 50 * d * d, rng)) for _ in range(cnt))
    for d in sorted(groups):
        m = (d - 1) // 2
        low, up = ((d + 1) / (d - 1)) ** m, upper_bound(d)
        cyc = rho_of(pairs(circ(d, set(range(1, m + 1))))[0], d)
        best = 0.0
        for kind, A in groups[d]:
            s, S = pairs(A)
            r = rho_of(s, d)
            best = max(best, r) if kind == "class" else best
            doubly = np.allclose(S @ S.T, d * np.eye(d) - np.ones((d, d)))
            ok = (abs(sum(s) - d * (d - 1) / 2) < 1e-6
                  and math.sqrt(s[0]) <= 1 / math.tan(math.pi / (2 * d)) + 1e-9
                  and r >= low - 1e-12 and (abs(r - low) < 1e-9) == doubly
                  and math.exp(1 / d) * r > math.e
                  and r <= up + 1e-12 and r <= cyc + 1e-12)
            bad += not ok
        n_cls = sum(1 for k, _ in groups[d] if k == "class")
        line = (f"d={d}: {len(groups[d])} tournaments ({n_cls} classes), lower bound {low:.6f}, "
                f"R_d {cyc:.6f}, upper bound U_d {up:.6f}")
        if n_cls:
            line += f", largest class value {best:.6f} ({'R_d' if abs(best - cyc) < 1e-12 else 'NOT R_d'})"
            bad += abs(best - cyc) > 1e-12
        print(line, flush=True)
    a, b = 4 / math.pi ** 2, 0.5 - 4 / math.pi ** 2
    print(f"lim U_d = (pi^2+4)(3pi^2-8)/((pi^2-4)(pi^2+8)) = {(1 + a) * (1 + b) / ((1 - a) * (1 - b)):.9f}; "
          f"cosh(1)/cos(1) = {math.cosh(1) / math.cos(1):.9f}")
    for d in (101, 1001, 10001):
        print(f"U_{d} = {upper_bound(d):.9f}")
    # Wormald's 9-vertex circulants
    print("Wormald's 9-vertex circulants (value of E[H(T(T*))]/g(dn) in the paper):")
    for D, paper in (({1, 3, 4, 5}, 1.4582), ({1, 2, 4, 5}, 1.6811), ({1, 2, 3, 5}, 2.4557),
                     ({2, 3, 4, 5}, 1.5618), ({1, 3, 5, 6}, 1.5529), ({1, 4, 5, 6}, 1.9325),
                     ({2, 3, 5, 6}, 1.8117)):
        A = circ(9, D)
        tour = all(A[i][j] + A[j][i] == 1 for i in range(9) for j in range(9) if i != j)
        both = sorted(r for r in D if 9 - r in D)
        print(f"  D={sorted(D)}: paper {paper}, "
              + ("a tournament" if tour else f"not a tournament (contains r and 9-r for r in {both})"))
    print("failures", bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
