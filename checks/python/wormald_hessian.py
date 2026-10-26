"""Wormald's Hessian determinant H for regular tournaments, computed two ways in exact arithmetic.

(1) From the definition (Wormald, "Tournaments with many Hamilton cycles", Theorem 4): with
    B(y) = I - (I/2 + A) diag(y) and K = det B, the function y_d = h(y_1, ..., y_{d-1}) solves K = 0
    near z = (2/d, ..., 2/d), and H is the determinant of the Hessian of
    f(theta) = log h(z e^{i theta_1}, ..., z e^{i theta_{d-1}}) at theta = 0.  Put G(u) = K(z e^u).
    Since K is affine in each y_j, G_j = z K_j, G_jj = G_j and G_jk = z^2 K_jk (j != k), where K_j and
    K_jk are determinants of B(z) with column j (and k) replaced by -(I/2 + A) e_j.  Implicit
    differentiation of u_d = psi(u_1, ..., u_{d-1}) gives psi_j = -G_j/G_d and
        psi_jk = -(G_jk + G_jd psi_k + G_kd psi_j + G_dd psi_j psi_k) / G_d,
    and f(theta) = log z + psi(i theta), so H = det(-Hess psi) = (-1)^(d-1) det(Hess psi).
(2) From Theorem 1.1 of the paper: with S = A - A^T,
        d / H = det(d^2 I - S^2) / det(d^2 I + S^2) = (q(d^2) / q(-d^2))^2,
    where det(tI - S) = t q(t^2).

    usage: python3 -I wormald_hessian.py CLASSES_FILE [CLASSES_FILE ...]

Each CLASSES_FILE is the output of regtour.c (one tournament per skew spectrum).  The script also
checks Wormald's table (R3, ..., R13, R9', the 7-vertex example) and random regular tournaments on
11 and 13 vertices, and exits with status 1 on any mismatch.
"""
import random
import sys
from fractions import Fraction as Fr
from math import isqrt


def det(M):
    M = [row[:] for row in M]
    n, d = len(M), Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            if f:
                for k in range(c, n):
                    M[r][k] -= f * M[c][k]
    return d


def hessian_H(A):
    """H from Wormald's definition (implicit differentiation, exact)."""
    d = len(A)
    z = Fr(2, d)
    col = lambda j: [(1 if i == j else 0) - z * (Fr(1, 2) * (i == j) + A[i][j]) for i in range(d)]
    dcol = lambda j: [-(Fr(1, 2) * (i == j) + A[i][j]) for i in range(d)]
    base = [col(j) for j in range(d)]

    def detcols(cols):
        return det([[cols[j][i] for j in range(d)] for i in range(d)])

    assert detcols(base) == 0
    K1 = []
    for j in range(d):
        c = base[:]
        c[j] = dcol(j)
        K1.append(detcols(c))
    K2 = [[Fr(0)] * d for _ in range(d)]
    for j in range(d):
        for k in range(j + 1, d):
            c = base[:]
            c[j], c[k] = dcol(j), dcol(k)
            K2[j][k] = K2[k][j] = detcols(c)
    G1 = [z * K1[j] for j in range(d)]
    G2 = [[(G1[j] if j == k else z * z * K2[j][k]) for k in range(d)] for j in range(d)]
    D = d - 1
    assert G1[D] != 0
    psi = [-G1[j] / G1[D] for j in range(D)]
    Hpsi = [[-(G2[j][k] + G2[j][D] * psi[k] + G2[k][D] * psi[j] + G2[D][D] * psi[j] * psi[k]) / G1[D]
             for k in range(D)] for j in range(D)]
    return (-1) ** D * det(Hpsi)


def closed_form(A):
    """(d/H, rho) from Theorem 1.1: d/H = det(d^2 I - S^2)/det(d^2 I + S^2) = rho^2."""
    d = len(A)
    S = [[Fr(A[i][j] - A[j][i]) for j in range(d)] for i in range(d)]
    S2 = [[sum(S[i][k] * S[k][j] for k in range(d)) for j in range(d)] for i in range(d)]
    minus = det([[d * d * (i == j) - S2[i][j] for j in range(d)] for i in range(d)])
    plus = det([[d * d * (i == j) + S2[i][j] for j in range(d)] for i in range(d)])
    rho_num = det([[d * (i == j) - S[i][j] for j in range(d)] for i in range(d)])
    rho_den = isqrt(int(plus))
    assert rho_den * rho_den == plus, "det(d^2 I + S^2) is not a square"
    assert rho_num * rho_num == minus
    return minus / plus, Fr(rho_num) / rho_den


def circ(d, D):
    return [[1 if (j - i) % d in D else 0 for j in range(d)] for i in range(d)]


def is_regular_tournament(A):
    d = len(A)
    return (all(A[i][i] == 0 for i in range(d))
            and all(A[i][j] + A[j][i] == 1 for i in range(d) for j in range(d) if i != j)
            and all(sum(r) == (d - 1) // 2 for r in A))


def random_regular(d, steps, rng):
    A = circ(d, set(range(1, (d + 1) // 2)))
    for _ in range(steps):  # reversing a directed triangle keeps every score
        i, j, k = rng.sample(range(d), 3)
        if A[i][j] and A[j][k] and A[k][i]:
            A[i][j] = A[j][k] = A[k][i] = 0
            A[j][i] = A[k][j] = A[i][k] = 1
    return A


def blowup_R3_R3():
    """R9' of Wormald's table: three copies of R3 arranged along R3."""
    A = [[0] * 9 for _ in range(9)]
    for a in range(3):
        for b in range(3):
            for i in range(3):
                for j in range(3):
                    u, v = 3 * a + i, 3 * b + j
                    if a == b:
                        A[u][v] = 1 if (j - i) % 3 == 1 else 0
                    else:
                        A[u][v] = 1 if (b - a) % 3 == 1 else 0
    return A


def wormald_7():
    """The 7-vertex tournament of Wormald's paper with d/H = (3088/1287)^2: cyclic triples abc and
    a'b'c', arcs aa', bb', cc', ad, bd, cd, the rest forced by regularity."""
    a, b, c, a1, b1, c1, dd = range(7)
    A = [[0] * 7 for _ in range(7)]
    for u, v in [(a, b), (b, c), (c, a), (a1, b1), (b1, c1), (c1, a1), (a, a1), (b, b1), (c, c1),
                 (a, dd), (b, dd), (c, dd)]:
        A[u][v] = 1
    # complete to a regular tournament: try all orientations of the remaining pairs
    rest = [(u, v) for u in range(7) for v in range(u + 1, 7) if A[u][v] == 0 and A[v][u] == 0]
    for mask in range(1 << len(rest)):
        B = [row[:] for row in A]
        for t, (u, v) in enumerate(rest):
            if mask >> t & 1:
                B[u][v] = 1
            else:
                B[v][u] = 1
        if is_regular_tournament(B):
            return B
    raise AssertionError("no completion")


TABLE = [  # Wormald's table: sqrt(d/H)
    ("R3", circ(3, {1}), Fr(2)),
    ("R5", circ(5, {1, 2}), Fr(44, 19)),
    ("R7", circ(7, {1, 2, 3}), Fr(3032, 1231)),
    ("R9", circ(9, {1, 2, 3, 4}), Fr(437584, 171769)),
    ("R9'", blowup_R3_R3(), Fr(5488, 2197)),
    ("R11", circ(11, {1, 2, 3, 4, 5}), Fr(108859552, 41833969)),
    ("R13", circ(13, {1, 2, 3, 4, 5, 6}), Fr(41637757888, 15769379963)),
    ("7-vertex example", wormald_7(), Fr(3088, 1287)),
]


def check(name, A, expect=None):
    assert is_regular_tournament(A), name
    d = len(A)
    H = hessian_H(A)
    dH, rho = closed_form(A)
    ok = H != 0 and Fr(d) / H == dH and dH == rho * rho and (expect is None or rho == expect)
    print(f"{name}: d={d} d/H={Fr(d) / H} rho={rho} {'ok' if ok else 'MISMATCH'}", flush=True)
    return ok


def classes(path):
    for line in open(path):
        t = line.split()
        if not t or t[0] != "poly":
            continue
        d = (len(t) - 5) // 2
        yield int(t[d + 3]), [[int(ch) for ch in r] for r in t[d + 5:]]


def main():
    bad = total = 0
    for name, A, rho in TABLE:
        total += 1
        bad += not check(name, A, rho)
    for path in sys.argv[1:]:
        for cnt, A in classes(path):
            total += 1
            bad += not check(f"class of {cnt} labelled tournaments", A)
    rng = random.Random(20261014)
    for d, cnt in ((11, 12), (13, 4)):
        for t in range(cnt):
            total += 1
            bad += not check(f"random regular #{t}", random_regular(d, 400, rng))
    print(f"checked {total} tournaments, mismatches {bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
