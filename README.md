# Two conjectures of Wormald on Hamilton cycles in tournaments

Deep Bhattacharjee. Final revision, 26 October 2026.

Wormald (*Tournaments with many Hamilton cycles*, arXiv:2610.16245, on his web page
since 2013) bounds the largest number h(n) of Hamilton cycles in a tournament on n
vertices from below by averaging over random blow-ups of a fixed regular tournament T on
d vertices. The bound is governed by a number H, the Hessian determinant of a function
defined implicitly by T, and he conjectured that H is never zero and that d/H is always
the square of a rational number (Conjecture 6), and that his constants for the cyclic
tournaments R_d have a limit (Conjecture 2). The paper proves both: with S the skew
adjacency matrix of T,

    H = d det(d^2 I + S^2) / det(d^2 I - S^2),   sqrt(d/H) = det(dI - S) / det(dI - iS),

and for R_d

    sqrt(d/H) = ((d+1)^d + (d-1)^d) / ((d+i)^d + (d-i)^d)  ->  cosh 1 / cos 1 = 2.8559578925...

which is the constant 2.855958 that Wormald estimated by fitting polynomials. The paper
also bounds sqrt(d/H) over all regular tournaments on d vertices. Every statement is
proved by hand; the programs here only check the identities on small cases, and no
statement of the paper rests on them.

The results were first made public on 14 October 2026 as release v1.0.0 of
[creelie/Breakthrough-problems-1](https://github.com/creelie/Breakthrough-problems-1)
(later versions archived on Zenodo as
[10.5281/zenodo.23567153](https://doi.org/10.5281/zenodo.23567153) and
[10.5281/zenodo.23631103](https://doi.org/10.5281/zenodo.23631103), concept DOI
[10.5281/zenodo.23567152](https://doi.org/10.5281/zenodo.23567152)); an earlier and
incomplete version of the paper appeared on Hexagon as
[hexagon:2610.00460](https://hexagonmath.org/2610.00460) on 22 October 2026.
Two later works overlap: Divoux and Wang (arXiv:2610.23755, 16 October 2026) prove
h(n) ~ (cosh 1 / cos 1) g(n) by a different method, and Liu and Ye
(github.com/LStar404/wormald-tournament-conjectures, 19 October 2026) treat Wormald's
Conjectures 1, 2 and 6 independently.

## Layout

- `paper/main.tex`, `paper/figures/`: the LaTeX source and the two figures; `paper/main.pdf` the compiled paper (12 pages, single column).
- `dist/`: the paper as `WormaldHamiltonCycles_final2.pdf`, its source zip and the arXiv tarball.
- `checks/c/regtour.c`: lists the regular tournaments on at most 9 vertices by the characteristic polynomial of S; its output is `checks/data/classes*.txt`, the source of Table 1.
- `checks/python/wormald_hessian.py`: computes H from Wormald's definition in exact rational arithmetic and compares it with Theorem 1.1 on 41 tournaments.
- `checks/python/hessian_matrix_bounds.py`: compares the whole Hessian matrix with E^T F E and checks the bounds of Corollary 3.3.
- `checks/python/spectral_bounds.py`: checks Theorem 1.5 and Remark 1.6, with random regular tournaments up to 21 vertices.
- `checks/gp/wormald.gp`: the same comparison in PARI/GP, and the closed form for R_d against the eigenvalue product to 80 digits for d <= 61.
- `run_all.sh`: runs all of them.

## Running the checks

    ./run_all.sh

Requirements: a C compiler, python3 with numpy, and PARI/GP.

## Releases

Releases are published by the workflow `release.yml` from the Actions tab.
Each release carries the paper files of `dist/`.

| Release | DOI |
| --- | --- |
| v1.0.0 | to be recorded |

## Author

Deep Bhattacharjee, Formerly, Electro-Gravitational Space Propulsion Laboratory (EGSPL), Bhubaneswar, Odisha 751030, India; Formerly, University of Calcutta, Senate House, 87/1 College Street, Kolkata 700073, West Bengal, India.
Emails: d.bhattacharjee@erl-forschung.de, itsdeep@live.com, egspl@proton.me (for questions about this repository).
ORCID: 0000-0003-0466-750X.

## Licence

MIT, see `LICENSE`.
