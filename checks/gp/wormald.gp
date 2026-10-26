\\ wormald.gp -- an independent PARI/GP check of Theorem 1.1 and of the cyclic product formula.
\\
\\ (a) H from Wormald's definition, by symbolic differentiation of K(y) = det(I - (I/2 + A) diag(y))
\\     and implicit differentiation at z = (2/d, ..., 2/d); compared with
\\     d/H = det(d^2 I - S^2) / det(d^2 I + S^2) and with (q(d^2)/q(-d^2))^2, det(tI - S) = t q(t^2).
\\ (b) For the cyclic tournament R_d (d odd, 3 <= d <= 61) the exact value q(d^2)/|q(-d^2)| is
\\     compared with prod_{j odd <= d-2} (d^2 + cot(pi j/2d)^2)/(d^2 - cot(pi j/2d)^2) at 80 digits.
\\
\\ usage: gp -q wormald.gp  (prints one line per check and "failures N" at the end)

hessH(A) =
{
  my(d = #A, v = vector(d, i, eval(Str("y", i))), z = 2/d, B, K, K1, G1, G2, D = d, psi, Hp);
  my(at = P -> substvec(P, v, vector(d, i, z)));
  B = matid(d) - (matid(d)/2 + A) * matdiagonal(v);
  K = matdet(B);
  if (at(K) != 0, error("K(z) is not zero"));
  K1 = vector(d, j, deriv(K, v[j]));
  G1 = vector(d, j, z * at(K1[j]));
  G2 = matrix(d, d, j, k, if (j == k, z * at(K1[j]) + z^2 * at(deriv(K1[j], v[j])),
                                       z^2 * at(deriv(K1[j], v[k]))));
  psi = vector(d - 1, j, -G1[j] / G1[D]);
  Hp = matrix(d - 1, d - 1, j, k,
              -(G2[j, k] + G2[j, D] * psi[k] + G2[k, D] * psi[j] + G2[D, D] * psi[j] * psi[k]) / G1[D]);
  (-1)^(d - 1) * matdet(Hp)
}

closedform(A) =
{
  my(d = #A, S = A - A~);
  matdet(d^2 * matid(d) - S^2) / matdet(d^2 * matid(d) + S^2)
}

rhopoly(A) =
{
  my(d = #A, S = A - A~, chi = charpoly(S, 't), q);
  q = substpol(chi / 't, 't^2, 's);
  subst(q, 's, d^2) / abs(subst(q, 's, -d^2))
}

circ(d, D) = matrix(d, d, i, j, if (setsearch(Set(D), (j - i) % d), 1, 0));

isregtour(A) =
{
  my(d = #A);
  for (i = 1, d, if (A[i, i] != 0, return(0));
    if (sum(j = 1, d, A[i, j]) != (d - 1) / 2, return(0));
    for (j = 1, d, if (i != j && A[i, j] + A[j, i] != 1, return(0))));
  1
}

failures = 0;
check(name, A, expect) =
{
  my(d = #A, H, c, r, ok);
  if (!isregtour(A), error(Str(name, " is not a regular tournament")));
  H = hessH(A); c = closedform(A); r = rhopoly(A);
  ok = (H != 0) && (d / H == c) && (c == r^2) && (expect == 0 || r == expect);
  if (!ok, failures++);
  print(name, ": d/H = ", d / H, ", rho = ", r, if (ok, " ok", " MISMATCH"));
}

\\ R9': three copies of R3 arranged along R3
{
  r9p = matrix(9, 9, u, v, my(a = (u - 1) \ 3, b = (v - 1) \ 3, i = (u - 1) % 3, j = (v - 1) % 3);
               if (a == b, if ((j - i) % 3 == 1, 1, 0), if ((b - a) % 3 == 1, 1, 0)));
}

check("R3", circ(3, [1]), 2);
check("R5", circ(5, [1, 2]), 44/19);
check("R7", circ(7, [1, 2, 3]), 3032/1231);
check("Paley 7", circ(7, [1, 2, 4]), 64/27);
check("R9", circ(9, [1, 2, 3, 4]), 437584/171769);
check("R9'", r9p, 5488/2197);
check("circulant 9 {1,2,3,5}", circ(9, [1, 2, 3, 5]), 0);
\\ Wormald's 7-vertex example: cyclic triples abc, a'b'c', arcs aa', bb', cc', ad, bd, cd, completed
\\ to a regular tournament (vertex order a, b, c, a', b', c', d)
w7 = [0, 1, 0, 1, 0, 0, 1; 0, 0, 1, 0, 1, 0, 1; 1, 0, 0, 0, 0, 1, 1; 0, 1, 1, 0, 1, 0, 0; 1, 0, 1, 0, 0, 1, 0; 1, 1, 0, 1, 0, 0, 0; 0, 0, 0, 1, 1, 1, 0];
check("Wormald's 7-vertex example", w7, 3088/1287);
check("R11", circ(11, [1, 2, 3, 4, 5]), 108859552/41833969);

\\ (b) the cyclic product formula
default(realprecision, 80);
{
  for (h = 1, 30,
    my(d = 2 * h + 1, exact, prodv);
    exact = rhopoly(circ(d, vector(h, k, k)));
    prodv = prod(k = 0, h - 1, my(c = cotan(Pi * (2 * k + 1) / (2 * d))^2); (d^2 + c) / (d^2 - c));
    if (abs(exact - prodv) > 10^-60, failures++; print("R", d, " product formula MISMATCH"),
        print("R", d, ": rho = ", exact, " = product formula to 60 digits")));
}

print("limit cosh(1)/cos(1) = ", cosh(1) / cos(1));
print("failures ", failures);
quit(failures > 0);
