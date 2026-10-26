"""Exact checks of the analytic Hessian matrix and its quantitative bounds.

Differentiate Wormald's determinant definition, compare every matrix entry
with E^T (d^2 I + S^2)(d^2 I - S^2)^(-1) E, and check the two rational
positive-semidefinite bounds by exact Schur elimination. Finite checks
complement the universal written argument; they are not its proof.
"""
from fractions import Fraction as Q
from pathlib import Path
import random
import wormald_hessian as original

def mul(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),Q(0))
             for j in range(len(b[0]))] for i in range(len(a))]

def inverse(a):
    n=len(a);c=[list(map(Q,row))+[Q(i==j) for j in range(n)] for i,row in enumerate(a)]
    for j in range(n):
        k=next(k for k in range(j,n) if c[k][j]);c[j],c[k]=c[k],c[j]
        q=c[j][j];c[j]=[v/q for v in c[j]]
        for i in range(n):
            if i!=j:
                q=c[i][j];c[i]=[u-q*v for u,v in zip(c[i],c[j])]
    return [row[n:] for row in c]

def psd(a):
    a=[list(map(Q,r)) for r in a];n=len(a)
    assert all(a[i][j]==a[j][i] for i in range(n) for j in range(n))
    for k in range(n):
        q=a[k][k]
        if q<0:return False
        if q==0:
            if any(a[k][j] for j in range(k+1,n)):return False
            continue
        for i in range(k+1,n):
            for j in range(k+1,n):a[i][j]-=a[i][k]*a[k][j]/q
    return True

def implicit_matrix(a):
    d=len(a);z=Q(2,d)
    cols=[[Q(i==j)-z*(Q(i==j)/2+a[i][j]) for i in range(d)] for j in range(d)]
    dc=lambda j:[-Q(i==j)/2-a[i][j] for i in range(d)]
    detcols=lambda c:original.det([[c[j][i] for j in range(d)] for i in range(d)])
    assert detcols(cols)==0
    g=[]
    for j in range(d):
        c=cols[:];c[j]=dc(j);g.append(z*detcols(c))
    g2=[[Q(0) for _ in range(d)] for _ in range(d)]
    for j in range(d):
        g2[j][j]=g[j]
        for k in range(j+1,d):
            c=cols[:];c[j]=dc(j);c[k]=dc(k);g2[j][k]=g2[k][j]=z*z*detcols(c)
    assert g[-1]!=0 and all(v==g[-1] for v in g)
    return [[(g2[j][k]-g2[j][-1]-g2[-1][k]+g2[-1][-1])/g[-1]
             for k in range(d-1)] for j in range(d-1)]

def check(a):
    assert original.is_regular_tournament(a)
    d=len(a);n=d-1;s=[[Q(a[i][j]-a[j][i]) for j in range(d)] for i in range(d)]
    s2=mul(s,s)
    plus=[[Q(d*d*(i==j))+s2[i][j] for j in range(d)] for i in range(d)]
    minus=[[Q(d*d*(i==j))-s2[i][j] for j in range(d)] for i in range(d)]
    f=mul(plus,inverse(minus))
    proposed=[[f[j][k]-f[j][-1]-f[-1][k]+f[-1][-1] for k in range(n)] for j in range(n)]
    actual=implicit_matrix(a);assert actual==proposed
    alpha=Q(2*d-1,2*d*d-2*d+1)
    gram=[[Q((i==j)+1) for j in range(n)] for i in range(n)]
    assert psd([[actual[i][j]-alpha*gram[i][j] for j in range(n)] for i in range(n)])
    assert psd([[gram[i][j]-actual[i][j] for j in range(n)] for i in range(n)])
    det=original.det(actual);assert Q(d)*alpha**n<=det<=d
    return n*n

def main():
    cases=[a for _,a,_ in original.TABLE];root=Path(__file__).resolve().parents[2]
    for d in (3,5,7,9):cases.extend(a for _,a in original.classes(root/'checks'/f'data/classes{d}.txt'))
    rng=random.Random(20261014)
    for d,count in ((11,12),(13,4)):
        cases.extend(original.random_regular(d,400,rng) for _ in range(count))
    entries=sum(check(a) for a in cases)
    print(f'EXACT MATRIX AUDIT: {len(cases)} tournaments, {entries} entries; both Loewner bounds and determinant bounds passed.')

if __name__=='__main__':main()
