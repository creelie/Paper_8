/* regtour.c: all labelled regular tournaments on d vertices (d odd <= 11) with N+(0) = {1..m},
 * m = (d-1)/2; prints each distinct characteristic polynomial of S = A - A^T once, with the number of
 * tournaments having it and one adjacency matrix (rows as 0/1 strings). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int d, m, A[16][16], outd[16], ind[16], pi[128], pj[128], np;
typedef struct { int64_t c[16]; long cnt; char mat[16][17]; } Rec;
static Rec *recs; static int nrec, caprec;
static void charpoly(int64_t *c) { /* Faddeev-LeVerrier in exact integers */
  int64_t S[16][16], Mk[16][16], T[16][16];
  for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) { S[i][j] = A[i][j] - A[j][i]; Mk[i][j] = 0; }
  c[0] = 1;
  for (int k = 1; k <= d; k++) {
    for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) {
      int64_t s = 0; for (int l = 0; l < d; l++) s += S[i][l] * Mk[l][j];
      T[i][j] = s + (i == j ? c[k - 1] : 0);
    }
    memcpy(Mk, T, sizeof Mk);
    int64_t tr = 0;
    for (int i = 0; i < d; i++) for (int l = 0; l < d; l++) tr += S[i][l] * Mk[l][i];
    if (tr % k) { fprintf(stderr, "inexact\n"); exit(1); }
    c[k] = -tr / k;
  }
}
static long total;
static void record(void) {
  int64_t c[16]; charpoly(c); total++;
  for (int r = 0; r < nrec; r++) if (!memcmp(recs[r].c, c, sizeof c)) { recs[r].cnt++; return; }
  if (nrec == caprec) { caprec = caprec ? 2 * caprec : 64; recs = realloc(recs, caprec * sizeof *recs); }
  Rec *R = &recs[nrec++]; memcpy(R->c, c, sizeof c); R->cnt = 1;
  for (int i = 0; i < d; i++) { for (int j = 0; j < d; j++) R->mat[i][j] = '0' + A[i][j]; R->mat[i][d] = 0; }
}
static void rec(int k) {
  if (k == np) { record(); return; }
  int i = pi[k], j = pj[k];
  /* i -> j */
  if (outd[i] < m && ind[j] < m) { A[i][j] = 1; outd[i]++; ind[j]++; rec(k + 1); A[i][j] = 0; outd[i]--; ind[j]--; }
  if (outd[j] < m && ind[i] < m) { A[j][i] = 1; outd[j]++; ind[i]++; rec(k + 1); A[j][i] = 0; outd[j]--; ind[i]--; }
}
int main(int argc, char **argv) {
  if (argc != 2 || (d = atoi(argv[1])) < 3 || d > 11 || d % 2 == 0) {
    fprintf(stderr, "usage: regtour d   (d odd, 3 <= d <= 11)\n");
    return 1;
  }
  m = (d - 1) / 2;
  for (int j = 1; j < d; j++) { if (j <= m) { A[0][j] = 1; outd[0]++; ind[j]++; } else { A[j][0] = 1; outd[j]++; ind[0]++; } }
  for (int i = 1; i < d; i++) for (int j = i + 1; j < d; j++) { pi[np] = i; pj[np] = j; np++; }
  rec(0);
  printf("total %ld distinct %d\n", total, nrec);
  for (int r = 0; r < nrec; r++) {
    printf("poly"); for (int k = 0; k <= d; k++) printf(" %lld", (long long)recs[r].c[k]);
    printf(" count %ld mat", recs[r].cnt); for (int i = 0; i < d; i++) printf(" %s", recs[r].mat[i]); printf("\n");
  }
  return 0;
}
