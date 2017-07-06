#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include <vector>
#include <utility>
#include <algorithm>
#include <inttypes.h>
#include <stdlib.h>

// Adjacency matrix of underlying graph.
int adj[888][888];
// Number of vertices in underlying graph.
int n;
// List of edges in underlying graph, written in both directions.  (For
// example, both (0,3) and (3,0) if there is an edge between 0 and 3.)
std::vector<std::pair<int, int> > arcs;

int lamb = -1;  // # common neighbours of two neighbours in SRG
int mu = -1;    // # common neighbours of two nonneighbours in SRG
int k = 0;      // degree of each vertex in SRG
uint64_t values[13];
bool setup_called = false;

void setup() {
  using std::max;
  if (setup_called) abort();

  lamb = -1;
  mu = -1;
  k = 0;

  for (int i = 0; i < n; i++) k += adj[0][i];
  for (int i = 0; i < n; i++) if (i && adj[0][i]) {
    lamb = 0;
    for (int j = 0; j < n; j++) if (adj[0][j] && adj[i][j]) lamb++;
  }
  for (int i = 0; i < n; i++) if (i && !adj[0][i]) {
    mu = 0;
    for (int j = 0; j < n; j++) if (adj[0][j] && adj[i][j]) mu++;
  }
  if (lamb == -1 || mu == -1) abort();

  values[0]  = (8*max(0, lamb-2) + 16 - 8*k) > 0;
  values[1]  = (8*max(0, lamb-1) + 4*(2-k)) > 0;
  values[2]  = (max(0, mu -2)*8 + 16 - 8*k) > 0;
  values[3]  = (8*max(mu-1, 0) + 8 - 4*k) > 0;
  values[4]  = (8*max(0, lamb-1) +8 - 4*k) > 0;
  values[5]  = 1;
  values[6]  = (8*max(mu-1, 0) + 8 - 4*k) > 0;
  values[7]  = 1;

  values[8]  = (8*max(0, lamb-1) + 8 - 4*k) > 0;
  values[9]  = (8*max(0, lamb-1) + 8 - 4*k) > 0;
  values[10] = (8*max(0, lamb-1) + 8 - 8*k + 2*k*k) > 0;
  values[11] = (8*max(0, mu-1) + 8 - 8*k + 2*k*k) > 0;
  values[12] = (k*k*(2-k)) > 0;

  setup_called = true;
}

// Returns A[i][j].
int aij(int i, int j) {
  int u = arcs[j].first;
  int v = arcs[j].second;
  int w = arcs[i].first;
  int x = arcs[i].second;

  int ndistinct = 2;
  ndistinct += (u != w) & (v != w);
  ndistinct += (u != x) & (v != x);

  int has = 0;
  if (ndistinct == 4) {  // common case
    has = values[4*!adj[u][w] + 2*!adj[v][w] + !adj[x][v]];
  } else if (ndistinct == 3) {
    if (v == w) {
    } else if (u == w) {
      if (adj[v][x]) has = values[8];
      else           has = 1;
    } else if (v == x) {
      if (adj[u][w]) has = values[9];
      else           has = 1;
    } else if (u == x) {
      if (adj[v][w]) has = values[10];
      else           has = values[11];
    }
  } else if (ndistinct == 2) {
    if (v == w) has = values[12];
    else has = 1;
  }

  return has;
}


// Compute M*v modulo P, returning the result.  Here, M is S^+(U^3(adj)).
// The entries of vv must be less than 2^64 / arcs.size() otherwise overflow
// modulo 2^64 will occur.  vv must be a row-major matrix with vstride columns.
template <uint64_t P, int vstride>
std::vector<uint64_t> hit(const std::vector<uint64_t> &vv) {
  if (!setup_called) abort();
  if (vv.size() != vstride * arcs.size()) abort();

  std::vector<uint64_t> ans(vv.size());

  static const int stride = 256;
  #pragma omp parallel for schedule(dynamic, 1)
  for (int ii = 0; ii < arcs.size(); ii += stride) {
    int imax = ii + stride;
    if (imax > arcs.size()) imax = arcs.size();
    for (int jj = 0; jj < arcs.size(); jj += stride) {
      int jmax = jj + stride;
      if (jmax > arcs.size()) jmax = arcs.size();
      for (int i = ii; i < imax; i++) {
        for (int j = jj; j < jmax; j++) {
          for (int k = 0; k < vstride; k++)
            ans[i*vstride+k] += aij(i, j) * vv[j*vstride+k];
        }
      }
    }
  }

  for (size_t i = 0; i < ans.size(); i++) ans[i] %= P;

  return ans;
}

// Returns x^k (mod P) provided P < 2^32 and x < 2^32.
template <uint64_t P>
uint64_t modpow(uint64_t x, uint64_t k) {
  if (!k) return 1;
  if (k == 1) return x;
  uint64_t y = modpow<P>((x*x)%P, k/2);
  if (k & 1) y = (x*y)%P;
  return y;
}

// Returns the modular inverse of x modulo P, provided P is prime, P < 2^32,
// and x < 2^32.
template <uint64_t P>
uint64_t modinv(uint64_t x) {
  return modpow<P>(x, P-2);
}

// Incrementally row-reduce vecs.  That is, subtract multiples of the first
// vecs.size()-1 rows of vecs from the last row of vecs so that the last row of
// vecs has a zero in every position where a previous row has its first
// non-zero entry.  Returns the index of the first nonzero in the last row of
// vecs.
//
// Requires that P, and all entries of vecs, are smaller than 2^32.
template <uint64_t P>
int reduce(std::vector<std::vector<uint64_t> > &vecs) {
  for (size_t i = 0; i < vecs.size() - 1; i++) {
    int piv = 0;
    while (!vecs[i][piv]) piv++;
    uint64_t mult = modinv<P>(vecs[i][piv]);
    if (vecs.back()[piv]) {
      uint64_t howmuch = P - (mult * vecs.back()[piv]) % P;
      for (size_t j = 0; j < vecs.back().size(); j++) {
        vecs.back()[j] = (vecs.back()[j] + howmuch * vecs[i][j]) % P;
      }
      if (vecs.back()[piv] != 0) abort();
    }
  }
  int idx = 0;
  while (idx < vecs.back().size() && !vecs.back()[idx]) idx++;
  return idx;
}

#endif
