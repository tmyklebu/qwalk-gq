#include "common.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
using namespace std;

#define FOR(i,n) for (int i=0;i<n;i++)

// Compute the minimal polynomial of S^+(U^3(adj)) on the Krylov space of a
// pseudorandom vector coming from rand() modulo P.  Prints the following line
// to stdout:
//
//   [P] [degree+1] [a_0] [a_1] ... [a_degree]
template <uint64_t P>
void doit() {
  vector<uint64_t> v(arcs.size());
  // Generate pseudorandom vector.
  for (size_t i = 0; i < arcs.size(); i++)
    v[i] = rand() % P;

  vector<vector<uint64_t> > vecs;
  while (1) {
    // Do row-reduction with a right-appended identity matrix so we can recover
    // the minimal polynomial from the first zero row that occurs.
    FOR(i,vecs.size()) vecs[i].push_back(0);
    vecs.push_back(v);
    while (vecs.back().size() < vecs[0].size()-1) vecs.back().push_back(0);
    vecs.back().push_back(1);
    int pivot = reduce<P>(vecs);
    fprintf(stderr, "%i vectors; pivot is %i\n", (int)vecs.size(), pivot);
    if (pivot >= n) break;
    v = hit<P,1>(v);
  }

  vector<uint64_t> poly;
  for (int i = arcs.size(); i < vecs.back().size(); i++)
    poly.push_back(vecs.back()[i]);

  printf("%lli %i", P, (int)poly.size());
  FOR(i, poly.size()) printf(" %lli", poly[i]);
  printf("\n"); fflush(stdout);
}

int main() {
  scanf("%i", &n);
  while (1) {
    int a,b;
    if (2 != scanf("%i %i", &a, &b)) break;
    adj[a][b] = adj[b][a] = 1;
    arcs.push_back(make_pair(a,b));
  }
  setup();

  doit<1999999003>();
  doit<1999999013>();
  doit<1999999049>();
  doit<1999999061>();
  doit<1999999081>();
  doit<1999999087>();
  doit<1999999093>();
  doit<1999999097>();
  doit<1999999117>();
  doit<1999999121>();
  doit<1999999151>();
  doit<1999999171>();
}
