#include "common.h"
#include <stdio.h>
using namespace std;

uint64_t tra2() {
  static const int stride = 256;
  uint64_t theans[arcs.size() / stride + 1] = {0};
  #pragma omp parallel for schedule(dynamic, 1)
  for (int ii = 0; ii < arcs.size(); ii += stride) {
    fprintf(stderr, "%i\n", ii);
    int imax = ii + stride;
    uint64_t ans = 0;
    if (imax > arcs.size()) imax = arcs.size();
    for (int jj = 0; jj < arcs.size(); jj += stride) {
      int jmax = jj + stride;
      if (jmax > arcs.size()) jmax = arcs.size();
      for (int i = ii; i < imax; i++) {
        for (int j = jj; j < jmax; j++) {
          ans += aij(i,j) && aij(j,i);
        }
      }
    }
    theans[ii/stride] = ans;
  }

  uint64_t ans = 0;
  for (int ii = 0; ii < arcs.size(); ii += stride) ans += theans[ii/stride];
  return ans;
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

  printf("%lli\n", tra2());
}
