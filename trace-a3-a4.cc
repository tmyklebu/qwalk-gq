#include "common.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <bitset>
#include <math.h>
#include <gmpxx.h>
#include <mutex>
using namespace std;

#define FOR(i,n) for (int i=0;i<n;i++)

int a2ij(int i, int j) {  // A^2[i][j]
  int ans = 0;
  FOR(k, arcs.size()) ans += aij(i,k) && aij(k,j);
  return ans;
}

static const int64_t p1 = 123456789012353ll;
static const int64_t p2 = 123456789012373ll;
static const int64_t p3 = 123456789012419ll;
static const int64_t p4 = 123456789012421ll;

static const int stride = 256;

// iblk[i][j] == A[i][j]
bitset<98304> *iblk = new bitset<98304>[100000];
// jblkt[i][j] == A[j][i]
bitset<98304> *jblkt = new bitset<98304>[100000];

void calc_blks() {
  #pragma omp parallel for schedule(dynamic, 256)
  FOR(i, arcs.size()) {
    if (!(i % 1024)) printf("iblk %i\n", i);
    FOR(k, arcs.size()) iblk[i][k] = aij(i, k);
  }
  #pragma omp parallel for schedule(dynamic, 256)
  FOR(j, arcs.size()) {
    if (!(j % 1024)) printf("jblk %i\n", j);
    FOR(k, arcs.size()) jblkt[j][k] = aij(k, j);
  }
}

// If imin == jmin:
//   tra3 += sum(imin <= i, j < imax) (A^2)_{ij} A_{ji}
//   tra4 += sum(imin <= i, j < imax) (A^2)_{ij} (A^2)_{ji}
// Otherwise:
//   tra3 += sum(i=imin..imax-1) sum(j=jmin..jmax-1)
//               (A^2)_{ij} A_{ji} + (A^2)_{ji} A_{ij}
//   tra4 += sum(i=imin..imax-1) sum(j=jmin..jmax-1)
//               (A^2)_{ij} (A^2)_{ji} + (A^2)_{ji} (A^2)_{ij}
void tra34_block(int imin, int imax, int jmin, int jmax,
    uint64_t &tra3, uint64_t &tra4) {
  // blk2[i-imin][j-jmin] = sum_k A[i][k] A[k][j] = (A^2)_{ij}
  int blk2[stride][stride];
  // klb2[j-jmin][i-imin] = sum_k A[i][k] A[k][j] = (A^2)_{ij}
  int klb2[stride][stride];

  if (imin == jmin) {
    FOR(i,stride) FOR(j,stride)
      blk2[i][j] = (iblk[imin+i] & jblkt[jmin+j]).count();
  
    tra3 = tra4 = 0;
    FOR(i,stride) FOR(j,stride) {
      tra3 += iblk[jmin+j][imin+i] * blk2[i][j];
      tra4 += (uint64_t)blk2[i][j] * (uint64_t)blk2[j][i];
    }
  } else {
    FOR(i,stride) FOR(j,stride) {
      blk2[i][j] = (iblk[imin+i] & jblkt[jmin+j]).count();
      klb2[i][j] = (iblk[jmin+j] & jblkt[imin+i]).count();
    }
  
    tra3 = tra4 = 0;
    FOR(i,stride) FOR(j,stride) {
      tra3 += iblk[jmin+j][imin+i] * blk2[i][j];
      tra3 += iblk[imin+i][jmin+j] * klb2[i][j];
      tra4 += 2 * (uint64_t)blk2[i][j] * (uint64_t)klb2[i][j];
    }
  }
}

// Compute tra3 = tr(A^3), tra4_p1 = tr(A^4) % p1, ..., tra4_p4 = tr(A^4) % p4.
// Uses the identities
//   tr(A^3) = <A^2, A^T) = sum_{ij} (A^2)_{ij} A_{ji}
//   tr(A^4) = <A^2, (A^2)^T) = sum_{ij} (A^2)_{ij} (A^2)_{ji}
// computing both in a tiled fashion.
void tra34(uint64_t &tra3, uint64_t &tra4_p1, uint64_t &tra4_p2, uint64_t &tra4_p3, uint64_t &tra4_p4) {
  uint64_t a3[384][384], a4[384][384];
  memset(a3, 0, sizeof(a3));
  memset(a4, 0, sizeof(a4));
  calc_blks();
  printf("done calc_blks\n");
  #pragma omp parallel for schedule(dynamic, 1)
  for (int ii = 0; ii < (int)arcs.size(); ii += stride) {
    int imax = min((int)arcs.size(), ii+stride);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int jj = 0; jj <= ii; jj += stride) {
      printf("%i %i\n", ii, jj);
      int jmax = min((int)arcs.size(), jj+stride);
      uint64_t la3=0, la4=0;
      tra34_block(ii, imax, jj, jmax, la3, la4);
      a3[ii/stride][jj/stride] = la3;
      a4[ii/stride][jj/stride] = la4;
    }
  }
  tra3 = tra4_p1 = tra4_p2 = tra4_p3 = tra4_p4 = 0;
  // 123456789054321 * 384^2 ~= 1.820e19 < 2^64 ~= 1.844e19.
  FOR(i,384) FOR(j,384) {
    tra3 += a3[i][j];
    tra4_p1 += a4[i][j] % p1;
    tra4_p2 += a4[i][j] % p2;
    tra4_p3 += a4[i][j] % p3;
    tra4_p4 += a4[i][j] % p4;
  }
  tra4_p1 %= p1;
  tra4_p2 %= p2;
  tra4_p3 %= p3;
  tra4_p4 %= p4;
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
  uint64_t tra3;
  uint64_t tra4_p1, tra4_p2, tra4_p3, tra4_p4;
  tra34(tra3, tra4_p1, tra4_p2, tra4_p3, tra4_p4);
  printf("tra3 = %lli\n", tra3);
  printf("tra4 = %lli mod %lli\n", tra4_p1, p1);
  printf("tra4 = %lli mod %lli\n", tra4_p2, p2);
  printf("tra4 = %lli mod %lli\n", tra4_p3, p3);
  printf("tra4 = %lli mod %lli\n", tra4_p4, p4);
}
