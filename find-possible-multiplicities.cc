// The enumeration step at the end of this program assumes that a 'long' is a
// 64-bit integer, as it is on amd64 Linux.  It may break on other platforms.

#include <stdio.h>
#include <gmpxx.h>

#define FOR(i,n) for (int i=0;i<n;i++)

/*
Find all positive integer solutions to the following system:

x1+x2+x3+x4+x5+x6+x7+x8+2x9+3x10 = 98280
x1+15x2+125x3+127x4+68005x5-25x6-23x7-9x8+5426x9-799x10 = 98280
x1+225x2+15625x3+16129x4+4624680025x5+625x6+529x7+81x8+15185018x9+392663x10 = 6670853280
x1+3375x2+1953125x3+2048383x4+314501365100125x5-15625x6-12167x7-729x8+43716137114x9-192667111x10 = 318986389121400
x1+50625x2+244140625x3+260144641x4+21387665333634000625x5+390625x6+279841x7+6561x8+128961474307442x9+99596332307x10 = 21401273663621790120

where x4 = 1, x8 = 105, and x9 = 650.
*/

const char *csineqs[5][11] = {
{"1","1","1","1","1","1","1","1","2","3","98280"},
{"1","15","125","127","68005","-25","-23","-9","5426","-799","98280"},
{"1","225","15625","16129","4624680025","625","529","81","15185018","392663","6670853280"},
{"1","3375","1953125","2048383","314501365100125","-15625","-12167","-729","43716137114","-192667111","318986389121400"},
{"1","50625","244140625","260144641","21387665333634000625","390625","279841","6561","128961474307442","99596332307","21401273663621790120"}
};

mpz_class gcd(mpz_class a, mpz_class b) {
  return b ? gcd(b, a%b) : a;
}

int main() {
  if (sizeof(long) < 8) abort();

  mpz_class ineqs[5][8];
  FOR(i,5) {
    ineqs[i][0] = mpz_class(csineqs[i][0]);
    ineqs[i][1] = mpz_class(csineqs[i][1]);
    ineqs[i][2] = mpz_class(csineqs[i][2]);
    ineqs[i][3] = mpz_class(csineqs[i][3]);
    ineqs[i][4] = mpz_class(csineqs[i][5]);
    ineqs[i][5] = mpz_class(csineqs[i][6]);
    ineqs[i][6] = mpz_class(csineqs[i][7]);
    ineqs[i][7] = mpz_class(csineqs[i][10]);
    ineqs[i][7] -= mpz_class(csineqs[i][4]);
    ineqs[i][7] -= 105*mpz_class(csineqs[i][8]);
    ineqs[i][7] -= 650*mpz_class(csineqs[i][9]);
  }
  FOR(i, 5) {
    FOR(j,5) if (j != i) {
      mpz_class g = gcd(ineqs[i][i], ineqs[j][i]);
      mpz_class lc = ineqs[i][i] * ineqs[j][i] / g;
      mpz_class imul = lc / ineqs[i][i];
      mpz_class jmul = lc / ineqs[j][i];
      FOR(k, 8) ineqs[j][k] = ineqs[j][k] * jmul - ineqs[i][k] * imul;
    }
  }
  FOR(i,5) {
    mpz_class g = 0;
    FOR(k,8) g = gcd(g, ineqs[i][k]);
    FOR(k,8) ineqs[i][k] /= g;
  }
  const char *varname[8] = {"x1", "x2", "x3", "x4", "x6", "x7", "x8", "rhs"};
  FOR(i, 5) {
    FOR(j, 8) if (ineqs[i][j] != 0) gmp_printf("%+6Zi%s ", ineqs[i][j].get_mpz_t(), varname[j]);
    printf("\n");
  }

  // The system now has pattern
  // x....xx=x
  // .x...xx=x
  // ..x..xx=x
  // ...x.xx=x
  // ....xxx=x

  for (long x8 = 1; x8 < 98280; x8++) {
    FOR(i,5) {
      // ineqs[i][i] * x_i + ineqs[i][5] * x_7 = rhs, so rhs must be an integer
      // linear combination of ineqs[i][i] and ineqs[i][5].
      mpz_class rhs = ineqs[i][7] - ineqs[i][6] * x8;
      mpz_class g = gcd(ineqs[i][i], ineqs[i][5]);
      if (rhs%g != 0) goto skip;
    }
    for (long x7 = 1; x7 < 98280; x7++) {
      mpz_class x[5];
      FOR(i,5) {
        x[i] = ineqs[i][7] - x8 * ineqs[i][6] - x7 * ineqs[i][5];
        if (x[i] % ineqs[i][i] != 0) goto innerskip;
        x[i] /= ineqs[i][i];
        if (x[i] < 0) goto innerskip;
      }
      FOR(i,5) printf("%li ", x[i].get_si());
      printf("%li %li\n", x7, x8);
      innerskip:;
    }
    skip:;
  }
}
