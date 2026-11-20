#ifdef MIINT_HAS_HMMER
// miint::hmmer, the layer between HMMER and read_hmm: Easel's exception handler,
// HMMER's one-time setup, and ParseProfile, which reads one model's header.
//
// The models are HMMER's own test files at the pinned fork (ext/hmmer), with
// trailing whitespace removed: testsuite/M1.hmm, XYPPX.hmm and 3box.hmm, and
// the two HMMER2 models written out by testsuite/i15-hmmconvert.pl. Expected
// header values are as written in the files. The HMMER2 statistics, which HMMER
// computes on read, are those the fork's hmmconvert prints for the same models.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <new>
#include <string>
#include <thread>
#include <vector>

#include "HmmerProfile.hpp"

extern "C" {
#include "hmmer.h"
}

using Catch::Matchers::WithinAbs;
using miint::hmmer::Initialize;
using miint::hmmer::ParseProfile;
using miint::hmmer::ProfileHeader;
using miint::hmmer::ProfileParseError;
using miint::hmmer::TakeExceptionMessage;
using miint::hmmer::detail::ThrowIfReadFailed;

namespace {

// HMMER3/e, amino, one node. No ACC, DESC, MAXL or cutoffs; CKSUM 0.
const char *const M1_HMM = R"hmm(HMMER3/e [3.0 | March 2010]
NAME  M1
LENG  1
ALPH  amino
RF    no
CONS  yes
CS    no
MAP   yes
DATE  Thu Jun 16 11:49:23 2011
NSEQ  2
EFFN  2.000000
CKSUM 0
STATS LOCAL MSV       -4.9046  1.46065
STATS LOCAL VITERBI   -5.1999  1.46065
STATS LOCAL FORWARD   -0.0517  1.46065
HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   0.35624  4.75100  4.29053  4.19828  4.88306  3.48560  5.11004  4.21322  4.22519  4.00099  5.04236  4.13232  4.23379  4.52899  4.38323  3.09073  3.42594  3.72730  6.17838  5.12044
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.01467  4.62483  5.34718  0.61958  0.77255  0.00000        *
      1   0.32372  4.76536  4.42980  4.32857  5.00499  3.55951  5.22620  4.27004  4.37081  4.10495  5.08789  4.22499  4.36948  4.63911  4.51684  3.12947  3.46009  3.76842  6.33337  5.25783      1 A - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.00990  4.62006        *  0.61958  0.77255  0.00000        *
//
)hmm";

// HMMER3/e, amino. ACC, DESC and all three cutoffs.
const char *const XYPPX_HMM = R"hmm(HMMER3/e [3.0 | March 2010]
NAME  XYPPX
ACC   PF02162.8
DESC  XYPPX repeat
LENG  5
ALPH  amino
RF    no
CONS  yes
CS    no
MAP   yes
DATE  Thu Jun 16 11:48:52 2011
NSEQ  7
EFFN  7.000000
CKSUM 149732247
GA    14.30 0.70
TC    14.40 0.70
NC    14.20 4.00
STATS LOCAL MSV       -5.1112  0.78852
STATS LOCAL VITERBI   -4.7526  0.78852
STATS LOCAL FORWARD   -2.5557  0.78852
HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   1.78915  6.00913  5.23249  5.11791  5.20653  2.22796  5.88594  5.77670  5.17236  5.27663  6.29854  5.19187  0.83219  2.40628  5.36874  4.33373  4.68171  5.20184  6.65125  2.09768
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.00529  5.64135  6.36369  0.61958  0.77255  0.00000        *
      1   0.91707  5.48468  5.65823  5.70004  6.37399  0.67196  6.42141  5.91783  5.81567  5.60504  6.35161  5.14615  5.14089  5.87323  5.80244  3.79498  4.17581  4.94514  7.64062  6.69238      1 G - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.00529  5.64135  6.36369  0.61958  0.77255  0.48576  0.95510
      2   1.29214  5.93727  5.64450  5.57968  3.89343  4.88067  5.29214  5.20374  5.53411  4.64869  5.82218  5.30271  5.60238  5.64395  5.57005  4.42103  4.74169  4.86331  5.54271  0.49105      2 Y - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.00529  5.64135  6.36369  0.61958  0.77255  0.48576  0.95510
      3   5.94672  7.11992  6.55030  6.62396  7.25663  5.73661  7.31745  7.48876  6.82943  6.72457  7.95279  6.80729  0.02290  7.12774  6.70449  6.20885  6.45690  6.95162  7.77702  7.45295      3 P - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.00529  5.64135  6.36369  0.61958  0.77255  0.48576  0.95510
      4   5.94672  7.11992  6.55030  6.62396  7.25663  5.73661  7.31745  7.48876  6.82943  6.72457  7.95279  6.80729  0.02290  7.12774  6.70449  6.20885  6.45690  6.95162  7.77702  7.45295      4 P - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.00529  5.64135  6.36369  0.61958  0.77255  0.48576  0.95510
      5   1.83526  5.65654  4.22930  4.02694  5.84734  4.24624  5.28601  5.33049  4.08818  4.90398  5.71607  4.35468  1.48546  0.81000  4.42096  3.73357  4.07553  4.71157  7.04231  5.78514      5 q - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.00356  5.63962        *  0.61958  0.77255  0.00000        *
//
)hmm";

// HMMER3/f, DNA, with MAXL.
const char *const BOX3_HMM = R"hmm(HMMER3/f [3.3 | Nov 2019]
NAME  3box
LENG  20
MAXL  75
ALPH  DNA
RF    yes
MM    no
CONS  yes
CS    yes
MAP   yes
DATE  Thu Mar 19 09:44:43 2020
NSEQ  22
EFFN  22.000000
CKSUM 1874497603
STATS LOCAL MSV       -6.0857  0.74408
STATS LOCAL VITERBI   -6.5394  0.74408
STATS LOCAL FORWARD   -3.2146  0.74408
HMM          A        C        G        T
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   0.85352  3.67699  1.53220  1.10043
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  0.00000        *
      1   0.00850  6.21185  5.49634  6.05242      1 A x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
      2   0.00850  6.21185  5.49634  6.05242      2 A x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
      3   0.07174  5.97461  2.75317  5.82342      3 A x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
      4   4.57582  4.38815  1.28255  0.35675      4 t x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
      5   5.04480  6.00326  0.01203  5.79333      5 G x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
      6   0.00850  6.21185  5.49634  6.05242      6 A x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
      7   0.16698  4.64734  5.03178  1.98292      7 a x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
      8   0.27480  2.21511  3.99743  2.18242      8 a x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
      9   0.81098  4.17232  4.48858  0.63691      9 t x - .
          1.38629  1.38629  1.38629  1.38629
          0.27817  5.48894  1.43256  1.46634  0.26236  1.09861  0.40547
     10   1.24764  1.81454  1.93201  0.90375     10 t . - .
          1.38629  1.38629  1.38629  1.38629
          0.01086  5.22163  5.22163  1.46634  0.26236  0.15006  1.97081
     11   2.62887  4.93134  0.09260  4.70307     11 G x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
     12   5.28311  4.88969  5.40686  0.01723     12 T x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
     13   2.55388  3.81030  2.68966  0.18371     13 t x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
     14   5.28311  4.88969  5.40686  0.01723     14 T x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
     15   5.04480  6.00326  0.01203  5.79333     15 G x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
     16   0.08235  4.92354  5.26649  2.70886     16 A x - .
          2.00737  2.00737  2.00737  0.51588
          0.09989  2.39766  5.48894  0.57943  0.82147  1.09861  0.40547
     17   2.48998  3.46889  1.86103  0.31413     19 t x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
     18   5.28311  4.88969  5.40686  0.01723     20 T x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
     19   0.00850  6.21185  5.49634  6.05242     21 A x - .
          1.38629  1.38629  1.38629  1.38629
          0.00830  5.48894  5.48894  1.46634  0.26236  1.09861  0.40547
     20   2.00283  3.72053  0.42449  1.67812     22 g x - .
          1.38629  1.38629  1.38629  1.38629
          0.00416  5.48480        *  1.46634  0.26236  0.00000        *
//
)hmm";

// HMMER2.0, nucleic (DNA to HMMER 3), from i15-hmmconvert.pl's h2_nt_file.
const char *const H2_NT_HMM = R"hmm(HMMER2.0  [2.3.2]
NAME  ecori
LENG  6
ALPH  Nucleic
RF    no
CS    no
MAP   yes
COM   hmmbuild ecori.hmm testsuite/ecori.sto
NSEQ  2
DATE  Tue May 18 13:30:52 2010
CKSUM 3108
XT      -9967     -1  -1000  -1000  -9967     -1  -9967     -1
NULT      -1  -9967
NULE       0      0      0      0
HMM        A      C      G      T
         m->m   m->i   m->d   i->m   i->i   d->m   d->d   b->m   m->e
         -415      *  -2000
     1   -585   -585   1000   -585     1
     -      0      0      0      0
     -    -21  -6672  -7714   -894  -1115   -701  -1378   -415      *
     2   1000   -585   -585   -585     2
     -      0      0      0      0
     -    -21  -6672  -7714   -894  -1115   -701  -1378      *      *
     3   1000   -585   -585   -585     3
     -      0      0      0      0
     -    -21  -6672  -7714   -894  -1115   -701  -1378      *      *
     4   -585   -585   -585   1000     4
     -      0      0      0      0
     -    -21  -6672  -7714   -894  -1115   -701  -1378      *      *
     5   -585   -585   -585   1000     5
     -      0      0      0      0
     -    -21  -6672  -7714   -894  -1115   -701  -1378      *      *
     6   -585   1000   -585   -585     6
     -      *      *      *      *
     -      *      *      *      *      *      *      *      *      0
//
)hmm";

// HMMER2.0, amino, from i15-hmmconvert.pl's h2_aa_file.
const char *const H2_AA_HMM = R"hmm(HMMER2.0  [2.3.2]
NAME  test
LENG  20
ALPH  Amino
RF    yes
CS    no
MAP   yes
COM   hmmbuild 20aa.hmm testsuite/20aa.sto
NSEQ  10
DATE  Tue May 18 13:38:05 2010
CKSUM 9020
XT      -8455     -4  -1000  -1000  -8455     -4  -8455     -4
NULT      -4  -8455
NULE     595  -1558     85    338   -294    453  -1158    197    249    902  -1085   -142    -21   -313     45    531    201    384  -1998   -644
HMM        A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y
         m->m   m->i   m->d   i->m   i->i   d->m   d->d   b->m   m->e
         -585      *  -1585
     1   2806   -444  -1270  -1265  -1732   -744  -1182  -1034  -1212  -1545  -1010   -869  -1380  -1102  -1325   -156   -259   -640  -2039  -1721     1
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378   -585      *
     2   -608   5196  -2368  -2362  -1635  -1253  -1702   -987  -2034  -1540  -1191  -1745  -1831  -1953  -1922   -934   -927   -770  -1855  -1692     2
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
     3  -1130  -2063   3441    243  -2502  -1157   -662  -2514   -843  -2578  -2029    -96  -1658   -421  -1434   -970  -1239  -2138  -2388  -1915     3
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
     4   -988  -1921    379   3052  -2309  -1196   -495  -2108   -350  -2194  -1612   -103  -1621   -210   -770   -862  -1042  -1799  -2207  -1742     4
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
     5  -1602  -1340  -2512  -2455   3740  -2315   -379   -478  -2225   -175   -190  -1874  -2503  -1781  -2079  -1892  -1647   -683    210   1257     5
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
     6   -759  -1104  -1314  -1488  -2461   3330  -1563  -2503  -1721  -2645  -2082  -1325  -1776  -1591  -1791   -959  -1094  -1922  -2199  -2280     6
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
     7  -1192  -1547   -812   -781   -520  -1557   4586  -1881   -457  -1789  -1357   -866  -1913   -669   -600  -1231  -1258  -1678   -918    -90     7
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
     8  -1016   -805  -2637  -2367   -591  -2470  -1853   3039  -2040    329    346  -2180  -2589  -1943  -2070  -1879  -1051   1184  -1698  -1261     8
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
     9  -1028  -1686   -800   -442  -2171  -1521   -263  -1873   3103  -1859  -1251   -547  -1750     74    599  -1004   -972  -1620  -1812  -1528     9
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
    10  -1369  -1199  -2542  -2295   -129  -2453  -1603    488  -1867   2607    742  -2168  -2539  -1724  -1835  -1992  -1384    130  -1284   -897    10
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
    11   -965   -982  -2125  -1844   -276  -2022  -1302    373  -1390    610   4141  -1701  -2229  -1356  -1425  -1444  -1028    126  -1279   -880    11
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
    12   -637  -1349     -5   -148  -1773  -1078   -624  -1965   -499  -2112  -1514   3456  -1574   -435   -810   -661   -818  -1560  -1900  -1327    12
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
    13   -934  -1220  -1352  -1451  -2187  -1314  -1456  -2160  -1488  -2268  -1827  -1359   3760  -1471  -1571  -1129  -1216  -1780  -2050  -2023    13
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
    14   -862  -1548   -321   -247  -1672  -1336   -419  -1721    130  -1660  -1161   -427  -1678   3540    -59   -855   -904  -1495  -1728  -1237    14
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
    15  -1181  -1614  -1431   -903  -2066  -1632   -446  -1944    559  -1867  -1326   -899  -1873   -156   3432  -1225  -1148  -1717  -1746  -1529    15
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
    16    116   -498   -840   -867  -1751   -671   -949  -1625   -916  -1896  -1218   -622  -1304   -824  -1097   2786   -220  -1038  -1992  -1530    16
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
    17    -52   -555  -1173  -1112  -1571   -872  -1041   -862   -950  -1353   -866   -843  -1455   -937  -1084   -285   3069   -557  -1898  -1539    17
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
    18   -700   -665  -2337  -2135   -873  -1965  -1727   1107  -1863   -105     22  -1851  -2286  -1794  -1920  -1368   -823   2927  -1816  -1362    18
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
    19  -1926  -1577  -2416  -2370    192  -2070   -795  -1602  -1882  -1347  -1210  -2018  -2383  -1854  -1777  -2150  -1992  -1643   5696    560    19
     x   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249
     -    -33  -6045  -7087   -894  -1115   -701  -1378      *      *
    20  -1667  -1473  -2001  -1964   1181  -2210    -25  -1224  -1658  -1018   -871  -1505  -2415  -1423  -1618  -1747  -1692  -1248    501   4149    20
     x      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *
     -      *      *      *      *      *      *      *      *      0
//
)hmm";

// The precision hmmconvert prints the statistics with: mu and tau to 4
// decimals, lambda to 5.
constexpr double MU_PRINTED = 0.00005;
constexpr double LAMBDA_PRINTED = 0.000005;

// HMMER reads a header number with atof() and stores it in a float: the
// expected value is the same conversion of the number as written.
float F(double value) {
	return static_cast<float>(value);
}

// <text> with its one occurrence of <from> replaced by <to>.
std::string Replace(std::string text, const std::string &from, const std::string &to) {
	auto at = text.find(from);
	REQUIRE(at != std::string::npos);
	REQUIRE(text.find(from, at + 1) == std::string::npos);
	return text.replace(at, from.size(), to);
}

// The first <n> lines of <text>, each with its newline.
std::string FirstLines(const std::string &text, int n) {
	size_t end = 0;
	for (int i = 0; i < n; i++) {
		end = text.find('\n', end);
		REQUIRE(end != std::string::npos);
		end++;
	}
	return text.substr(0, end);
}

// Raises an Easel exception on this thread, and returns the status the raising
// call returned: p7_hmm_CreateBody refuses a model length outside
// 1..p7_MAXM.
int RaiseModelLengthException(int length) {
	ESL_ALPHABET *abc = esl_alphabet_Create(eslAMINO);
	P7_HMM *hmm = p7_hmm_CreateShell();
	int status = p7_hmm_CreateBody(hmm, length, abc);
	p7_hmm_Destroy(hmm);
	esl_alphabet_Destroy(abc);
	return status;
}

void ExpectParseError(const std::function<void()> &call, const std::string &message, int64_t line) {
	try {
		call();
	} catch (const ProfileParseError &e) {
		CHECK(std::string(e.what()) == message);
		CHECK(e.Line() == line);
		return;
	}
	FAIL("no ProfileParseError was thrown");
}

void ExpectRecordError(const std::string &record, const std::string &message, int64_t line) {
	ExpectParseError([&] { ParseProfile(record); }, message, line);
}

void CheckNoCutoffs(const ProfileHeader &h) {
	CHECK_FALSE(h.gathering_threshold_sequence.has_value());
	CHECK_FALSE(h.gathering_threshold_domain.has_value());
	CHECK_FALSE(h.trusted_cutoff_sequence.has_value());
	CHECK_FALSE(h.trusted_cutoff_domain.has_value());
	CHECK_FALSE(h.noise_cutoff_sequence.has_value());
	CHECK_FALSE(h.noise_cutoff_domain.has_value());
}

} // namespace

TEST_CASE("HMMER profile - an Easel exception is recorded for the thread, not fatal", "[hmmer]") {
	Initialize();
	TakeExceptionMessage(); // anything an earlier test left
	// Without a handler that returns, this call aborts the process.
	REQUIRE(RaiseModelLengthException(0) == eslEINVAL);
	CHECK(TakeExceptionMessage() == "model length M=0 out of range 1..100000");
	// Taking the message clears it.
	CHECK(TakeExceptionMessage().empty());
}

#if !defined(__EMSCRIPTEN__) || defined(__EMSCRIPTEN_PTHREADS__)
TEST_CASE("HMMER profile - each thread sees only its own exception", "[hmmer]") {
	Initialize();
	TakeExceptionMessage();
	REQUIRE(RaiseModelLengthException(0) == eslEINVAL);

	std::string seen_before, seen_after;
	std::thread other([&] {
		seen_before = TakeExceptionMessage();
		RaiseModelLengthException(100001);
		seen_after = TakeExceptionMessage();
	});
	other.join();

	CHECK(seen_before.empty());
	CHECK(seen_after == "model length M=100001 out of range 1..100000");
	CHECK(TakeExceptionMessage() == "model length M=0 out of range 1..100000");
}
#endif

TEST_CASE("HMMER profile - Initialize fills HMMER's logsum table", "[hmmer]") {
	Initialize();
	// log(e^0 + e^0) = log 2. An unfilled table holds zeros, so p7_FLogsum
	// would return 0.
	CHECK_THAT(p7_FLogsum(0.0f, 0.0f), WithinAbs(std::log(2.0), 1e-6));
}

TEST_CASE("HMMER profile - M1: the header, and the fields it does not have", "[hmmer]") {
	Initialize();
	ProfileHeader h = ParseProfile(M1_HMM);
	CHECK(h.name == "M1");
	CHECK(h.length == 1);
	CHECK(h.alphabet == "amino");
	CHECK(h.format == "HMMER3/e");
	CHECK(h.sequence_count == 2);
	CHECK(h.effective_sequence_count == F(2.000000));
	// CKSUM 0 is a checksum of 0, not a missing one.
	CHECK(h.checksum == 0u);
	CHECK(h.consensus == "A");
	CHECK(h.msv_mu == F(-4.9046));
	CHECK(h.msv_lambda == F(1.46065));
	CHECK(h.viterbi_mu == F(-5.1999));
	CHECK(h.viterbi_lambda == F(1.46065));
	CHECK(h.forward_tau == F(-0.0517));
	CHECK(h.forward_lambda == F(1.46065));

	CHECK_FALSE(h.accession.has_value());
	CHECK_FALSE(h.description.has_value());
	CHECK_FALSE(h.max_instance_length.has_value());
	CheckNoCutoffs(h);
}

TEST_CASE("HMMER profile - XYPPX: accession, description and cutoffs", "[hmmer]") {
	Initialize();
	ProfileHeader h = ParseProfile(XYPPX_HMM);
	CHECK(h.name == "XYPPX");
	CHECK(h.accession == "PF02162.8");
	CHECK(h.description == "XYPPX repeat");
	CHECK(h.length == 5);
	CHECK(h.alphabet == "amino");
	CHECK(h.format == "HMMER3/e");
	CHECK_FALSE(h.max_instance_length.has_value());
	CHECK(h.sequence_count == 7);
	CHECK(h.effective_sequence_count == F(7.000000));
	CHECK(h.checksum == 149732247u);
	CHECK(h.gathering_threshold_sequence == F(14.30));
	CHECK(h.gathering_threshold_domain == F(0.70));
	CHECK(h.trusted_cutoff_sequence == F(14.40));
	CHECK(h.trusted_cutoff_domain == F(0.70));
	CHECK(h.noise_cutoff_sequence == F(14.20));
	CHECK(h.noise_cutoff_domain == F(4.00));
	CHECK(h.consensus == "GYPPq");
	CHECK(h.msv_mu == F(-5.1112));
	CHECK(h.msv_lambda == F(0.78852));
	CHECK(h.viterbi_mu == F(-4.7526));
	CHECK(h.viterbi_lambda == F(0.78852));
	CHECK(h.forward_tau == F(-2.5557));
	CHECK(h.forward_lambda == F(0.78852));
}

TEST_CASE("HMMER profile - 3box: a HMMER3/f DNA model with MAXL", "[hmmer]") {
	Initialize();
	ProfileHeader h = ParseProfile(BOX3_HMM);
	CHECK(h.name == "3box");
	CHECK_FALSE(h.accession.has_value());
	CHECK_FALSE(h.description.has_value());
	CHECK(h.length == 20);
	CHECK(h.alphabet == "DNA");
	CHECK(h.format == "HMMER3/f");
	CHECK(h.max_instance_length == 75);
	CHECK(h.sequence_count == 22);
	CHECK(h.effective_sequence_count == F(22.000000));
	CHECK(h.checksum == 1874497603u);
	CheckNoCutoffs(h);
	CHECK(h.consensus == "AAAtGAaattGTtTGAtTAg");
	CHECK(h.msv_mu == F(-6.0857));
	CHECK(h.msv_lambda == F(0.74408));
	CHECK(h.viterbi_mu == F(-6.5394));
	CHECK(h.viterbi_lambda == F(0.74408));
	CHECK(h.forward_tau == F(-3.2146));
	CHECK(h.forward_lambda == F(0.74408));
}

TEST_CASE("HMMER profile - with CONS no, HMMER computes the consensus", "[hmmer]") {
	Initialize();
	// HMMER writes "-" in the consensus column of a model without one.
	// Computed: A is node 1's most likely residue (p = 0.72), upper case at
	// p >= 0.5.
	std::string record = Replace(Replace(M1_HMM, "CONS  yes\n", "CONS  no\n"), "      1 A - -\n", "      1 - - -\n");
	CHECK(ParseProfile(record).consensus == "A");
}

TEST_CASE("HMMER profile - HMMER2 models: format, header, and calibrated statistics", "[hmmer]") {
	Initialize();

	SECTION("nucleic") {
		ProfileHeader h = ParseProfile(H2_NT_HMM);
		CHECK(h.name == "ecori");
		CHECK(h.length == 6);
		CHECK(h.alphabet == "DNA");
		CHECK(h.format == "HMMER2.0");
		CHECK(h.sequence_count == 2);
		// HMMER2 has no EFFN or MAXL, and HMMER ignores a HMMER2 CKSUM,
		// which is a different checksum.
		CHECK_FALSE(h.effective_sequence_count.has_value());
		CHECK_FALSE(h.max_instance_length.has_value());
		CHECK_FALSE(h.checksum.has_value());
		CHECK_FALSE(h.accession.has_value());
		CHECK_FALSE(h.description.has_value());
		CheckNoCutoffs(h);
		CHECK(h.consensus == "gaattc");
		REQUIRE(h.msv_mu.has_value());
		CHECK_THAT(*h.msv_mu, WithinAbs(-8.4900, MU_PRINTED));
		CHECK_THAT(*h.msv_lambda, WithinAbs(1.84961, LAMBDA_PRINTED));
		CHECK_THAT(*h.viterbi_mu, WithinAbs(-8.7980, MU_PRINTED));
		CHECK_THAT(*h.viterbi_lambda, WithinAbs(1.84961, LAMBDA_PRINTED));
		CHECK_THAT(*h.forward_tau, WithinAbs(-0.8137, MU_PRINTED));
		CHECK_THAT(*h.forward_lambda, WithinAbs(1.84961, LAMBDA_PRINTED));
	}

	SECTION("amino") {
		ProfileHeader h = ParseProfile(H2_AA_HMM);
		CHECK(h.name == "test");
		CHECK(h.length == 20);
		CHECK(h.alphabet == "amino");
		CHECK(h.format == "HMMER2.0");
		CHECK(h.sequence_count == 10);
		CHECK_FALSE(h.checksum.has_value());
		CHECK(h.consensus == "ACDEFGHiKLmnPqRstvWY");
		REQUIRE(h.msv_mu.has_value());
		CHECK_THAT(*h.msv_mu, WithinAbs(-6.8245, MU_PRINTED));
		CHECK_THAT(*h.msv_lambda, WithinAbs(0.73768, LAMBDA_PRINTED));
		CHECK_THAT(*h.viterbi_mu, WithinAbs(-7.0309, MU_PRINTED));
		CHECK_THAT(*h.viterbi_lambda, WithinAbs(0.73768, LAMBDA_PRINTED));
		CHECK_THAT(*h.forward_tau, WithinAbs(-3.5306, MU_PRINTED));
		CHECK_THAT(*h.forward_lambda, WithinAbs(0.73768, LAMBDA_PRINTED));
	}
}

#if !defined(__EMSCRIPTEN__) || defined(__EMSCRIPTEN_PTHREADS__)
// read_hmm reads files on several threads at once. A HMMER2 read calibrates,
// which runs the SIMD filters with their own random generator: each thread
// must get what one thread alone gets.
TEST_CASE("HMMER profile - reads on several threads at once agree", "[hmmer]") {
	Initialize();
	constexpr int THREADS = 4;
	std::vector<ProfileHeader> results(THREADS);
	std::vector<std::string> errors(THREADS);
	std::vector<std::thread> threads;
	for (int i = 0; i < THREADS; i++) {
		threads.emplace_back([&, i] {
			try {
				results[i] = ParseProfile(H2_AA_HMM);
			} catch (const std::exception &e) {
				errors[i] = e.what();
			}
		});
	}
	for (auto &t : threads) {
		t.join();
	}

	ProfileHeader alone = ParseProfile(H2_AA_HMM);
	REQUIRE(alone.msv_mu.has_value());
	for (int i = 0; i < THREADS; i++) {
		INFO("thread " << i);
		CHECK(errors[i].empty());
		CHECK(results[i].msv_mu == alone.msv_mu);
		CHECK(results[i].msv_lambda == alone.msv_lambda);
		CHECK(results[i].viterbi_mu == alone.viterbi_mu);
		CHECK(results[i].viterbi_lambda == alone.viterbi_lambda);
		CHECK(results[i].forward_tau == alone.forward_tau);
		CHECK(results[i].forward_lambda == alone.forward_lambda);
	}
}
#endif

TEST_CASE("HMMER profile - a model HMMER rejects is an error with HMMER's message and line", "[hmmer]") {
	Initialize();
	const std::string m1 = M1_HMM;

	// A cutoff line before ALPH: HMMER needs the alphabet to know how many
	// values the line has.
	SECTION("GA before ALPH") {
		ExpectRecordError(Replace(m1, "ALPH  amino\n", "GA    25.00 25.00\nALPH  amino\n"), "ALPH must precede GA", 4);
	}
	SECTION("TC before ALPH") {
		ExpectRecordError(Replace(m1, "ALPH  amino\n", "TC    25.00 25.00\nALPH  amino\n"), "ALPH must precede TC", 4);
	}
	SECTION("NC before ALPH") {
		ExpectRecordError(Replace(m1, "ALPH  amino\n", "NC    25.00 25.00\nALPH  amino\n"), "ALPH must precede NC", 4);
	}
	SECTION("ALPH custom") {
		ExpectRecordError(Replace(m1, "ALPH  amino\n", "ALPH  custom\n"), "Unrecognized alphabet type custom", 4);
	}
	// Found at the HMM line, where the header ends.
	SECTION("no LENG line") {
		ExpectRecordError(Replace(m1, "LENG  1\n", ""), "No LENG found for HMM (or LENG <= 0)", 15);
	}
	SECTION("HMMER2 without a NULE line") {
		ExpectRecordError(Replace(H2_NT_HMM, "NULE       0      0      0      0\n", ""),
		                  "No NULE line (null model emissions)", 14);
	}
	SECTION("MM line in a file older than HMMER3/f") {
		ExpectRecordError(Replace(m1, "LENG  1\n", "MM    yes\nLENG  1\n"), "MM line in a file older than HMMER3/f", 3);
	}
	SECTION("ends in the middle of the model") {
		ExpectRecordError(FirstLines(m1, 20), "Premature end of data in main model section, at node 1 (expected 1)",
		                  20);
	}
	SECTION("no closing //") {
		ExpectRecordError(FirstLines(m1, 23), "Premature end of data: missing //?", 23);
	}
}

TEST_CASE("HMMER profile - a record that is not a model is an error", "[hmmer]") {
	Initialize();

	SECTION("not HMMER text") {
		ExpectRecordError(">seq1 not a model\nACGT\n",
		                  "Format tag is '>seq1': unrecognized (expected HMMER3/a to HMMER3/f, or HMMER2.0)", 1);
	}
	SECTION("the format tag's line counts blank and comment lines before it") {
		ExpectRecordError("\n# a comment\n>seq1\nACGT\n",
		                  "Format tag is '>seq1': unrecognized (expected HMMER3/a to HMMER3/f, or HMMER2.0)", 3);
	}
	// There is no line to point at; the record's first is reported.
	SECTION("empty") {
		ExpectRecordError("", "Empty record: no HMMER format tag", 1);
	}
	SECTION("only blank and comment lines") {
		ExpectRecordError("\n   \n# just a comment\n", "Empty record: no HMMER format tag", 1);
	}
	// HMMER itself reports no error here: at the end of a file, it takes a
	// header with no model after it as the end of the file.
	SECTION("a header and no model") {
		ExpectRecordError(FirstLines(M1_HMM, 15), "Premature end of data in the header: no HMM line", 15);
	}
}

TEST_CASE("HMMER profile - how a read's outcome becomes an error", "[hmmer]") {
	Initialize();
	TakeExceptionMessage();

	SECTION("a model") {
		CHECK_NOTHROW(ThrowIfReadFailed(eslOK, true, "", 24));
	}
	// HMMER returns this when allocating a HMMER2 model's background fails.
	SECTION("success, but no model") {
		ExpectParseError([] { ThrowIfReadFailed(eslOK, false, "", 24); },
		                 "HMMER reported success but returned no model", 24);
	}
	SECTION("HMMER's error buffer is the message") {
		ExpectParseError([] { ThrowIfReadFailed(eslEFORMAT, false, "No NAME found for HMM", 24); },
		                 "No NAME found for HMM", 24);
	}
	SECTION("with nothing in the error buffer, the exception's message is") {
		REQUIRE(RaiseModelLengthException(0) == eslEINVAL);
		ExpectParseError([] { ThrowIfReadFailed(eslEINVAL, false, "", 16); }, "model length M=0 out of range 1..100000",
		                 16);
	}
	SECTION("the error buffer comes first, and the exception is taken all the same") {
		REQUIRE(RaiseModelLengthException(0) == eslEINVAL);
		ExpectParseError([] { ThrowIfReadFailed(eslEFORMAT, false, "Failed to allocate body of the new HMM", 16); },
		                 "Failed to allocate body of the new HMM", 16);
		CHECK(TakeExceptionMessage().empty());
	}
	SECTION("with neither, the status") {
		ExpectParseError([] { ThrowIfReadFailed(eslECORRUPT, false, "", 5); },
		                 "HMMER failed with status " + std::to_string(eslECORRUPT) + " and no message", 5);
	}
	SECTION("out of memory") {
		CHECK_THROWS_AS(ThrowIfReadFailed(eslEMEM, false, "", 1), std::bad_alloc);
	}
	// HMMER's parser ignores what some calls it makes return (p7_hmm_SetName, SetAccession,
	// SetDescription), so a read can succeed without the field such a call failed to set. An
	// exception raised during the read says so.
	SECTION("success, with an exception raised during the read") {
		REQUIRE(RaiseModelLengthException(0) == eslEINVAL);
		ExpectParseError([] { ThrowIfReadFailed(eslOK, true, "", 24); }, "model length M=0 out of range 1..100000", 24);
		CHECK(TakeExceptionMessage().empty());
	}
	SECTION("success, with an allocation failure raised during the read") {
		esl_exception(eslEMEM, FALSE, const_cast<char *>(__FILE__), __LINE__,
		              const_cast<char *>("malloc of size %d failed"), 8);
		CHECK_THROWS_AS(ThrowIfReadFailed(eslOK, true, "", 1), std::bad_alloc);
		CHECK(TakeExceptionMessage().empty());
	}
}

// Easel's buffer reader sizes its line buffer as twice a line's length, as an int
// (esl_fileparser.c, nextline()), so a line of INT_MAX / 2 bytes or more would overflow it.
// A record no longer than that holds no such line.
TEST_CASE("HMMER profile - the largest record holds no line Easel cannot buffer", "[hmmer]") {
	CHECK(miint::hmmer::MAX_RECORD_BYTES <= static_cast<size_t>(std::numeric_limits<int>::max() / 2));
}

// What a caller splitting a file into records relies on, as HMMER reads lines (p7_hmmfile.c's
// format tags, Easel's fileparser with '#' as its comment character).
TEST_CASE("HMMER profile - lines as HMMER reads them", "[hmmer]") {
	using miint::hmmer::FirstWord;
	using miint::hmmer::IsBlankOrComment;
	using miint::hmmer::IsFormatTag;
	using miint::hmmer::IsModelEnd;

	CHECK(FirstWord("HMMER3/f [3.1b2 | February 2015]") == "HMMER3/f");
	CHECK(FirstWord(" \tNAME  fn3") == "NAME");
	CHECK(FirstWord("//") == "//");
	CHECK(FirstWord("").empty());
	CHECK(FirstWord(" \t ").empty());
	CHECK(FirstWord("NAME\rfn3") == "NAME");
	// Easel splits on space, tab, CR and LF only, so a vertical tab or form feed is part of a
	// word, and a line holding one is not blank.
	CHECK(FirstWord("\vNAME  fn3") == "\vNAME");
	CHECK(FirstWord("NAME\f fn3") == "NAME\f");
	CHECK_FALSE(IsBlankOrComment(FirstWord("\v")));

	CHECK(IsBlankOrComment(""));
	CHECK(IsBlankOrComment("#"));
	CHECK(IsBlankOrComment("#comment"));
	CHECK_FALSE(IsBlankOrComment("NAME"));
	CHECK_FALSE(IsBlankOrComment("//"));

	for (const char *tag : {"HMMER3/a", "HMMER3/b", "HMMER3/c", "HMMER3/d", "HMMER3/e", "HMMER3/f", "HMMER2.0"}) {
		CHECK(IsFormatTag(tag));
	}
	for (const char *not_tag : {"HMMER3/g", "HMMER3", "HMMER2", "HMMER", "hmmer3/f", "HMMER3/F", ">seq1", "//", ""}) {
		CHECK_FALSE(IsFormatTag(not_tag));
	}

	CHECK(IsModelEnd("//"));
	CHECK_FALSE(IsModelEnd("///"));
	CHECK_FALSE(IsModelEnd("#//"));
	CHECK_FALSE(IsModelEnd(""));
}
#endif
