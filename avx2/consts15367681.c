#include <stdint.h>
#include "params.h"
#include "consts1536.h"

#define P 7681
#define PINV -7679 // p^-1 mod 2^16
#define MONT -3593 // 2^16 mod p
#define MONT_PINV -9
#define V 17474 // floor(2^27/p + 0.5)
#define SHIFT 2048
#define F 956 // mont^2/512
#define F_PINV -1092 // pinv*FHI

__attribute__((aligned(32)))
const int16_t pdata7681[] = {
#define _16XP 0
  P, P, P, P, P, P, P, P, P, P, P, P, P, P, P, P,

#define _16XPINV 16
  PINV, PINV, PINV, PINV, PINV, PINV, PINV, PINV,
  PINV, PINV, PINV, PINV, PINV, PINV, PINV, PINV,

#define _16XV 32
  V, V, V, V, V, V, V, V, V, V, V, V, V, V, V, V,

#define _16XSHIFT 48
  SHIFT, SHIFT, SHIFT, SHIFT, SHIFT, SHIFT, SHIFT, SHIFT,
  SHIFT, SHIFT, SHIFT, SHIFT, SHIFT, SHIFT, SHIFT, SHIFT,

#define _16XMONT_PINV 64
  MONT_PINV, MONT_PINV, MONT_PINV, MONT_PINV, MONT_PINV, MONT_PINV, MONT_PINV, MONT_PINV,
  MONT_PINV, MONT_PINV, MONT_PINV, MONT_PINV, MONT_PINV, MONT_PINV, MONT_PINV, MONT_PINV,

#define _16XMONT 80
  MONT, MONT, MONT, MONT, MONT, MONT, MONT, MONT,
  MONT, MONT, MONT, MONT, MONT, MONT, MONT, MONT,

#define _16XF_PINV 96
  F_PINV, F_PINV, F_PINV, F_PINV, F_PINV, F_PINV, F_PINV, F_PINV,
  F_PINV, F_PINV, F_PINV, F_PINV, F_PINV, F_PINV, F_PINV, F_PINV,

#define _16XF 112
  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,

#define _ZETAS_PINV 128
   28865,  28865,  10350,  10350, -16425, -16425,   7244,   7244,
    4496,   4496,  -4974,  -4974, -14744, -14744, -22593, -22593,

#define _ZETAS 144
    3777,   3777,   3182,   3182,  -3625,  -3625,   1100,   1100,
   -3696,  -3696,   2194,   2194,  -2456,  -2456,   2495,   2495,

#define _TWIST96 160
      -9,  28165,  23361,  20469,    102, -10298, -18191,  16518,
   -1229,  -7491,  21680,  -1613,  14744,  24359,   1979,  19351,
   -3593,   1029,  -1727,   1525,  -2970,   3014,  -2319,  -2938,
   -2765,   2237,  -2896,  -3149,   2456,  -3801,  -3653,   -617,
   19684, -30170, -23754,  29931,  25938, -31177,  22900, -31492,
   16425, -19095, -12662, -15307,   -486,  32533,  20870, -12926,
    1252,   -474,  -2250,   -277,    338,  -1993,  -3724,   3324,
    3625,    873,  -1398,  -1483,   2586,  -2795,   1414,   2434,
    5836,   2816,  11706,  24044,  -4496,  31748,  -9403, -26382,
  -11578,  12235, -18242, -11100,   7867, -15750,  22295,   2133,
    -308,   2816,  -1606,   1516,   3696,  -3068,  -3771,  -2830,
    1734,  -1589,   -834,   3236,   2235,   3706,   2327,   -427,
  -28865,  -7602,  -5392, -25597,  18694,  25691,   -828, -20520,
  -27721,  19394,   9932, -15904,   4974,  29419,  11894,  -5759,
   -3777,   1614,   2800,  -2557,   -762,   3675,  -2876,    -40,
    1463,   1986,   3788,    480,  -2194,   -789,    630,   1921,
    5845, -25349, -11655,   3575,  -4599, -23489,   8788,  22636,
  -10350,  19726,  25614,  -9488,  -6877,  25426,  20315, -17218,
    3285,   1787,    121,     -9,  -1015,   1599,  -1452,    108,
   -3182,  -3826,   2062,  -1296,   -221,   -174,  -1701,    190,
   16988,  22568,  18361,  10008,  -7244,  -8669, -23728,  10972,
   21390, -27047,  22593,   -597,   5461,  -3114,  -8976,   7167,
    2652,   2088,  -2631,  -2280,  -1100,  -2013,    848,  -3364,
   -2162,   1113,  -2495,   1963,   2901,   2006,   -784,   -513,
      -9,  -7167,   8976,   3114,  -5461,    597, -22593,  27047,
  -21390, -10972,  23728,   8669,   7244, -10008, -18361, -22568,
   -3593,    513,    784,  -2006,  -2901,  -1963,   2495,  -1113,
    2162,   3364,   -848,   2013,   1100,   2280,   2631,  -2088,
  -16988,  17218, -20315, -25426,   6877,   9488, -25614, -19726,
   10350, -22636,  -8788,  23489,   4599,  -3575,  11655,  25349,
   -2652,   -190,   1701,    174,    221,   1296,  -2062,   3826,
    3182,   -108,   1452,  -1599,   1015,      9,   -121,  -1787,
   -5845,   5759, -11894, -29419,  -4974,  15904,  -9932, -19394,
   27721,  20520,    828, -25691, -18694,  25597,   5392,   7602,
   -3285,  -1921,   -630,    789,   2194,   -480,  -3788,  -1986,
   -1463,     40,   2876,  -3675,    762,   2557,  -2800,  -1614,
   28865,  -2133, -22295,  15750,  -7867,  11100,  18242, -12235,
   11578,  26382,   9403, -31748,   4496, -24044, -11706,  -2816,
    3777,    427,  -2327,  -3706,  -2235,  -3236,    834,   1589,
   -1734,   2830,   3771,   3068,  -3696,  -1516,   1606,  -2816,
   -5836,  12926, -20870, -32533,    486,  15307,  12662,  19095,
  -16425,  31492, -22900,  31177, -25938, -29931,  23754,  30170,
     308,  -2434,  -1414,   2795,  -2586,   1483,   1398,   -873,
   -3625,  -3324,   3724,   1993,   -338,    277,   2250,    474,
  -19684, -19351,  -1979, -24359, -14744,   1613, -21680,   7491,
    1229, -16518,  18191,  10298,   -102, -20469, -23361, -28165,
   -1252,    617,   3653,   3801,  -2456,   3149,   2896,  -2237,
    2765,   2938,   2319,  -3014,   2970,  -1525,   1727,  -1029,
      -9,  18455,  -7167, -29667,   8976, -16211,   3114,  14061,
   -5461,  14846,    597,  18856, -22593, -20494,  27047,  20674,
   -3593,  -1513,    513,  -3555,    784,   1709,  -2006,   -787,
   -2901,   -514,  -1963,  -1624,   2495,   3058,  -1113,   3266,
  -21390,  -1237, -10972,   9351,  23728,  -9215,   8669,  -7184,
    7244,  -5358, -10008,  21066, -18361,    768, -22568, -10324,
    2162,   1323,   3364,  -2425,   -848,  -1535,   2013,   1008,
    1100,   1810,   2280,   -438,   2631,    768,  -2088,    -84,
  -16988, -21399,  17218, -18140, -20315,  32704, -25426,  11783,
    6877,   1783,   9488,  -9411, -25614,  24581, -19726,   4479,
   -2652,  -1431,   -190,  -3804,   1701,    -64,    174,      7,
     221,  -1801,   1296,    317,  -2062,  -2555,   3826,  -3201,
   10350,  27158, -22636,   6246,  -8788,   3413,  23489, -16757,
    4599,  14121,  -3575,   -520,  11655,  21561,  25349,  23242,
    3182,   -490,   -108,   3174,   1452,    853,  -1599,   2187,
    1015,   1321,      9,   3576,   -121,    569,  -1787,   1738,
   -5845,  -6638,   5759,  -5418, -11894, -29103, -29419,  -1937,
   -4974,  11476,  15904,   5913,  -9932, -13959, -19394,   5623,
   -3285,    530,  -1921,   -298,   -630,   3153,    789,  -1425,
    2194,   1236,   -480,   1305,  -3788,  -2183,  -1986,   2039,
   27721,   4505,  20520, -16877,    828, -20682, -25691,  -5930,
  -18694,   -375,  25597,  12329,   5392,   1724,   7602,   5955,
   -1463,   -103,     40,  -2029,   2876,    822,  -3675,   -810,
     762,   3209,   2557,   -471,  -2800,   3772,  -1614,  -3773,
      -9,  27875,  28165,  -2073,  23361,  28318,  20469, -24854,
     102,  -6817, -10298,  24880, -18191, -12141,  16518, -29428,
   -3593,   1763,   1029,   2535,  -1727,  -3426,   1525,   2794,
   -2970,   1887,   3014,    304,  -2319,   2707,  -2938,  -2804,
   -1229,  16271,  -7491,  29120,  21680,  14624,  -1613,  25452,
   14744,   1357,  24359, -21766,   1979,  21117,  19351,  22261,
   -2765,    399,   2237,  -3648,  -2896,  -1760,  -3149,   2924,
    2456,   2893,  -3801,  -2310,  -3653,  -1923,   -617,   3317,
   19684, -16279, -30170,   -956, -23754,   8737,  29931,  -4983,
   25938,  -1254, -31177,  11467,  22900,  26228, -31492,  -5742,
    1252,   3689,   -474,  -3004,  -2250,     33,   -277,  -1399,
     338,   1818,  -1993,  -2357,  -3724,   -396,   3324,   1426,
   16425,  15051, -19095,  -6536, -12662,  12943, -15307,   3370,
    -486,  15998,  32533,  12892,  20870, -24249, -12926,  25093,
    3625,   1227,    873,  -2440,  -1398,  -2929,  -1483,  -1750,
    2586,    638,  -2795,  -1444,   1414,  -3257,   2434,  -2043,
    5836,   4633,   2816, -23634,  11706,  28839,  24044,  26561,
   -4496,   9940,  31748,  21467,  -9403, -18387, -26382,   8950,
    -308,     25,   2816,   1966,  -1606,    679,   1516,   1473,
    3696,   -300,  -3068,   -549,  -3771,   -467,  -2830,  -2314,
  -11578,  11792,  12235,   4539, -18242,  24035, -11100,  23668,
    7867, -10426, -15750,  11066,  22295, -26279,   2133, -21877,
    1734,   3600,  -1589,  -1093,   -834,  -2077,   3236,  -2956,
    2235,   2886,   3706,  -2246,   2327,   1881,   -427,  -2933,
      -9,  -1442,  18455, -26800,  -7167, -31586, -29667,  32627,
    8976,   9599, -16211,  12815,   3114,   2022,  14061,  17030,
   -3593,   -418,  -1513,  -2224,    513,   2206,  -3555,  -1677,
     784,   1919,   1709,  -3057,  -2006,  -1050,   -787,  -2426,
   -5461,  -5341,  14846,  18617,    597, -30136,  18856, -24564,
  -22593, -17184, -20494,   9855,  27047, -16553,  20674,  20426,
   -2901,   1315,   -514,  -2375,  -1963,  -1464,  -1624,   2060,
    2495,   -800,   3058,   2175,  -1113,  -3753,   3266,  -1078,
  -21390, -10478,  -1237, -12474, -10972, -19334,   9351,   7508,
   23728,  12355,  -9215,  15563,   8669,  23225,  -7184,  31066,
    2162,  -3310,   1323,    838,   3364,    122,  -2425,  -2732,
    -848,   2627,  -1535,   1739,   2013,   2233,   1008,   1370,
    7244, -20972,  -5358,  11962, -10008,  12534,  21066,  21220,
  -18361, -22875,    768,  -1297, -22568, -12858, -10324,  24718,
    1100,   1556,   1810,  -1350,   2280,   1270,   -438,   2788,
    2631,   -859,    768,   -785,  -2088,    454,    -84,   1166,
  -16988,  12670, -21399, -11920,  17218,  26262, -18140,  14616,
  -20315, -25400,  32704,  21953, -25426,  28378,  11783,  25247,
   -2652,  -2690,  -1431,  -3728,   -190,  -1386,  -3804,   2328,
    1701,   3272,    -64,  -3135,    174,  -1318,      7,   1183,
    6877, -17440,   1783, -26313,   9488,  30579,  -9411, -17602,
  -25614,  -3345,  24581,  25477, -19726,   8558,   4479, -29411,
     221,  -1056,  -1801,   2871,   1296,  -3725,    317,   -194,
   -2062,  -2833,  -2555,  -1659,   3826,   1390,  -3201,  -3299,
      -9, -24266,  -2073,   1689,  20469,  17303,  -6817,  10409,
  -18191, -26441, -29428,  -2722,  -7491,  -6032,  14624, -24052,
   -3593,  -2762,   2535,  -2919,   1525,  -2665,   1887,  -2391,
   -2319,   2743,  -2804,  -1698,   2237,   2160,  -1760,   2572,
   14744, -12013, -21766,  29871,  19351, -15768, -16279, -30426,
  -23754,  12073,  -4983, -15358, -31177,   3558,  26228,  12611,
    2456,   2835,  -2310,  -2385,   -617,  -3480,   3689,   -730,
   -2250,   -727,  -1399,      2,  -1993,    486,   -396,   2883,
   16425, -15810,  -6536,  24743, -15307, -16655,  15998,  16092,
   20870, -21860,  25093,  -3456,   2816,  12269,  28839,  32329,
    3625,   1598,  -2440,  -3417,  -1483,   -783,    638,   1756,
    1414,  -3428,  -2043,  -3456,   2816,  -2579,    679,   3145,
   -4496,  -8472,  21467, -27201, -26382,   9360,  11792, -19317,
  -18242,  24573,  23668,   7415, -15750,  32252, -26279, -27132,
    3696,   3816,   -549,  -2113,  -2830,   1168,   3600,   -373,
    -834,   1533,  -2956,   3831,   3706,   1532,   1881,   3588,
  -28865,  25946,  -1724,  13541, -25597,  13575,   5930,  21868,
    -828,   5529,  -4505, -32738,  19394, -25511,  -5913,  26663,
   -3777,  -3750,  -3772,   2789,  -2557,   1799,    810,   -660,
   -2876,    921,    103,   1054,   1986,   2649,  -1305,  -1497,
    4974,  -8908,  29103,  -1869,  -5759,   4693, -23242,  26220,
  -11655,  14351, -14121,  13933, -23489, -22124,  -6246,  -2193,
   -2194,  -2764,  -3153,  -3405,   1921,   2133,  -1738,   3692,
     121,  -1521,  -1321,   -915,   1599,    404,  -3174,  -1681,
      -9, -26800, -29667,   9599,   3114,  17030,  14846, -30136,
  -22593,   9855,  20674, -10478, -10972,   7508,  -9215,  23225,
   -3593,  -2224,  -3555,   1919,  -2006,  -2426,   -514,  -1464,
    2495,   2175,   3266,  -3310,   3364,  -2732,  -1535,   2233,
    7244,  11962,  21066, -22875, -22568,  24718, -21399,  26262,
  -20315,  21953,  11783, -17440,   9488, -17602,  24581,   8558,
    1100,  -1350,   -438,   -859,  -2088,   1166,  -1431,  -1386,
    1701,  -3135,      7,  -1056,   1296,   -194,  -2555,   1390,
   10350,   2193,   6246,  22124,  23489, -13933,  14121, -14351,
   11655, -26220,  23242,  -4693,   5759,   1869, -29103,   8908,
    3182,   1681,   3174,   -404,  -1599,    915,   1321,   1521,
    -121,  -3692,   1738,  -2133,  -1921,   3405,   3153,   2764,
   -4974, -26663,   5913,  25511, -19394,  32738,   4505,  -5529,
     828, -21868,  -5930, -13575,  25597, -13541,   1724, -25946,
    2194,   1497,   1305,  -2649,  -1986,  -1054,   -103,   -921,
    2876,    660,   -810,  -1799,   2557,  -2789,   3772,   3750,
   28865,  27132,  26279, -32252,  15750,  -7415, -23668, -24573,
   18242,  19317, -11792,  -9360,  26382,  27201, -21467,   8472,
    3777,  -3588,  -1881,  -1532,  -3706,  -3831,   2956,  -1533,
     834,    373,  -3600,  -1168,   2830,   2113,    549,  -3816,
    4496, -32329, -28839, -12269,  -2816,   3456, -25093,  21860,
  -20870, -16092, -15998,  16655,  15307, -24743,   6536,  15810,
   -3696,  -3145,   -679,   2579,  -2816,   3456,   2043,   3428,
   -1414,  -1756,   -638,    783,   1483,   3417,   2440,  -1598,
      -9,  -7756,  27875, -24266,  28165, -22713,  -2073,  15887,
   23361,   1689,  28318, -14181,  20469,  -6083, -24854,  17303,
   -3593,  -1612,   1763,  -2762,   1029,  -1721,   2535,     15,
   -1727,  -2919,  -3426,  -3429,   1525,   3645,   2794,  -2665,
     102,  27533,  -6817,  29044, -10298,  10409,  24880,   5964,
  -18191, -20273, -12141, -26441,  16518,   7466, -29428, -11032,
   -2970,  -3699,   1887,   2420,   3014,  -2391,    304,   -180,
   -2319,  -3377,   2707,   2743,  -2938,   2346,  -2804,   1256,
   -1229,  -2722,  16271, -20844,  -7491,   6160,  29120,  -6032,
   21680, -18873,  14624, -10384,  -1613, -24052,  25452,   1314,
   -2765,  -1698,    399,   1684,   2237,  -2032,  -3648,   2160,
   -2896,   2119,  -1760,  -2192,  -3149,   2572,   2924,    290,
   14744,  32661,   1357, -12013,  24359,  -8387, -21766,   6851,
    1979,  29871,  21117,  -6467,  19351,  26484,  22261, -15768,
    2456,  -2667,   2893,   2835,  -3801,   1341,  -2310,  -2877,
   -3653,  -2385,  -1923,   3261,   -617,   -140,   3317,  -3480,
   19684,   1280, -16279,  13088, -30170, -30426,   -956, -16680,
  -23754, -30776,   8737,  12073,  29931,   9872,  -4983,  -7397,
    1252,   1280,   3689,  -3296,   -474,   -730,  -3004,   3800,
   -2250,  -2104,     33,   -727,   -277,   1680,  -1399,   3355,
   25938, -15358,  -1254, -25989, -31177, -28105,  11467,   3558,
   22900, -23907,  26228, -13805, -31492,  12611,  -5742,  23233,
     338,      2,   1818,   1147,  -1993,   1079,  -2357,    486,
   -3724,   2205,   -396,   1043,   3324,   2883,   1426,  -1855,
      -9,   -111,  -1442, -18745,  18455, -22226, -26800, -20716,
   -7167, -27636, -31586, -17406, -29667,   7551,  32627,  30938,
   -3593,   -623,   -418,   2247,  -1513,   3374,  -2224,   1812,
     513,  -1012,   2206,  -2046,  -3555,   -129,  -1677,   1242,
    8976, -14385,   9599,  -6288, -16211, -14138,  12815, -30008,
    3114, -25051,   2022,  26288,  14061, -13814,  17030,  24786,
     784,   2511,   1919,   1904,   1709,   -826,  -3057,  -1336,
   -2006,  -3035,  -1050,   1712,   -787,  -2550,  -2426,   -814,
   -5461,  -5452,  -5341,  -3899,  14846,  -3609,  18617, -20119,
     597,   7764, -30136,   1450,  18856, -17013, -24564,   8345,
   -2901,    692,   1315,   1733,   -514,    999,  -2375,   -151,
   -1963,  -2476,  -1464,  -3670,  -1624,   1931,   2060,   3737,
  -22593, -31569, -17184, -26783, -20494,  -4283,   9855,  -2961,
   27047,  23933, -16553, -18575,  20674,   6612,  20426,   3396,
    2495,   1711,   -800,  -2719,   3058,   1349,   2175,  -2449,
   -1113,    893,  -3753,  -2703,   3266,  -3628,  -1078,   1348,
  -21390, -15930, -10478,  -5136,  -1237, -16083, -12474, -31091,
  -10972, -11570, -19334,  10802,   9351,  -9505,   7508,  32073,
    2162,  -2618,  -3310,   3056,   1323,   1837,    838,   3213,
    3364,  -2354,    122,   1586,  -2425,   -801,  -2732,   2889,
   23728, -19215,  12355,  29539,  -9215,  11280,  15563,   5708,
    8669, -18378,  23225, -25759,  -7184, -27858,  31066,  10640,
    -848,  -3343,   2627,   3427,  -1535,   3088,   1739,   -436,
    2013,   3126,   2233,  -1695,   1008,  -2258,   1370,   2448,
      -9, -27073,  15887,  12261, -24854, -23088,  10409,  -1399,
   16518, -15972, -20844, -12884,  14624,   3191,  32661,  22039,
   -3593,  -1985,     15,   1509,   2794,   1488,  -2391,   2185,
   -2938,   2460,   1684,  -2644,  -1760,   -393,  -2667,   2071,
    1979, -10580, -15768, -26458,   -956, -17491,   9872,  -3046,
  -31177, -29974, -13805, -25648,  15051, -19428,  24743,  -1177,
   -3653,   -340,  -3480,   3238,  -3004,    429,   1680,     26,
   -1993,  -2326,   1043,  -1072,   1227,   -996,  -3417,   3431,
    -486,  29718, -11945, -22030,  25093,  -5307,   3507, -14223,
   24044,   7082,  -8472, -13472, -18387, -14718,  26689,  11024,
    2586,   2070,    855,   1522,  -2043,    325,   1971,   1649,
    1516,   1962,   3816,   2912,   -467,    642,   1601,   2832,
  -18242, -13233,  18754,   -802,  11066, -13950, -27132,  22986,
   -7602,  -4565,   -461, -20136,   5930,   6723, -31416,  -1578,
    -834,   3663,   1346,    222,  -2246,   1410,   3588,   1482,
    1614,  -2005,  -1997,    344,    810,  -3005,  -2744,   3542,
  -27721, -10025, -25511, -10537, -11476,  25179,   3276, -24291,
   -5759,  10443, -24180,  18506,    520,  13037,  13933, -27013,
    1463,   2775,   2649,   2263,  -1236,   3163,  -2868,   1821,
    1921,  -3381,   2444,  -2998,  -3576,  -1811,   -915,    123,
    8788,  32124,  20392,  19820, -24581,  -8728,  26313,   -529,
   25426,   1954, -26262,  31893,  10324,  -9983, -21220, -24436,
   -1452,   1404,    -88,  -2708,   2555,   3560,  -2871,    -17,
    -174,    930,   1386,  -3435,     84,  -2303,  -2788,   2188,
      -9, -22226, -31586,  30938, -16211, -25051,  17030,  -3899,
     597, -17013, -17184,  -2961,  20674, -15930, -12474,  10802,
   -3593,   3374,   2206,   1242,   1709,  -3035,  -2426,   1733,
   -1963,   1931,   -800,  -2449,   3266,  -2618,    838,   1586,
   23728,  11280,  23225,  10640,  -5358,    964,  21220,  30306,
  -22568,  -3140,  12670, -23882, -18140,  -1954,  21953, -24300,
    -848,   3088,   2233,   2448,   1810,  -1084,   2788,  -3486,
   -2088,  -1092,  -2690,  -2378,  -3804,   -930,  -3135,  -1772,
    6877,  23182,  30579, -32218,  24581,   5700, -29411,  -2952,
  -22636,  15657,  22124,  27013, -16757,  -5751,  27124,  10042,
     221,   -370,  -3725,  -2522,  -2555,   3652,  -3299,   1144,
    -108,   2857,   -404,   -123,   2187,  -2167,    500,  -3270,
   11655,  18148,  24180,   9633,  -6638,   9334,   1869,  17926,
  -29419, -25179,  11314, -18942,   5913,   1962,    273,  -9650,
    -121,   -284,  -2444,    929,    530,  -1930,   3405,  -1530,
     789,  -3163,   2098,  -3582,   1305,  -3158,   -239,   -434,
   27721,  -6971,  -5529,  15187, -20682,  -6297, -19112,  20136,
   25597,  29206,  -6194, -14539,   5955, -17977,  27132,  32269,
   -1463,  -1339,   -921,  -2733,    822,  -1689,   1368,   -344,
    2557,   1558,   3022,   -715,  -3773,   3015,  -3588,   1037,
  -22295, -12790, -25144, -30853, -23668,  13233,   1271, -30597,
  -12235, -22218,  -9360,  -3012,  18387,  -8839, -23455, -20930,
   -2327,  -1526,   3528,  -3717,   2956,  -3663,  -2313,  -3461,
    1589,   -714,  -1168,   -964,    467,   2937,    609,  -3522,
      -9, -30844, -24266,  32414,  -2073, -23984,   1689,  12261,
   20469,   4573,  17303,  30255,  -6817,  -2807,  10409,  11996,
   -3593,   -124,  -2762,    670,   2535,    592,  -2919,   1509,
    1525,  -2083,  -2665,  -2001,   1887,    777,  -2391,  -2340,
  -18191, -26766, -26441,  31518, -29428, -15972,  -2722,  -8831,
   -7491, -14650,  -6032,  16791,  14624, -20964, -24052,  17082,
   -2319,  -3214,   2743,   -226,  -2804,   2460,  -1698,  -1151,
    2237,  -1338,   2160,  -3177,  -1760,  -2532,   2572,   3770,
   14744,  17636, -12013,  22039, -21766,  25733,  29871, -18541,
   19351,  27243, -15768,  16638, -16279,   1007, -30426, -20238,
    2456,   -796,   2835,   2071,  -2310,  -1403,  -2385,  -3693,
    -617,  -2965,  -3480,   1278,   3689,   1519,   -730,   3314,
  -23754, -17491,  12073,  -2739,  -4983,   9513, -15358, -10179,
  -31177,  18003,   3558,  16885,  26228, -16186,  12611, -25648,
   -2250,    429,   -727,  -1203,  -1399,  -3287,      2,   -451,
   -1993,     83,    486,  -2059,   -396,  -2874,   2883,  -1072,
   16425,   -947, -15810,  -6510,  -6536,  32004,  24743,  -9087,
  -15307, -21723, -16655,  20128,  15998,  29718,  16092, -24214,
    3625,    589,   1598,    658,  -2440,  -2812,  -3417,  -1407,
   -1483,    293,   -783,   3744,    638,   2070,   1756,   3434,
   20870,  12449, -21860,  14129,  25093,  10333,  -3456,  25563,
    2816,  20435,  12269, -14223,  28839, -15111,  32329,  17167,
    1414,   3745,  -3428,  -2767,  -2043,   3677,  -3456,   3547,
    2816,   2515,  -2579,   1649,    679,  -3335,   3145,   1295,
      -9, -18745, -26800, -27636, -29667,  30938,   9599, -14138,
    3114,  26288,  17030,  -5452,  14846, -20119, -30136, -17013,
   -3593,   2247,  -2224,  -1012,  -3555,   1242,   1919,   -826,
   -2006,   1712,  -2426,    692,   -514,   -151,  -1464,   1931,
  -22593, -26783,   9855,  23933,  20674,   3396, -10478, -16083,
  -10972,  10802,   7508, -19215,  -9215,   5708,  23225, -27858,
    2495,  -2719,   2175,    893,   3266,   1348,  -3310,   1837,
    3364,   1586,  -2732,  -3343,  -1535,   -436,   2233,  -2258,
    7244, -10495,  11962,    964,  21066,  13711, -22875,   9983,
  -22568,  29453,  24718, -24232, -21399, -23882,  26262,  26330,
    1100,  -2815,  -1350,  -1084,   -438,  -2161,   -859,   2303,
   -2088,  -1779,   1166,  -3752,  -1431,  -2378,  -1386,  -3366,
  -20315,  -2526,  21953,  -2858,  11783,    529, -17440,  23182,
    9488,   4317, -17602,  -5299,  24581,   3524,   8558,  -7304,
    1701,  -3550,  -3135,   2262,      7,     17,  -1056,   -370,
    1296,  -2339,   -194,  -3763,  -2555,   1476,   1390,  -3208,
   10350,  -2952,   2193, -32124,   6246,  24530,  22124, -21168,
   23489,  28958, -13933,  -5751,  14121,  24931, -14351,  -6766,
    3182,   1144,   1681,  -1404,   3174,  -1070,   -404,   3408,
   -1599,  -2786,    915,  -2167,   1321,  -1181,   1521,    402,
   11655, -18506, -26220,   1860,  23242,   9633,  -4693, -20759,
    5759,   4624,   1869, -23549, -29103,  22943,   8908, -25179,
    -121,   2998,  -3692,   -188,   1738,    929,  -2133,   -791,
   -1921,  -3568,   3405,   -509,   3153,  -1121,   2764,  -3163,
      -9,  12227, -22713, -23984,  28318, -13549,  17303, -23088,
  -10298,  11996, -20273,  18131, -29428,  30153, -20844, -14650,
   -3593,   2499,  -1721,    592,  -3426,   1299,  -2665,   1488,
    3014,  -2340,  -3377,    211,  -2804,    969,   1684,  -1338,
   21680,  -3916, -24052,  -4940,   1357,  22039,   6851,  12380,
   19351,  -8370,   1280,   1007,   -956,  -6868,  12073,    759,
   -2896,   2228,   2572,   1204,   2893,   2071,  -2877,  -1956,
    -617,    846,   1280,   1519,  -3004,   3372,   -727,  -2825,
   25938, -10179, -28105, -29974,  26228, -32670,  23233,   -947,
  -19095, -30750,  24743,  -2380,   3370,  20128,  -6886, -29010,
     338,   -451,   1079,  -2326,   -396,   -926,  -1855,    589,
     873,    994,  -3417,   3764,  -1750,   3744,  -3814,  -3410,
   20870, -22030, -19078,  10333,   4633, -19948,  12269, -18310,
   24044,  17167, -24069,  19513,  21467, -25929,  -2048, -14718,
    1414,   1522,    378,   3677,     25,   2580,  -2579,   1146,
    1516,   1295,   3067,  -1479,   -549,   3255,  -2048,    642,
  -11578,  10990, -19317,  24999,  24035,  29394,  18754,  -4471,
  -15750,  26057,  15256,  27730, -21877,  22986,  25946, -22406,
    1734,   3822,   -373,  -3161,  -2077,   3794,   1346,   -887,
    3706,  -3127,   2968,   2130,  -2933,   1482,  -3750,  -2950,
   -5392,  -5990,  -2090, -19121,   5930,  22141,  -8797,  22790,
  -20520, -20511, -32738, -10025,  13959, -14360,   -811, -18114,
    2800,  -2918,   3030,  -2225,    810,   -899,  -2141,   3334,
     -40,   3553,   1054,   2775,   2183,  -2072,  -3371,   -706,
      -9, -20716,  32627, -14138,  14061,  -3899, -30136, -31569,
   27047,   3396, -12474,  -9505,  -9215, -25759, -20972,    964,
   -3593,   1812,  -1677,   -826,   -787,   1733,  -1464,   1711,
   -1113,   1348,    838,   -801,  -1535,  -1695,   1556,  -1084,
  -18361, -16860,  24718, -16041, -18140,  -2526,  28378,  23865,
    9488, -32218,  25477,  -7304,  27158,  10563,  22124, -22320,
    2631,  -2524,   1166,  -3241,  -3804,  -3550,  -1318,   2873,
    1296,  -2522,  -1659,  -3208,   -490,    835,   -404,   2256,
    4599,  24931, -22423,  18148,  23242,   4531,  -9727, -23549,
  -29419,   4659, -26663,  11331, -13959,  -9650,  31825,   4616,
    1015,  -1181,  -2455,   -284,   1738,   2995,  -2047,   -509,
     789,   3123,   1497,   1603,  -2183,   -434,  -1455,    520,
     828, -22141, -19112,  -4880,  12329, -14991, -25946, -17977,
   -2133,  -1715,  30341,   4471, -23668,   8234,   2688, -27986,
    2876,    899,   1368,   3312,   -471,    881,   3750,   3015,
     427,   -179,   3205,    887,   2956,   3114,   2688,  -2386,
   11578,  -3012,  27201, -16928,  -9940, -17167,   -171, -21109,
   -2816, -20614,  19078, -12449, -12892,  23976,  16655,   2380,
   -1734,   -964,   2113,   -544,    300,  -1295,  -2731,  -2165,
   -2816,  -1158,   -378,  -3745,   1444,   3496,    783,  -3764,
   12662,   6510,  12312,   9112, -26228,  19283,  25989,  -9513,
  -29931,   6868,  30426,  15025, -22261,  18541,  -6851,  11007,
    1398,   -658,     24,  -3176,    396,   1363,  -1147,   3287,
     277,  -3372,    730,  -1871,  -3317,   3693,   2877,   3327,
      -9,  30247,  -7756, -30844,  27875,  12227, -24266, -27073,
   28165,  32414, -22713, -26953,  -2073,   9923,  15887, -23984,
   -3593,   2087,  -1612,   -124,   1763,   2499,  -2762,  -1985,
    1029,    670,  -1721,   2231,   2535,    195,     15,    592,
   23361,  21962,   1689, -25076,  28318,  12261, -14181,   3950,
   20469, -13549,  -6083,   4573, -24854,  28336,  17303,   1331,
   -1727,    458,  -2919,   1548,  -3426,   1509,  -3429,  -3218,
    1525,   1299,   3645,  -2083,   2794,   3760,  -2665,   -205,
     102,  30255,  27533, -23088,  -6817, -15648,  29044,  -2807,
  -10298,   4249,  10409,  -4241,  24880,  11996,   5964,  25665,
   -2970,  -2001,  -3699,   1488,   1887,    736,   2420,    777,
    3014,   -359,  -2391,  -3729,    304,  -2340,   -180,    577,
  -18191,  -1399, -20273, -26766, -12141, -16058, -26441,  18131,
   16518,  31518,   7466,  10657, -29428, -12346, -11032, -15972,
   -2319,   2185,  -3377,  -3214,   2707,  -2746,   2743,    211,
   -2938,   -226,   2346,   1953,  -2804,    966,   1256,   2460,
   -1229,  30153,  -2722,  14914,  16271,  -8831, -20844, -31851,
   -7491,  14547,   6160, -14650,  29120, -12884,  -6032,  19701,
   -2765,    969,  -1698,  -2494,    399,  -1151,   1684,  -1643,
    2237,  -3373,  -2032,  -1338,  -3648,  -2644,   2160,    757,
   21680,  16791, -18873,  -6493,  14624,  -3916, -10384, -20964,
   -1613,  15000, -24052,   3191,  25452,  17082,   1314,  -4940,
   -2896,  -3177,   2119,    163,  -1760,   2228,  -2192,  -2532,
   -3149,   2712,   2572,   -393,   2924,   3770,    290,   1204,

#define _TWIST12 3232
      -9,  -1229,  19684,  16425,     -9,  -5461, -21390,   7244,
   -3593,  -2765,   1252,   3625,  -3593,  -2901,   2162,   1100,
    5836, -11578, -28865, -27721, -16988,   6877,  10350,   4599,
    -308,   1734,  -3777,   1463,  -2652,    221,   3182,   1015,
    5845, -10350,  16988,  21390,  -5845,  -4974,  27721, -18694,
    3285,  -3182,   2652,  -2162,  -3285,   2194,  -1463,    762,
      -9, -21390, -16988,  10350,     -9,    102,  -1229,  14744,
   -3593,   2162,  -2652,   3182,  -3593,  -2970,  -2765,   2456,
   -5845,  27721,  28865,  11578,  19684,  25938,  16425,   -486,
   -3285,  -1463,   3777,  -1734,   1252,    338,   3625,   2586,
   -5836, -16425, -19684,   1229,   5836,  -4496, -11578,   7867,
     308,  -3625,  -1252,   2765,   -308,   3696,   1734,   2235,
      -9,   8976,  -5461, -22593,     -9, -22593,   7244, -20315,
   -3593,    784,  -2901,   2495,  -3593,   2495,   1100,   1701,
  -21390,  23728,   7244, -18361,  10350,  11655,  -4974,    828,
    2162,   -848,   1100,   2631,   3182,   -121,   2194,   2876,
  -16988, -20315,   6877, -25614,  28865,  18242,   4496, -20870,
   -2652,   1701,    221,  -2062,   3777,    834,  -3696,  -1414,
      -9, -18191,  14744, -23754,     -9,  23361,    102, -18191,
   -3593,  -2319,   2456,  -2250,  -3593,  -1727,  -2970,  -2319,
   16425,  20870,  -4496, -18242,  -1229,  21680,  14744,   1979,
    3625,   1414,   3696,   -834,  -2765,  -2896,   2456,  -3653,
  -28865,   -828,   4974, -11655,  19684, -23754,  25938,  22900,
   -3777,  -2876,  -2194,    121,   1252,  -2250,    338,  -3724,
};
