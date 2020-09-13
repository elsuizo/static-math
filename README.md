![MIT](https://img.shields.io/badge/license-MIT-blue.svg)
[![Documentation](https://docs.rs/static-math/badge.svg)](https://docs.rs/static-math)
[![crates.io](https://img.shields.io/crates/v/static-math.svg)](https://crates.io/crates/static-math)

# Static Math in Rust programming language

"*Simple things should be simple, complex things should be possible*" Alan Kay.

- This crate take advantage of the static arrays in Rust for fast operations in
stack memory.

- We use a tuple to indexing elements: `m[(i, j)]` allowing nice interface with the `match` feature of Rust

- No `unsafe` code :ballot_box_with_check:

- Could be optimize more with the use of SIMD

- This crate could be used in an `no-std` environment.

- You can visualize the matrices

```text
inverse:
|-0.54    0.58    0.67    -0.08   -0.17    -1.18|
|2.16     -1.53   -2.44   0.44    0.32      3.77|
|0.21     -0.42   -0.39   0.15    0.20      0.62|
|0.70     -0.24   -0.53   0.20    -0.21     0.73|
|0.85     -0.47   -0.83   0.11    0.11      1.20|
|-3.91    2.47    4.17    -0.87   -0.31    -6.08|
```

- The determinant of the matrices are evaluated "in-place" without loops and code
bifurcations

- The use cases can be: Robotics, Game programming, Simulations ...etc.

The matrix types `Mnn` (where `n=2..6`) implements the Methods from the
`LinearAlgebra` trait:

 - `det()`: Determinant of the matrix
 - `inverse()`: Inverse of the matrix
 - `qr()`: QR decomposition of the matrix
 - `norm2()`: norm of the matrix
 - `transpose()`: transpose of the matrix
 - `trace()`: trace of the matrix
 - `shape()`: shape of the matrix

## Benchmarks

Using the criterion crate:

https://github.com/bheisler/criterion.rs

this are the results for matrixs inverse operations(in a very old machine)


```text
inverse 6x6             time:   [9.6090 us 9.6128 us 9.6172 us]
                        change: [-3.2723% -3.0278% -2.8038%] (p = 0.00 < 0.05)
                        Performance has improved.
Found 5 outliers among 100 measurements (5.00%)
  1 (1.00%) low mild
  1 (1.00%) high mild
  3 (3.00%) high severe

inverse 4x4             time:   [98.560 ns 98.605 ns 98.677 ns]
                        change: [-5.4359% -3.2101% -1.4680%] (p = 0.00 < 0.05)
                        Performance has improved.
Found 15 outliers among 100 measurements (15.00%)
  5 (5.00%) high mild
  10 (10.00%) high severe
```

you can look the bench here: [bench](benches/bench_inverse.rs)

Others benches comparing the performance with others crates are in this repo: https://github.com/bitshifter/mathbench-rs

with the following results:


| benchmark                  |          glam   |        cgmath   |      nalgebra   |       euclid   |           vek   |    pathfinder   |   static-math   |   ultraviolet   |
|:---------------------------|----------------:|----------------:|----------------:|---------------:|----------------:|----------------:|----------------:|----------------:|
| euler 2d x10000            |    __7.555 us__ |    __7.521 us__ |      16.38 us   |     11.86 us   |    __7.513 us__ |      9.806 us   |      11.83 us   |    __7.499 us__ |
| euler 3d x10000            |     __16.2 us__ |      25.04 us   |      106.3 us   |     25.05 us   |      25.16 us   |      16.76 us   |      25.03 us   |      25.05 us   |
| matrix2 determinant        |   __2.0332 ns__ |   __2.0311 ns__ |   __2.0250 ns__ |      N/A       |   __2.0209 ns__ |   __2.0323 ns__ |   __2.0254 ns__ |       N/A       |
| matrix2 inverse            |   __2.6114 ns__ |     3.0331 ns   |     2.9792 ns   |      N/A       |       N/A       |     2.7550 ns   |     3.0132 ns   |       N/A       |
| matrix2 mul matrix2        |     2.6047 ns   |   __2.5346 ns__ |   __2.5426 ns__ |      N/A       |     8.7573 ns   |   __2.5381 ns__ |     2.6028 ns   |     2.9668 ns   |
| matrix2 mul vector2 x1     |     2.6592 ns   |     2.6104 ns   |     2.6214 ns   |      N/A       |     4.2512 ns   |   __2.0663 ns__ |     2.8674 ns   |     2.6172 ns   |
| matrix2 mul vector2 x100   |   245.2897 ns   |   233.7149 ns   |   238.7395 ns   |      N/A       |   399.3148 ns   | __218.4107 ns__ |   260.6645 ns   |   234.7099 ns   |
| matrix2 return self        |   __2.4740 ns__ |     2.5994 ns   |     2.5968 ns   |      N/A       |     2.5969 ns   |   __2.4607 ns__ |     2.5928 ns   |     2.5974 ns   |
| matrix2 transpose          |   __2.0852 ns__ |   __2.0814 ns__ |     2.3426 ns   |      N/A       |   __2.1053 ns__ |       N/A       |   __2.0829 ns__ |       N/A       |
| matrix3 determinant        |   __3.3675 ns__ |   __3.4261 ns__ |   __3.3780 ns__ |      N/A       |   __3.4479 ns__ |       N/A       |   __3.4375 ns__ |       N/A       |
| matrix3 inverse            |    11.4209 ns   |   __8.3701 ns__ |     9.4315 ns   |      N/A       |       N/A       |       N/A       |     9.1710 ns   |    20.1731 ns   |
| matrix3 mul matrix3        |   __5.8501 ns__ |     6.5350 ns   |     9.8196 ns   |      N/A       |    47.9203 ns   |       N/A       |     9.5170 ns   |     6.5211 ns   |
| matrix3 mul vector3 x1     |   __3.9266 ns__ |     4.3876 ns   |     4.3333 ns   |      N/A       |    16.0858 ns   |       N/A       |     4.4220 ns   |     4.3304 ns   |
| matrix3 mul vector3 x100   |   __0.4372 us__ |   __0.4416 us__ |     0.4594 us   |      N/A       |       1.59 us   |       N/A       |      0.454 us   |   __0.4425 us__ |
| matrix3 return self        |   __4.8566 ns__ |   __4.8401 ns__ |   __4.8226 ns__ |      N/A       |   __4.8340 ns__ |       N/A       |   __4.8303 ns__ |   __4.8383 ns__ |
| matrix3 transpose          |   __5.7688 ns__ |   __5.6980 ns__ |     8.1508 ns   |      N/A       |   __5.6910 ns__ |       N/A       |   __5.6936 ns__ |   __5.6766 ns__ |
| matrix4 determinant        |   __8.3724 ns__ |    11.1604 ns   |    52.8697 ns   |   16.0723 ns   |    17.5301 ns   |       N/A       |    16.1402 ns   |       N/A       |
| matrix4 inverse            |  __21.3281 ns__ |    38.5833 ns   |    64.5172 ns   |   61.2347 ns   |   275.5253 ns   |       N/A       |    48.0641 ns   |    37.1436 ns   |
| matrix4 mul matrix4        |   __7.5043 ns__ |     8.3723 ns   |     9.4094 ns   |   10.1761 ns   |    90.7185 ns   |       N/A       |    20.6424 ns   |     8.4072 ns   |
| matrix4 mul vector4 x1     |   __3.3645 ns__ |     3.7273 ns   |     3.7251 ns   |      N/A       |    24.2185 ns   |       N/A       |     6.1311 ns   |     3.7524 ns   |
| matrix4 mul vector4 x100   |   __0.6105 us__ |   __0.6237 us__ |   __0.6202 us__ |      N/A       |      2.402 us   |       N/A       |     0.7044 us   |   __0.6202 us__ |
| matrix4 return self        |     6.8863 ns   |     7.1298 ns   |   __6.6961 ns__ |      N/A       |   __6.7079 ns__ |       N/A       |   __6.6772 ns__ |   __6.7079 ns__ |
| matrix4 transpose          |   __5.7312 ns__ |    10.1612 ns   |    14.9424 ns   |      N/A       |    10.2015 ns   |       N/A       |    10.1996 ns   |    10.2391 ns   |
| rotation3 inverse          |   __2.1867 ns__ |     2.9086 ns   |     2.8853 ns   |    2.9092 ns   |     2.8987 ns   |       N/A       |       N/A       |     2.9064 ns   |
| rotation3 mul rotation3    |   __3.3422 ns__ |     4.3602 ns   |     7.0680 ns   |    7.7111 ns   |     8.9616 ns   |       N/A       |       N/A       |    18.4088 ns   |
| rotation3 mul vector3      |   __6.6977 ns__ |   __6.7831 ns__ |     6.9924 ns   |    6.9801 ns   |    32.8778 ns   |       N/A       |       N/A       |    13.5267 ns   |
| rotation3 return self      |   __2.4622 ns__ |     2.5983 ns   |     2.6021 ns   |      N/A       |     2.5989 ns   |       N/A       |       N/A       |     2.5980 ns   |
| transform point2 x1        |     3.8946 ns   |     2.8843 ns   |     4.6543 ns   |    3.2271 ns   |    17.0089 ns   |   __2.3608 ns__ |       N/A       |       N/A       |
| transform point2 x100      |     0.4265 us   |     0.3677 us   |     0.4632 us   |   __0.322 us__ |      1.712 us   |   __0.3206 us__ |       N/A       |       N/A       |
| transform point3 x1        |     4.9958 ns   |     6.3712 ns   |     6.6426 ns   |    6.1114 ns   |    24.8255 ns   |   __3.1011 ns__ |       N/A       |       N/A       |
| transform point3 x100      |   __0.6261 us__ |     0.7418 us   |     0.7447 us   |    0.7296 us   |      2.507 us   |   __0.6295 us__ |       N/A       |       N/A       |
| transform vector2 x1       |   __2.7159 ns__ |       N/A       |     3.9917 ns   |    2.8070 ns   |    16.8257 ns   |       N/A       |       N/A       |       N/A       |
| transform vector2 x100     |     0.3463 us   |       N/A       |     0.4018 us   |  __0.2893 us__ |      1.709 us   |       N/A       |       N/A       |       N/A       |
| transform vector3 x1       |   __3.9868 ns__ |     5.5573 ns   |     8.4892 ns   |    4.4068 ns   |    25.0274 ns   |       N/A       |       N/A       |       N/A       |
| transform vector3 x100     |   __0.5905 us__ |     0.6584 us   |     0.8936 us   |    0.6365 us   |      2.513 us   |       N/A       |       N/A       |       N/A       |
| transform2 inverse         |       N/A       |       N/A       |     9.4094 ns   |    4.6388 ns   |       N/A       |   __3.9983 ns__ |       N/A       |       N/A       |
| transform2 mul transform2  |       N/A       |       N/A       |     9.8173 ns   |    6.2162 ns   |       N/A       |   __3.8699 ns__ |       N/A       |       N/A       |
| transform2 return self     |       N/A       |       N/A       |     4.8447 ns   |  __3.5091 ns__ |       N/A       |     4.1391 ns   |       N/A       |       N/A       |
| transform3 inverse         |       N/A       |       N/A       |    65.3982 ns   |   52.6160 ns   |       N/A       |  __32.0466 ns__ |       N/A       |       N/A       |
| transform3 mul transform3d |       N/A       |       N/A       |    10.9731 ns   |    9.9741 ns   |       N/A       |   __7.6754 ns__ |       N/A       |       N/A       |
| transform3 return self     |       N/A       |       N/A       |     7.1596 ns   |  __6.6096 ns__ |       N/A       |     7.0148 ns   |       N/A       |       N/A       |
| vector3 cross              |   __2.4542 ns__ |     3.5894 ns   |     3.2434 ns   |    3.4923 ns   |     3.5150 ns   |       N/A       |     3.2947 ns   |     7.1968 ns   |
| vector3 dot                |   __2.1001 ns__ |     2.3025 ns   |     2.2986 ns   |    2.3030 ns   |     2.3084 ns   |       N/A       |     2.3072 ns   |     3.7322 ns   |
| vector3 length             |   __2.1722 ns__ |   __2.1747 ns__ |     2.3414 ns   |  __2.1716 ns__ |   __2.2151 ns__ |       N/A       |   __2.2063 ns__ |     3.4787 ns   |
| vector3 normalize          |   __4.4248 ns__ |   __4.3266 ns__ |     8.1124 ns   |    8.0704 ns   |     8.0747 ns   |       N/A       |       N/A       |     8.0778 ns   |
| vector3 return self        |   __2.4642 ns__ |     2.9591 ns   |     2.9586 ns   |      N/A       |     2.9579 ns   |       N/A       |     2.9633 ns   |     2.9572 ns   |


## TODOS:

 - [X] `Quaternion` type and methods
 - [ ] `expm()`: Exponential matrix implementation
 - [X] Eigenvalues
 - [X] QR decomposition


