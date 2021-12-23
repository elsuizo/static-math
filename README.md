![MIT](https://img.shields.io/badge/license-MIT-blue.svg)
[![Documentation](https://docs.rs/static-math/badge.svg)](https://docs.rs/static-math)
[![crates.io](https://img.shields.io/crates/v/static-math.svg)](https://crates.io/crates/static-math)

# Static Math in Rust programming language

"*Simple things should be simple, complex things should be possible*" Alan Kay.

- This crate take advantage of the static arrays in Rust for fast operations in
stack memory.

- We use a tuple to indexing elements: `m[(i, j)]` allowing nice interface with the `match` feature of Rust

- You could index the rows of the matrix with simply: `m[i]`

- No `unsafe` code :ballot_box_with_check:

- Could be optimize more with the use of SIMD

- This crate could be used in an `no-std` environment.

   by enabling the feature `no-std`, for example in your `Cargo.toml`:

   ```toml
   [dependencies.static-math]
   default-features = false
   version = "0.2.0"
   features = ["no-std"]
   ```

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

 - We have implemented `Quaternion`s (and all the most used methods)
 - We have implemented `DualQuaternion`s (and all the most used methods in Robotics and graphics like *Screw Linear Interpolation*)
 - We have implemented in the `transformations.rs` module a wide variety of functions used in Robotics (which conforms to the screw theory)

## Benchmarks

Using the criterion crate:

https://github.com/bheisler/criterion.rs

run with: `cargo bench`

Others benches comparing the performance with others crates are in this repo: https://github.com/bitshifter/mathbench-rs

*NOTE*: this is the only crate that not have unsafe code

with the following results:


| benchmark                      |          glam   |        cgmath   |      nalgebra   |       euclid   |           vek   |    pathfinder   |   static-math   |   ultraviolet   |
|--------------------------------|-----------------|-----------------|-----------------|----------------|-----------------|-----------------|-----------------|-----------------|
| euler 2d x10000                |      16.23 us   |      16.13 us   |    __9.954 us__ |     16.18 us   |       16.2 us   |      10.42 us   |     __9.97 us__ |      16.17 us   |
| euler 3d x10000                |    __15.95 us__ |      32.11 us   |      32.13 us   |     32.13 us   |      32.13 us   |    __16.27 us__ |      32.16 us   |      32.11 us   |
| matrix2 determinant            |   __2.0386 ns__ |     2.0999 ns   |     2.1018 ns   |      N/A       |     2.0997 ns   |     2.0987 ns   |     2.0962 ns   |     2.1080 ns   |
| matrix2 inverse                |   __2.8226 ns__ |     8.4418 ns   |     7.6303 ns   |      N/A       |       N/A       |     3.3459 ns   |     9.4636 ns   |     5.8796 ns   |
| matrix2 mul matrix2            |   __2.6036 ns__ |     5.0007 ns   |     4.8172 ns   |      N/A       |     9.3814 ns   |   __2.5516 ns__ |     4.7274 ns   |     4.9428 ns   |
| matrix2 mul vector2 x1         |     2.4904 ns   |     2.6144 ns   |     2.8714 ns   |      N/A       |     4.2139 ns   |   __2.0839 ns__ |     2.8873 ns   |     2.6250 ns   |
| matrix2 mul vector2 x100       |   227.5271 ns   |   243.3579 ns   |   265.1698 ns   |      N/A       |   400.6940 ns   | __219.7127 ns__ |   267.8780 ns   |   243.9880 ns   |
| matrix2 return self            |   __2.4235 ns__ |     2.8841 ns   |     2.8756 ns   |      N/A       |     2.8754 ns   |   __2.4147 ns__ |     2.8717 ns   |     2.8697 ns   |
| matrix2 transpose              |   __2.2887 ns__ |     3.0645 ns   |     7.9154 ns   |      N/A       |     2.9635 ns   |       N/A       |     3.0637 ns   |     3.0652 ns   |
| matrix3 determinant            |     3.9129 ns   |   __3.8107 ns__ |   __3.8191 ns__ |      N/A       |   __3.8180 ns__ |       N/A       |   __3.8151 ns__ |     8.9368 ns   |
| matrix3 inverse                |    17.5373 ns   |    18.6931 ns   |  __12.3183 ns__ |      N/A       |       N/A       |       N/A       |    12.8195 ns   |    21.9098 ns   |
| matrix3 mul matrix3            |     9.9578 ns   |    13.3648 ns   |     7.8154 ns   |      N/A       |    35.5802 ns   |       N/A       |   __6.4938 ns__ |    10.0527 ns   |
| matrix3 mul vector3 x1         |     4.8090 ns   |     4.9339 ns   |   __4.5046 ns__ |      N/A       |    12.5518 ns   |       N/A       |     4.8002 ns   |     4.8118 ns   |
| matrix3 mul vector3 x100       |   __0.4836 us__ |   __0.4808 us__ |   __0.4755 us__ |      N/A       |      1.247 us   |       N/A       |   __0.4816 us__ |   __0.4755 us__ |
| matrix3 return self            |   __5.4421 ns__ |   __5.4469 ns__ |   __5.4526 ns__ |      N/A       |   __5.4656 ns__ |       N/A       |   __5.4718 ns__ |   __5.4043 ns__ |
| matrix3 transpose              |   __9.9567 ns__ |  __10.0794 ns__ |    10.9704 ns   |      N/A       |   __9.9257 ns__ |       N/A       |    10.7350 ns   |    10.5334 ns   |
| matrix4 determinant            |   __6.2050 ns__ |    11.1041 ns   |    69.2549 ns   |   17.1809 ns   |    18.5233 ns   |       N/A       |    16.5331 ns   |     8.2704 ns   |
| matrix4 inverse                |  __16.4386 ns__ |    47.0674 ns   |    71.8174 ns   |   64.1356 ns   |   284.3703 ns   |       N/A       |    52.6993 ns   |    41.1780 ns   |
| matrix4 mul matrix4            |   __7.7715 ns__ |    26.7308 ns   |     8.6500 ns   |   10.4414 ns   |    86.1501 ns   |       N/A       |    21.7985 ns   |    26.8056 ns   |
| matrix4 mul vector4 x1         |   __3.0303 ns__ |     7.7400 ns   |     3.4091 ns   |      N/A       |    21.0968 ns   |       N/A       |     6.2971 ns   |     6.2537 ns   |
| matrix4 mul vector4 x100       |   __0.6136 us__ |     0.9676 us   |    __0.627 us__ |      N/A       |      2.167 us   |       N/A       |     0.7893 us   |     0.8013 us   |
| matrix4 return self            |     7.1741 ns   |   __6.8838 ns__ |     7.5030 ns   |      N/A       |     7.0410 ns   |       N/A       |   __6.7768 ns__ |     6.9508 ns   |
| matrix4 transpose              |   __6.6826 ns__ |    12.4966 ns   |    15.3265 ns   |      N/A       |    12.6386 ns   |       N/A       |    15.2657 ns   |    12.3396 ns   |
| ray-sphere intersection x10000 |       56.2 us   |       55.7 us   |    __15.32 us__ |     55.45 us   |      56.02 us   |       N/A       |       N/A       |      50.94 us   |
| rotation3 inverse              |   __2.3113 ns__ |     3.1752 ns   |     3.3292 ns   |    3.3311 ns   |     3.1808 ns   |       N/A       |     8.7109 ns   |     3.6535 ns   |
| rotation3 mul rotation3        |   __3.6584 ns__ |     7.5255 ns   |     7.4808 ns   |    8.1393 ns   |    14.1636 ns   |       N/A       |     6.8044 ns   |     7.6386 ns   |
| rotation3 mul vector3 x1       |   __6.4950 ns__ |     7.6808 ns   |     7.5784 ns   |    7.5746 ns   |    18.2547 ns   |       N/A       |     7.2727 ns   |     8.9732 ns   |
| rotation3 mul vector3 x100     |   __0.6465 us__ |     0.7844 us   |     0.7573 us   |    0.7533 us   |      1.769 us   |       N/A       |     0.7317 us   |     0.9416 us   |
| rotation3 return self          |   __2.4928 ns__ |     2.8740 ns   |     2.8687 ns   |      N/A       |     2.8724 ns   |       N/A       |     4.7868 ns   |     2.8722 ns   |
| transform point2 x1            |     2.7854 ns   |     2.8878 ns   |     4.4207 ns   |    2.8667 ns   |    11.9427 ns   |   __2.3601 ns__ |       N/A       |     4.1770 ns   |
| transform point2 x100          |     0.3316 us   |     0.3574 us   |     0.4445 us   |  __0.3008 us__ |      1.212 us   |     0.3184 us   |       N/A       |     0.4332 us   |
| transform point3 x1            |   __2.9619 ns__ |    10.6812 ns   |     6.1037 ns   |    7.7051 ns   |    13.2607 ns   |     3.0934 ns   |       N/A       |     6.8419 ns   |
| transform point3 x100          |   __0.6095 us__ |       1.27 us   |     0.8064 us   |    0.7674 us   |      1.446 us   |   __0.6189 us__ |       N/A       |     0.8899 us   |
| transform vector2 x1           |   __2.4944 ns__ |       N/A       |     3.7174 ns   |    2.6273 ns   |    11.9424 ns   |       N/A       |       N/A       |     3.0458 ns   |
| transform vector2 x100         |     0.3125 us   |       N/A       |     0.3871 us   |  __0.2817 us__ |      1.213 us   |       N/A       |       N/A       |     0.3649 us   |
| transform vector3 x1           |   __2.8091 ns__ |     7.7343 ns   |     5.5064 ns   |    4.4810 ns   |    15.4097 ns   |       N/A       |       N/A       |     4.8819 ns   |
| transform vector3 x100         |   __0.6035 us__ |     0.9439 us   |     0.7573 us   |    0.6327 us   |       1.63 us   |       N/A       |       N/A       |     0.6703 us   |
| transform2 inverse             |   __9.0256 ns__ |       N/A       |    12.2614 ns   |    9.4803 ns   |       N/A       |   __8.9047 ns__ |       N/A       |       N/A       |
| transform2 mul transform2      |     4.5111 ns   |       N/A       |     8.1434 ns   |    5.8677 ns   |       N/A       |   __3.8513 ns__ |       N/A       |       N/A       |
| transform2 return self         |   __4.1707 ns__ |       N/A       |     5.4356 ns   |    4.2775 ns   |       N/A       |   __4.1117 ns__ |       N/A       |       N/A       |
| transform3 inverse             |  __10.9869 ns__ |       N/A       |    71.4437 ns   |   56.0136 ns   |       N/A       |    23.0392 ns   |       N/A       |       N/A       |
| transform3 mul transform3d     |   __6.5903 ns__ |       N/A       |     8.5673 ns   |   10.1802 ns   |       N/A       |     7.6587 ns   |       N/A       |       N/A       |
| transform3 return self         |   __7.1828 ns__ |       N/A       |   __7.2619 ns__ |  __7.2407 ns__ |       N/A       |   __7.3214 ns__ |       N/A       |       N/A       |
| vector3 cross                  |   __2.4257 ns__ |     3.6842 ns   |     3.7945 ns   |    3.6821 ns   |     3.8323 ns   |       N/A       |     3.8622 ns   |     3.6927 ns   |
| vector3 dot                    |   __2.1055 ns__ |     2.3179 ns   |     2.3174 ns   |    2.3190 ns   |     2.3195 ns   |       N/A       |     2.3204 ns   |     2.3160 ns   |
| vector3 length                 |   __2.5020 ns__ |   __2.5002 ns__ |     2.5986 ns   |  __2.5013 ns__ |   __2.5021 ns__ |       N/A       |   __2.5036 ns__ |   __2.5017 ns__ |
| vector3 normalize              |   __4.0454 ns__ |     5.8411 ns   |     8.4069 ns   |    8.0679 ns   |     8.8137 ns   |       N/A       |       N/A       |     5.8440 ns   |
| vector3 return self            |   __2.4087 ns__ |     3.1021 ns   |     3.1061 ns   |      N/A       |     3.1052 ns   |       N/A       |     3.1136 ns   |     3.1071 ns   |


## TODOS:

 - [X] `Quaternion` type and methods
 - [ ] `expm()`: Exponential matrix implementation
 - [X] Eigenvalues
 - [X] QR decomposition


