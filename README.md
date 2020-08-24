![MIT](https://img.shields.io/badge/license-MIT-blue.svg)
[![Documentation](https://docs.rs/static-math/badge.svg)](https://docs.rs/static-math)
[![crates.io](https://img.shields.io/crates/v/static-math.svg)](https://crates.io/crates/static-math)

# Static Math in Rust programming language

- This crate take advantage of the static arrays in Rust for fast operations in
stack memory.

- This crate could be used in an `no-std` environment.

- The determinant of the matrixs are evaluated "in-place" without loops and code
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

this are the results for one operation(6x6 matrix inverse):


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


The same Matrix and test but in Julia language:

```text
BenchmarkTools.Trial:
  memory estimate:  33.48 KiB
  allocs estimate:  455
  --------------
  minimum time:     1.536 ms (0.00% GC)
  median time:      1.566 ms (0.00% GC)
  mean time:        1.643 ms (0.62% GC)
  maximum time:     20.027 ms (78.89% GC)
  --------------
  samples:          3040
  evals/sample:     1
```

## TODOS:

 - [ ] `Quaternion` type and methods
 - [ ] `expm()`: Exponential matrix implementation
 - [X] Eigenvalues
 - [X] QR decomposition


