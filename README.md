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
inverse 6x6             time:   [9.9747 us 9.9786 us 9.9827 us]
                        change: [-0.0532% +0.3519% +1.0973%] (p = 0.44 > 0.05)
                        No change in performance detected.
Found 4 outliers among 100 measurements (4.00%)
  1 (1.00%) high mild
  3 (3.00%) high severe

inverse 4x4             time:   [118.78 ns 118.97 ns 119.17 ns]
                        change: [-0.3000% -0.0129% +0.2796%] (p = 0.93 > 0.05)
                        No change in performance detected.
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


