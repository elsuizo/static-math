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
inverse 6x6             time:   [19.912 us 20.047 us 20.193 us]
                        change: [-32.374% -30.094% -28.425%] (p = 0.00 < 0.05)
                        Performance has improved.
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

 - [ ] Eigenvalues and Eigenvectors
 - [ ] `expm()`: Exponential matrix implementation
 - [X] QR decomposition
 - [ ] `Quaternion` type and methods


