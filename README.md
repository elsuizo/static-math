# Static Math in Rust programming language

- This crate take advantage of the static arrays in Rust for fast operations in
stack memory.

- This crate could be used in an `no-std` environment.

- The determinant of the matrixs are evaluated "in-place" without loops and code
bifurcations

types:
 - Matrix from 2x2 to 6x6
 - Vector from 2 elements to 6

Methods:

 - `det()`: Determinant of the matrix
 - `inverse()`: Inverse of the matrix
 - `norm2`: norm of the matrix

## Benchmarks

Using the criterion crate:

https://github.com/bheisler/criterion.rs

this are the results for one operation(6x6 matrix inverse and printing):


```text
inverse 6x6             time:   [108.15 us 113.65 us 120.04 us]
                        change: [-71.614% -70.047% -68.452%] (p = 0.00 < 0.05)
                        Performance has improved.
```

you can look the bench here: [bench](benches/bench_inverse.rs)


The same Matrix and test but in Julia language:

```text
BenchmarkTools.Trial:
  memory estimate:  80.33 KiB
  allocs estimate:  1298
  --------------
  minimum time:     2.039 ms (0.00% GC)
  median time:      2.267 ms (0.00% GC)
  mean time:        2.348 ms (0.81% GC)
  maximum time:     18.015 ms (76.06% GC)
  --------------
  samples:          2118
  evals/sample:     1

```


## TODOS:

 - [ ] Eigenvalues and Eigenvectors
 - [ ] `expm()`: Exponential matrix implementation
 - [ ] QR decomposition
 - [ ] `Quaternion` type and methods
