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


TODOS:

 - [ ] Eigenvalues and Eigenvectors
 - [ ] `expm()`: Exponential matrix implementation
 - [ ] QR decomposition
 - [ ] `Quaternion` type and methods
