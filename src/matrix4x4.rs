//-------------------------------------------------------------------------
// @file matrix4x4.rs
//
// @date 06/02/20 19:54:42
// @author Martin Noblia
// @email mnoblia@disroot.org
//
// @brief
//
// @detail
//
// Licence MIT:
// Copyright <2020> <Martin Noblia>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.  THE SOFTWARE IS PROVIDED
// "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
// LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//-------------------------------------------------------------------------
use core::fmt;
use core::ops::{Add, Mul, Div, Sub, AddAssign, SubAssign, Neg};
use core::ops::{Deref, DerefMut, Index, IndexMut};

use num::{Float, One, Zero, Num, Signed};
use crate::matrix3x3::*;
use crate::slices_methods::*;
use crate::traits::LinearAlgebra;
use crate::vector4::V4;
use crate::utils::nearly_zero;

//-------------------------------------------------------------------------
//                        code
//-------------------------------------------------------------------------
/// A static Matrix of 4x4 shape
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct M44<T>([[T; 4]; 4]);

impl<T> M44<T> {
    pub const fn new(data_input: [[T; 4]; 4]) -> M44<T> {
        M44(data_input)
    }

    pub const fn rows(&self) -> usize {
        self.0.len()
    }

    pub const fn cols(&self) -> usize {
        self.rows()
    }
}

impl<T: Float + core::iter::Sum> LinearAlgebra<T> for M44<T> {
    fn rows(&self) -> usize {
        self.0.len()
    }

    fn cols(&self) -> usize {
        self.rows()
    }

    #[inline(always)]
    fn transpose(&self) -> M44<T> {
        M44::new([
            [self[(0, 0)], self[(1, 0)], self[(2, 0)], self[(3, 0)]],
            [self[(0, 1)], self[(1, 1)], self[(2, 1)], self[(3, 1)]],
            [self[(0, 2)], self[(1, 2)], self[(2, 2)], self[(3, 2)]],
            [self[(0, 3)], self[(1, 3)], self[(2, 3)], self[(3, 3)]],
        ])
    }

    #[inline(always)]
    fn trace(&self) -> T {
        self[(0, 0)] + self[(1, 1)] + self[(2, 2)] + self[(3, 3)]
    }

    fn norm2(&self) -> T {
        let a1 = self[(0, 0)];
        let a2 = self[(0, 1)];
        let a3 = self[(0, 2)];
        let a4 = self[(0, 3)];
        let a5 = self[(1, 0)];
        let a6 = self[(1, 1)];
        let a7 = self[(1, 2)];
        let a8 = self[(1, 3)];
        let a9 = self[(2, 0)];
        let a10 = self[(2, 1)];
        let a11 = self[(2, 2)];
        let a12 = self[(2, 3)];
        let a13 = self[(3, 0)];
        let a14 = self[(3, 1)];
        let a15 = self[(3, 2)];
        let a16 = self[(3, 3)];

        T::sqrt(a1 * a1
                + a2 * a2
                + a3 * a3
                + a4 * a4
                + a5 * a5
                + a6 * a6
                + a7 * a7
                + a8 * a8
                + a9 * a9
                + a10 * a10
                + a11 * a11
                + a12 * a12
                + a13 * a13
                + a14 * a14
                + a15 * a15
                + a16 * a16,
        )
    }

    #[inline(always)]
    fn det(&self) -> T {
        let a1 = self[(0, 0)];
        let a2 = self[(0, 1)];
        let a3 = self[(0, 2)];
        let a4 = self[(0, 3)];
        let a5 = self[(1, 0)];
        let a6 = self[(1, 1)];
        let a7 = self[(1, 2)];
        let a8 = self[(1, 3)];
        let a9 = self[(2, 0)];
        let a10 = self[(2, 1)];
        let a11 = self[(2, 2)];
        let a12 = self[(2, 3)];
        let a13 = self[(3, 0)];
        let a14 = self[(3, 1)];
        let a15 = self[(3, 2)];
        let a16 = self[(3, 3)];

        a1 * a10 * a15 * a8
        - a1 * a10 * a16 * a7
        - a1 * a11 * a14 * a8
        + a1 * a11 * a16 * a6
        + a1 * a12 * a14 * a7
        - a1 * a12 * a15 * a6
        - a10 * a13 * a3 * a8
        + a10 * a13 * a4 * a7
        - a10 * a15 * a4 * a5
        + a10 * a16 * a3 * a5
        + a11 * a13 * a2 * a8
        - a11 * a13 * a4 * a6
        + a11 * a14 * a4 * a5
        - a11 * a16 * a2 * a5
        - a12 * a13 * a2 * a7
        + a12 * a13 * a3 * a6
        - a12 * a14 * a3 * a5
        + a12 * a15 * a2 * a5
        + a14 * a3 * a8 * a9
        - a14 * a4 * a7 * a9
        - a15 * a2 * a8 * a9
        + a15 * a4 * a6 * a9
        + a16 * a2 * a7 * a9
        - a16 * a3 * a6 * a9
    }

    /// Calculate the inverse
    fn inverse(&self) -> Option<Self> {
        let det = self.det();
        if !nearly_zero(det) {
            let det_recip = det.recip();
            let a1 = self.get_submatrix((0, 0)).det() * det_recip;
            let a2 = -self.get_submatrix((1, 0)).det() * det_recip;
            let a3 = self.get_submatrix((2, 0)).det() * det_recip;
            let a4 = -self.get_submatrix((3, 0)).det() * det_recip;

            let a5 = -self.get_submatrix((0, 1)).det() * det_recip;
            let a6 = self.get_submatrix((1, 1)).det() * det_recip;
            let a7 = -self.get_submatrix((2, 1)).det() * det_recip;
            let a8 = self.get_submatrix((3, 1)).det() * det_recip;

            let a9 = self.get_submatrix((0, 2)).det() * det_recip;
            let a10 = -self.get_submatrix((1, 2)).det() * det_recip;
            let a11 = self.get_submatrix((2, 2)).det() * det_recip;
            let a12 = -self.get_submatrix((3, 2)).det() * det_recip;

            let a13 = -self.get_submatrix((0, 3)).det() * det_recip;
            let a14 = self.get_submatrix((1, 3)).det() * det_recip;
            let a15 = -self.get_submatrix((2, 3)).det() * det_recip;
            let a16 = self.get_submatrix((3, 3)).det() * det_recip;

            Some(M44::new([
                [a1, a2 , a3 , a4],
                [a5, a6 , a7 , a8],
                [a9, a10, a11, a12],
                [a13, a14, a15, a16],
            ]))
        } else {
            None
        }
    }

    /// Calculate de QR factorization of the M44 via gram-schmidt
    /// orthogonalization process
    fn qr(&self) -> Option<(Self, Self)> {
        if !nearly_zero(self.det()) {
            let cols = self.get_cols();
            let mut q: [V4<T>; 4] = *M44::zeros().get_cols();
            for i in 0..q.len() {
                let mut q_tilde = cols[i];
                for elem in q.iter().take(i) {
                    q_tilde -= *elem * project_x_over_y(&*cols[i], &**elem);
                }
                normalize(&mut *q_tilde);
                q[i] = q_tilde;
            }
            let basis = V4::new_from(q[0], q[1], q[2], q[3]);
            let q     = M44::new_from_vecs(basis);
            let r     = q.transpose() * (*self);
            Some((q, r))
        } else {
            None
        }
    }

}

impl<T: Float> M44<T> {
    /// compute the LU factorization
    pub fn lu(&self) -> (Self, T, V4<usize>) {
        const N: usize = 4;
        let tiny = T::from(1e-40).unwrap();
        let mut indx: V4<usize> = V4::zeros();
        let mut lu = *self;
        let mut vv = V4::zeros();
        let mut d = T::one();
        let mut big = T::zero();
        for i in 0..N {
            for j in 0..N {
                let temp = T::abs(lu[(i, j)]);
                if temp > big {
                    big = temp;
                }
            }
            if big == T::zero() {
                panic!("the matrix should be non zero")
            }
            vv[i] = big.recip();
        }
        for k in 0..N {
            big = T::zero();
            let mut i_max = k;
            for i in k..N {
                let temp = vv[i] * T::abs(lu[(i, k)]);
                if temp > big {
                    big = temp;
                    i_max = i;
                }
            }
            if k != i_max {
                for j in 0..N {
                    let temp = lu[(i_max, j)];
                    lu[(i_max, j)] = lu[(k, j)];
                    lu[(k, j)] = temp;
                }
                d = -d;
                vv[i_max] = vv[k];
            }
            indx[k] = i_max;
            if lu[(k, k)] == T::zero() {
                lu[(k, k)] = tiny;
            }
            for i in (k + 1)..N {
                lu[(i, k)] = lu[(i, k)] / lu[(k, k)];
                for j in (k + 1)..N {
                    lu[(i, j)] = lu[(i, j)] - lu[(i, k)] * lu[(k, j)];
                }
            }
        }
        let det = d * lu[(0, 0)] * lu[(1, 1)] * lu[(2, 2)] * lu[(3, 3)];
        // return
        (lu, det, indx)
    }
}

impl<T: Num + Copy> M44<T> {
    /// contruct identity matrix
    pub fn identity() -> M44<T> {
        <M44<T> as One>::one()
    }

    /// construct the matrix with all zeros
    pub fn zeros() -> M44<T> {
        <M44<T> as Zero>::zero()
    }

    /// transform the matrix to a flatten vector
    pub fn as_vec(&self) -> [T; 16] {
        let mut result = [T::zero(); 16];
        for (index, element) in self.iter().flatten().enumerate() {
            result[index] = *element;
        }
        result
    }

    pub fn new_from_vecs(cols: V4<V4<T>>) -> Self {
        let mut result = Self::zeros();

        for i in 0..result.cols() {
            result[(i, 0)] = cols[0][i];
            result[(i, 1)] = cols[1][i];
            result[(i, 2)] = cols[2][i];
            result[(i, 3)] = cols[3][i];
        }
        result
    }

    /// get the diagonal of the matrix
    pub fn get_diagonal(&self) -> V4<T> {
        let mut result = V4::zeros();
        let mut index: usize = 0;
        for i in 0..self.rows() {
            for j in 0..self.cols() {
                if i == j {
                    result[index] = self[(i, j)];
                    index += 1;
                }
            }
        }
        result
    }

    // NOTE(elsuizo:2021-04-30): with this inline the performance its worse
    // #[inline(always)]
    pub fn get_submatrix(&self, selected: (usize, usize)) -> M33<T> {
        let mut values: [T; 9] = [T::zero(); 9];
        let mut result: M33<T> = M33::zeros();
        let mut index: usize = 0;
        for i in 0..4 {
            for j in 0..4 {
                if !(i == selected.0 || j == selected.1) {
                    // get the values from the M33
                    values[index] = self[(i, j)];
                    index += 1;
                }
            }
        }
        index = 0;
        for r in 0..3 {
            for c in 0..3 {
                result[(r, c)] = values[index];
                index += 1;
            }
        }
        result
    }

    pub fn get_rows(self) -> V4<V4<T>> {
        let mut r0 = V4::zeros();
        let mut r1 = V4::zeros();
        let mut r2 = V4::zeros();
        let mut r3 = V4::zeros();

        for j in 0..self.rows() {
            r0[j] = self[(0, j)];
            r1[j] = self[(1, j)];
            r2[j] = self[(2, j)];
            r3[j] = self[(3, j)];
        }
        V4::new([r0, r1, r2, r3])
    }

    pub fn get_cols(self) -> V4<V4<T>> {
        let mut c0 = V4::zeros();
        let mut c1 = V4::zeros();
        let mut c2 = V4::zeros();
        let mut c3 = V4::zeros();

        for i in 0..self.cols() {
            c0[i] = self[(i, 0)];
            c1[i] = self[(i, 1)];
            c2[i] = self[(i, 2)];
            c3[i] = self[(i, 3)];
        }
        V4::new([c0, c1, c2, c3])
    }

    pub fn get_upper_triagular(&self) -> [T; 6] {
        let zero = T::zero();
        let mut result: [T; 6] = [zero; 6];
        let mut index = 0;
        for i in 0..self.rows() {
            for j in 0..self.cols() {
                if i < j {
                    result[index] = self[(i, j)];
                    index += 1;
                }
            }
        }
        result
    }

    pub fn get_lower_triangular(&self) -> [T; 6] {
        let zero = T::zero();
        let mut result: [T; 6] = [zero; 6];
        let mut index = 0;
        for i in 0..self.rows() {
            for j in 0..self.cols() {
                if i > j {
                    result[index] = self[(i, j)];
                    index += 1;
                }
            }
        }
        result
    }

    /// Applies `f` of each element in the M44
    pub fn for_each(&self, f: impl Fn(T) -> T) -> Self {
        let mut result = Self::zeros();
        for i in 0..self.rows() {
            for j in 0..self.cols() {
                result[(i, j)] = f(self[(i, j)]);
            }
        }
        result
    }

}

impl<T: Num + Copy> Add for M44<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let a1 = self[(0, 0)];
        let a2 = self[(0, 1)];
        let a3 = self[(0, 2)];
        let a4 = self[(0, 3)];
        let a5 = self[(1, 0)];
        let a6 = self[(1, 1)];
        let a7 = self[(1, 2)];
        let a8 = self[(1, 3)];
        let a9 = self[(2, 0)];
        let a10 = self[(2, 1)];
        let a11 = self[(2, 2)];
        let a12 = self[(2, 3)];
        let a13 = self[(3, 0)];
        let a14 = self[(3, 1)];
        let a15 = self[(3, 2)];
        let a16 = self[(3, 3)];
        let b1 = rhs[(0, 0)];
        let b2 = rhs[(0, 1)];
        let b3 = rhs[(0, 2)];
        let b4 = rhs[(0, 3)];
        let b5 = rhs[(1, 0)];
        let b6 = rhs[(1, 1)];
        let b7 = rhs[(1, 2)];
        let b8 = rhs[(1, 3)];
        let b9 = rhs[(2, 0)];
        let b10 = rhs[(2, 1)];
        let b11 = rhs[(2, 2)];
        let b12 = rhs[(2, 3)];
        let b13 = rhs[(3, 0)];
        let b14 = rhs[(3, 1)];
        let b15 = rhs[(3, 2)];
        let b16 = rhs[(3, 3)];
        M44::new([
            [a1 + b1, a2 + b2, a3 + b3, a4 + b4],
            [a5 + b5, a6 + b6, a7 + b7, a8 + b8],
            [a9 + b9, a10 + b10, a11 + b11, a12 + b12],
            [a13 + b13, a14 + b14, a15 + b15, a16 + b16],
        ])
    }
}

// M44 += M44
impl<T: Num + Copy> AddAssign for M44<T> {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other
    }
}

// M44 - M44
impl<T: Num + Copy> Sub for M44<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let a1  = self[(0, 0)];
        let a2  = self[(0, 1)];
        let a3  = self[(0, 2)];
        let a4  = self[(0, 3)];
        let a5  = self[(1, 0)];
        let a6  = self[(1, 1)];
        let a7  = self[(1, 2)];
        let a8  = self[(1, 3)];
        let a9  = self[(2, 0)];
        let a10 = self[(2, 1)];
        let a11 = self[(2, 2)];
        let a12 = self[(2, 3)];
        let a13 = self[(3, 0)];
        let a14 = self[(3, 1)];
        let a15 = self[(3, 2)];
        let a16 = self[(3, 3)];

        let b1 = rhs[(0, 0)];
        let b2 = rhs[(0, 1)];
        let b3 = rhs[(0, 2)];
        let b4 = rhs[(0, 3)];
        let b5 = rhs[(1, 0)];
        let b6 = rhs[(1, 1)];
        let b7 = rhs[(1, 2)];
        let b8 = rhs[(1, 3)];
        let b9 = rhs[(2, 0)];
        let b10 = rhs[(2, 1)];
        let b11 = rhs[(2, 2)];
        let b12 = rhs[(2, 3)];
        let b13 = rhs[(3, 0)];
        let b14 = rhs[(3, 1)];
        let b15 = rhs[(3, 2)];
        let b16 = rhs[(3, 3)];
        M44::new([
            [a1 - b1, a2 - b2, a3 - b3, a4 - b4],
            [a5 - b5, a6 - b6, a7 - b7, a8 - b8],
            [a9 - b9, a10 - b10, a11 - b11, a12 - b12],
            [a13 - b13, a14 - b14, a15 - b15, a16 - b16],
        ])
    }
}

// M44 -= M44
impl<T: Num + Copy> SubAssign for M44<T> {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other
    }
}

// M44 * constant
impl<T: Num + Copy> Mul<T> for M44<T> {
    type Output = M44<T>;

    fn mul(self, rhs: T) -> M44<T> {
        let a1 = self[(0, 0)] * rhs;
        let a2 = self[(0, 1)] * rhs;
        let a3 = self[(0, 2)] * rhs;
        let a4 = self[(0, 3)] * rhs;
        let a5 = self[(1, 0)] * rhs;
        let a6 = self[(1, 1)] * rhs;
        let a7 = self[(1, 2)] * rhs;
        let a8 = self[(1, 3)] * rhs;
        let a9 = self[(2, 0)] * rhs;
        let a10 = self[(2, 1)] * rhs;
        let a11 = self[(2, 2)] * rhs;
        let a12 = self[(2, 3)] * rhs;
        let a13 = self[(3, 0)] * rhs;
        let a14 = self[(3, 1)] * rhs;
        let a15 = self[(3, 2)] * rhs;
        let a16 = self[(3, 3)] * rhs;

        M44::new([
            [a1, a2, a3, a4],
            [a5, a6, a7, a8],
            [a9, a10, a11, a12],
            [a13, a14, a15, a16],
        ])
    }
}

// M44 / constant
impl<T: Num + Copy> Div<T> for M44<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        let a1 = self[(0, 0)] / rhs;
        let a2 = self[(0, 1)] / rhs;
        let a3 = self[(0, 2)] / rhs;
        let a4 = self[(0, 3)] / rhs;
        let a5 = self[(1, 0)] / rhs;
        let a6 = self[(1, 1)] / rhs;
        let a7 = self[(1, 2)] / rhs;
        let a8 = self[(1, 3)] / rhs;
        let a9 = self[(2, 0)] / rhs;
        let a10 = self[(2, 1)] / rhs;
        let a11 = self[(2, 2)] / rhs;
        let a12 = self[(2, 3)] / rhs;
        let a13 = self[(3, 0)] / rhs;
        let a14 = self[(3, 1)] / rhs;
        let a15 = self[(3, 2)] / rhs;
        let a16 = self[(3, 3)] / rhs;

        M44::new([
            [a1, a2, a3, a4],
            [a5, a6, a7, a8],
            [a9, a10, a11, a12],
            [a13, a14, a15, a16],
        ])
    }
}

// M44 * M44
impl<T: Num + Copy> Mul for M44<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let a1 = self[(0, 0)];
        let a2 = self[(0, 1)];
        let a3 = self[(0, 2)];
        let a4 = self[(0, 3)];
        let a5 = self[(1, 0)];
        let a6 = self[(1, 1)];
        let a7 = self[(1, 2)];
        let a8 = self[(1, 3)];
        let a9 = self[(2, 0)];
        let a10 = self[(2, 1)];
        let a11 = self[(2, 2)];
        let a12 = self[(2, 3)];
        let a13 = self[(3, 0)];
        let a14 = self[(3, 1)];
        let a15 = self[(3, 2)];
        let a16 = self[(3, 3)];
        let b1 = rhs[(0, 0)];
        let b2 = rhs[(0, 1)];
        let b3 = rhs[(0, 2)];
        let b4 = rhs[(0, 3)];
        let b5 = rhs[(1, 0)];
        let b6 = rhs[(1, 1)];
        let b7 = rhs[(1, 2)];
        let b8 = rhs[(1, 3)];
        let b9 = rhs[(2, 0)];
        let b10 = rhs[(2, 1)];
        let b11 = rhs[(2, 2)];
        let b12 = rhs[(2, 3)];
        let b13 = rhs[(3, 0)];
        let b14 = rhs[(3, 1)];
        let b15 = rhs[(3, 2)];
        let b16 = rhs[(3, 3)];
        M44::new([
            [
                a1 * b1 + a2 * b5 + a3 * b9 + a4 * b13,
                a1 * b2 + a2 * b6 + a3 * b10 + a4 * b14,
                a1 * b3 + a2 * b7 + a3 * b11 + a4 * b15,
                a1 * b4 + a2 * b8 + a3 * b12 + a4 * b16,
            ],
            [
                a5 * b1 + a6 * b5 + a7 * b9 + a8 * b13,
                a5 * b2 + a6 * b6 + a7 * b10 + a8 * b14,
                a5 * b3 + a6 * b7 + a7 * b11 + a8 * b15,
                a5 * b4 + a6 * b8 + a7 * b12 + a8 * b16,
            ],
            [
                a10 * b5 + a11 * b9 + a12 * b13 + a9 * b1,
                a10 * b6 + a11 * b10 + a12 * b14 + a9 * b2,
                a10 * b7 + a11 * b11 + a12 * b15 + a9 * b3,
                a10 * b8 + a11 * b12 + a12 * b16 + a9 * b4,
            ],
            [
                a13 * b1 + a14 * b5 + a15 * b9 + a16 * b13,
                a13 * b2 + a14 * b6 + a15 * b10 + a16 * b14,
                a13 * b3 + a14 * b7 + a15 * b11 + a16 * b15,
                a13 * b4 + a14 * b8 + a15 * b12 + a16 * b16,
            ],
        ])
    }
}

// -M44
impl<T: Num + Copy + Signed> Neg for M44<T> {
    type Output = Self;

    fn neg(self) -> Self {
        let a1  = -self[(0, 0)];
        let a2  = -self[(0, 1)];
        let a3  = -self[(0, 2)];
        let a4  = -self[(0, 3)];
        let a5  = -self[(1, 0)];
        let a6  = -self[(1, 1)];
        let a7  = -self[(1, 2)];
        let a8  = -self[(1, 3)];
        let a9  = -self[(2, 0)];
        let a10 = -self[(2, 1)];
        let a11 = -self[(2, 2)];
        let a12 = -self[(2, 3)];
        let a13 = -self[(3, 0)];
        let a14 = -self[(3, 1)];
        let a15 = -self[(3, 2)];
        let a16 = -self[(3, 3)];

        M44::new([
            [a1, a2, a3, a4],
            [a5, a6, a7, a8],
            [a9, a10, a11, a12],
            [a13, a14, a15, a16],
        ])
    }
}
// M44 * V4
impl<T: Num + Copy> Mul<V4<T>> for M44<T> {
    type Output = V4<T>;

    #[inline]
    fn mul(self, rhs: V4<T>) -> V4<T> {
        let a_00 = self[(0, 0)];
        let a_01 = self[(0, 1)];
        let a_02 = self[(0, 2)];
        let a_03 = self[(0, 3)];
        let a_10 = self[(1, 0)];
        let a_11 = self[(1, 1)];
        let a_12 = self[(1, 2)];
        let a_13 = self[(1, 3)];
        let a_20 = self[(2, 0)];
        let a_21 = self[(2, 1)];
        let a_22 = self[(2, 2)];
        let a_23 = self[(2, 3)];
        let a_30 = self[(3, 0)];
        let a_31 = self[(3, 1)];
        let a_32 = self[(3, 2)];
        let a_33 = self[(3, 3)];

        let a = rhs[0];
        let b = rhs[1];
        let c = rhs[2];
        let d = rhs[3];

        V4::new_from(a_00 * a + a_01 * b + a_02 * c + a_03 * d,
                     a_10 * a + a_11 * b + a_12 * c + a_13 * d,
                     a_20 * a + a_21 * b + a_22 * c + a_23 * d,
                     a_30 * a + a_31 * b + a_32 * c + a_33 * d)
    }

}

impl<T: Num + Copy> Zero for M44<T> {
    fn zero() -> M44<T> {
        M44::new([[T::zero(); 4]; 4])
    }

    fn is_zero(&self) -> bool {
        *self == M44::zero()
    }
}

impl<T: Num + Copy> One for M44<T> {
    /// Create an identity matrix
    fn one() -> M44<T> {
        let one = T::one();
        let zero = T::zero();
        M44::new([
            [one, zero, zero, zero],
            [zero, one, zero, zero],
            [zero, zero, one, zero],
            [zero, zero, zero, one],
        ])
    }
}
//
impl<T> Deref for M44<T> {
    type Target = [[T; 4]; 4];
    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> DerefMut for M44<T> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T> From<[[T; 4]; 4]> for M44<T> {
    fn from(data: [[T; 4]; 4]) -> M44<T> {
        M44(data)
    }
}

//-------------------------------------------------------------------------
//                        index elements
//-------------------------------------------------------------------------
impl<T> Index<(usize, usize)> for M44<T> {
    type Output = T;
    #[inline(always)]
    fn index(&self, index: (usize, usize)) -> &T {
        &self.0[index.0][index.1]
    }
}

impl<T> IndexMut<(usize, usize)> for M44<T> {
    #[inline(always)]
    fn index_mut(&mut self, index: (usize, usize)) -> &mut T {
        &mut self.0[index.0][index.1]
    }
}

//-------------------------------------------------------------------------
//                        index rows
//-------------------------------------------------------------------------
impl<T> Index<usize> for M44<T> {
    type Output = [T; 4];
    #[inline(always)]
    fn index(&self, index: usize) -> &[T; 4] {
        &self.0[index]
    }
}

impl<T> IndexMut<usize> for M44<T> {
    #[inline(always)]
    fn index_mut(&mut self, index: usize) -> &mut [T; 4] {
        &mut self.0[index]
    }
}

// TODO(elsuizo:2020-03-26): hay que hacerlo mas "inteligente" para que cuando
// ponemos un numero de mas de 1 cifra no se rompa
//-------------------------------------------------------------------------
//                        Display
//-------------------------------------------------------------------------
impl<T: Num + fmt::Display> fmt::Display for M44<T> {
    fn fmt(&self, dest: &mut fmt::Formatter) -> fmt::Result {
        println!();
        writeln!(
            dest,
            "|{0:<7.2} {1:^7.2} {2:^7.2} {3:>7.2}|",
            self[(0, 0)],
            self[(0, 1)],
            self[(0, 2)],
            self[(0, 3)]
        )?;
        writeln!(
            dest,
            "|{0:<7.2} {1:^7.2} {2:^7.2} {3:>7.2}|",
            self[(1, 0)],
            self[(1, 1)],
            self[(1, 2)],
            self[(1, 3)]
        )?;
        writeln!(
            dest,
            "|{0:<7.2} {1:^7.2} {2:^7.2} {3:>7.2}|",
            self[(2, 0)],
            self[(2, 1)],
            self[(2, 2)],
            self[(2, 3)]
        )?;
        writeln!(
            dest,
            "|{0:<7.2} {1:^7.2} {2:^7.2} {3:>7.2}|",
            self[(3, 0)],
            self[(3, 1)],
            self[(3, 2)],
            self[(3, 3)]
        )
    }
}

//-------------------------------------------------------------------------
//                        macros
//-------------------------------------------------------------------------
#[macro_export]
macro_rules! m44_new {
    ($($first_row:expr),*;
     $($second_row:expr),*;
     $($third_row:expr),*;
     $($fourth_row:expr),*
     ) => {
        M44::new([[$($first_row),*],
                 [$($second_row),*],
                 [$($third_row),*],
                 [$($fourth_row),*]])
    }
}

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------


#[cfg(test)]
mod test_matrix4x4 {
    use crate::traits::LinearAlgebra;
    use crate::matrix3x3::M33;
    use crate::matrix4x4::M44;
    use crate::vector4::V4;
    use crate::utils::nearly_equal;
    use crate::utils::compare_vecs;

    const EPS: f32 = 1e-8;

    #[test]
    fn matrix4x4_create_matrix4x4_test() {
        let m = M44::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0, 16.0],
        ]);

        assert_eq!(m[(1, 1)], 6.0);
    }

    #[test]
    fn matrix4x4_identity_creation_test() {

        use super::test_matrix4x4::EPS;

        let expected = M44::new([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);
        let result: M44<f32> = M44::identity();
        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn matrix4x4_add_test() {

        use super::test_matrix4x4::EPS;

        let m1 = M44::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0, 16.0],
        ]);

        let m2 = M44::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0, 16.0],
        ]);

        let expected = M44::new([
            [2.0, 4.0, 6.0, 8.0],
            [10.0, 12.0, 14.0, 16.0],
            [18.0, 20.0, 22.0, 24.0],
            [26.0, 28.0, 30.0, 32.0],
        ]);
        let result = m1 + m2;
        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn mul_vector_rhs() {
        let m = m44_new!(1.0,  2.0,  3.0,  4.0;
                         5.0,  6.0,  7.0,  8.0;
                         9.0, 10.0, 11.0, 12.0;
                        13.0, 14.0, 15.0, 16.0);

        let v = V4::new([0.0, 1.0, 2.0, 3.0]);

        let result = m * v;
        let expected = V4::new([20.0, 44.0, 68.0, 92.0]);

        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn sub_test() {
        use super::test_matrix4x4::EPS;

        let m1 = m44_new!( 1.0,  2.0,  3.0,  4.0;
                           5.0,  6.0,  7.0,  8.0;
                           9.0, 10.0, 11.0, 12.0;
                          13.0, 14.0, 15.0, 16.0);

        let m2 = m44_new!( 1.0,  2.0,  3.0,  4.0;
                           5.0,  6.0,  7.0,  8.0;
                           9.0, 10.0, 11.0, 12.0;
                          13.0, 14.0, 15.0, 16.0);
        let result = m1 - m2;
        let expected = m44_new!( 0.0,  0.0,  0.0,  0.0;
                                 0.0,  0.0,  0.0,  0.0;
                                 0.0,  0.0,  0.0,  0.0;
                                 0.0,  0.0,  0.0,  0.0);
        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn matrix4x4_product_test() {

        use super::test_matrix4x4::EPS;

        let m1 = M44::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0, 16.0],
        ]);

        let m2 = M44::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0, 16.0],
        ]);

        let expected = M44::new([
            [90.0, 100.0, 110.0, 120.0],
            [202.0, 228.0, 254.0, 280.0],
            [314.0, 356.0, 398.0, 440.0],
            [426.0, 484.0, 542.0, 600.0],
        ]);
        let result = m1 * m2;
        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn matrix4x4_det_test() {
        let m1 = M44::new([
            [1.0, 2.0, 3.0, 1.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 0.0, 11.0, 12.0],
            [13.0, 1.0, 15.0, 16.0],
        ]);

        let expected = 168.0;
        let result = m1.det();
        assert_eq!(result, expected);
    }

    #[test]
    fn matrix4x4_norm_test() {

        use super::test_matrix4x4::EPS;

        let m1 = M44::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0, 16.0],
        ]);
        // NOTE(elsuizo:2019-08-08): el resultado lo calculo con Julia
        let expected = 38.67815921162743;
        let result = m1.norm2();
        assert!(nearly_equal(result, expected, EPS));
    }

    #[test]
    fn matrix4x4_transpose_test() {

        use super::test_matrix4x4::EPS;

        let m1 = M44::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0, 16.0],
        ]);

        let expected = M44::new([
            [1.0, 5.0, 9.0, 13.0],
            [2.0, 6.0, 10.0, 14.0],
            [3.0, 7.0, 11.0, 15.0],
            [4.0, 8.0, 12.0, 16.0],
        ]);
        let result = m1.transpose();
        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn matrix4x4_zeros_test() {

        use super::test_matrix4x4::EPS;

        let expected = M44::new([
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]);
        let result: M44<f32> = M44::zeros();
        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn matrix4x4_get_submatrix_test() {

        use super::test_matrix4x4::EPS;

        let m = M44::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0, 16.0],
        ]);

        let result1 = m.get_submatrix((0, 0));

        let expected1 = M33::new([[6.0, 7.0, 8.0], [10.0, 11.0, 12.0], [14.0, 15.0, 16.0]]);

        assert!(compare_vecs(&result1.as_vec(), &expected1.as_vec(), EPS));

        let result2 = m.get_submatrix((0, 1));

        let expected2 = M33::new([[5.0, 7.0, 8.0], [9.0, 11.0, 12.0], [13.0, 15.0, 16.0]]);

        assert!(compare_vecs(&result2.as_vec(), &expected2.as_vec(), EPS));

        let result3 = m.get_submatrix((0, 2));

        let expected3 = M33::new([[5.0, 6.0, 8.0], [9.0, 10.0, 12.0], [13.0, 14.0, 16.0]]);

        assert!(compare_vecs(&result3.as_vec(), &expected3.as_vec(), EPS));
    }

    #[test]
    fn matrix4x4_inverse_test() {

        use super::test_matrix4x4::EPS;

        let m = M44::new([
            [1.0, 1.0, 1.0, -1.0],
            [1.0, 1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0, 1.0],
        ]);

        let expected = M44::new([
            [1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, -1.0 / 4.0],
            [1.0 / 4.0, 1.0 / 4.0, -1.0 / 4.0, 1.0 / 4.0],
            [1.0 / 4.0, -1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0],
            [-1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0],
        ]);

        if let Some(result) = m.inverse() {
            assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
        }
    }

    #[test]
    fn inverse_fail() {
        let m = M44::new([
            [1.0, 1.0, 1.0, -1.0],
            [1.0, 1.0, 1.0, -1.0],
            [1.0, -1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0, 1.0],
        ]);

        let result = m.inverse();
        let expected = None;
        assert_eq!(result, expected);
    }

    #[test]
    fn get_cols_test() {
        let m = m44_new!( 0.0,  1.0,  2.0,  3.0;
                          4.0,  5.0,  6.0,  7.0;
                          8.0,  9.0, 10.0, 11.0;
                         12.0, 13.0, 14.0, 15.0);
        let result = m.get_cols();

        let expected0 = V4::new([0.0, 4.0, 8.0, 12.0]);
        let expected1 = V4::new([1.0, 5.0, 9.0, 13.0]);
        let expected2 = V4::new([2.0, 6.0, 10.0, 14.0]);
        let expected3 = V4::new([3.0, 7.0, 11.0, 15.0]);

        let expected = V4::new([expected0, expected1, expected2, expected3]);

        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );

    }

    #[test]
    fn get_rows_test() {
        let m = m44_new!( 0.0,  1.0,  2.0,  3.0;
                          4.0,  5.0,  6.0,  7.0;
                          8.0,  9.0, 10.0, 11.0;
                         12.0, 13.0, 14.0, 15.0);
        let result = m.get_rows();

        let expected0 = V4::new([0.0, 1.0, 2.0, 3.0]);
        let expected1 = V4::new([4.0, 5.0, 6.0, 7.0]);
        let expected2 = V4::new([8.0, 9.0, 10.0, 11.0]);
        let expected3 = V4::new([12.0, 13.0, 14.0, 15.0]);

        let expected = V4::new([expected0, expected1, expected2, expected3]);

        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );

    }

    #[test]
    fn new_from_vecs_test() {
        let expected = m44_new!( 0.0,  1.0,  2.0,  3.0;
                                 4.0,  5.0,  6.0,  7.0;
                                 8.0,  9.0, 10.0, 11.0;
                                12.0, 13.0, 14.0, 15.0);

        let cols = expected.get_cols();

        let result = M44::new_from_vecs(cols);

        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn qr_test() {
        let expected = m44_new!( 0.0,  1.0,  2.0,  3.0;
                                 4.0,  5.0,  6.0,  7.0;
                                 8.0,  9.0, 10.0, 11.0;
                                12.0, 13.0, 14.0, 15.0);
        if let Some((q, r)) = expected.qr() {
            let result = q * r;
            assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
            assert!(nearly_equal(q.det().abs(), 1.0, EPS));
        }
    }

    #[test]
    fn get_diagonal() {
        let m = m44_new!( 0.0,  1.0,  2.0,  3.0;
                          4.0,  5.0,  6.0,  7.0;
                          8.0,  9.0, 10.0, 11.0;
                         12.0, 13.0, 14.0, 15.0);
        let result = m.get_diagonal();
        let expected = V4::new([0.0, 5.0, 10.0, 15.0]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn get_upper_triangular_test() {
        let m = m44_new!( 0.0,  1.0,  2.0,  3.0;
                          4.0,  5.0,  6.0,  7.0;
                          8.0,  9.0, 10.0, 11.0;
                         12.0, 13.0, 14.0, 15.0);
        let result = m.get_upper_triagular();
        let expected = [1.0, 2.0, 3.0, 6.0, 7.0, 11.0];
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn get_lower_triangular_test() {
        let m = m44_new!( 0.0,  1.0,  2.0,  3.0;
                          4.0,  5.0,  6.0,  7.0;
                          8.0,  9.0, 10.0, 11.0;
                         12.0, 13.0, 14.0, 15.0);

        let result = m.get_lower_triangular();
        let expected = [4.0, 8.0, 9.0, 12.0, 13.0, 14.0];
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn for_each_test() {
        let m = m44_new!( 0.0,  1.0,  2.0,  3.0;
                          4.0,  5.0,  6.0,  7.0;
                          8.0,  9.0, 10.0, 11.0;
                         12.0, 13.0, 14.0, 15.0);
        let result = m.for_each(|element| element + 1.0);
        let expected = m44_new!( 1.0,  2.0,  3.0,  4.0;
                                 5.0,  6.0,  7.0,  8.0;
                                 9.0,  10.0, 11.0, 12.0;
                                13.0, 14.0, 15.0, 16.0);
        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }
}
