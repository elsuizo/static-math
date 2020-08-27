//-------------------------------------------------------------------------
// @file matrix5x5.rs
//
// @date 06/02/20 20:39:15
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
use num::{Float, One, Zero, Num};
use std::fmt;
use std::ops::{Add, Mul, Sub};
use std::ops::{Deref, DerefMut, Index, IndexMut};

use crate::slices_methods::*;
use crate::traits::LinearAlgebra;
use crate::matrix4x4::M44;
use crate::vector5::V5;

//-------------------------------------------------------------------------
//                        code
//-------------------------------------------------------------------------
/// A static Matrix of 5x5 shape
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct M55<T>([[T; 5]; 5]);

impl<T: Float + std::iter::Sum> LinearAlgebra<T> for M55<T> {
    fn rows(&self) -> usize {
        self.0.len()
    }

    fn cols(&self) -> usize {
        self.rows()
    }

    fn transpose(&self) -> Self {
        let a_00 = self[(0, 0)];
        let a_01 = self[(0, 1)];
        let a_02 = self[(0, 2)];
        let a_03 = self[(0, 3)];
        let a_04 = self[(0, 4)];
        let a_10 = self[(1, 0)];
        let a_11 = self[(1, 1)];
        let a_12 = self[(1, 2)];
        let a_13 = self[(1, 3)];
        let a_14 = self[(1, 4)];
        let a_20 = self[(2, 0)];
        let a_21 = self[(2, 1)];
        let a_22 = self[(2, 2)];
        let a_23 = self[(2, 3)];
        let a_24 = self[(2, 4)];
        let a_30 = self[(3, 0)];
        let a_31 = self[(3, 1)];
        let a_32 = self[(3, 2)];
        let a_33 = self[(3, 3)];
        let a_34 = self[(3, 4)];
        let a_40 = self[(4, 0)];
        let a_41 = self[(4, 1)];
        let a_42 = self[(4, 2)];
        let a_43 = self[(4, 3)];
        let a_44 = self[(4, 4)];

        M55::new([
            [a_00, a_10, a_20, a_30, a_40],
            [a_01, a_11, a_21, a_31, a_41],
            [a_02, a_12, a_22, a_32, a_42],
            [a_03, a_13, a_23, a_33, a_43],
            [a_04, a_14, a_24, a_34, a_44],
        ])
    }

    fn trace(&self) -> T {
        self[(0, 0)] + self[(1, 1)] + self[(2, 2)] + self[(3, 3)] + self[(4, 4)]
    }

    fn norm2(&self) -> T {
        T::sqrt(
            self.iter()
                .flatten()
                .cloned()
                .map(|element| element * element)
                .sum(),
        )
    }

    fn det(&self) -> T {
        let a_00 = self[(0, 0)];
        let a_01 = self[(0, 1)];
        let a_02 = self[(0, 2)];
        let a_03 = self[(0, 3)];
        let a_04 = self[(0, 4)];
        let a_10 = self[(1, 0)];
        let a_11 = self[(1, 1)];
        let a_12 = self[(1, 2)];
        let a_13 = self[(1, 3)];
        let a_14 = self[(1, 4)];
        let a_20 = self[(2, 0)];
        let a_21 = self[(2, 1)];
        let a_22 = self[(2, 2)];
        let a_23 = self[(2, 3)];
        let a_24 = self[(2, 4)];
        let a_30 = self[(3, 0)];
        let a_31 = self[(3, 1)];
        let a_32 = self[(3, 2)];
        let a_33 = self[(3, 3)];
        let a_34 = self[(3, 4)];
        let a_40 = self[(4, 0)];
        let a_41 = self[(4, 1)];
        let a_42 = self[(4, 2)];
        let a_43 = self[(4, 3)];
        let a_44 = self[(4, 4)];

        a_00 * a_11 * a_22 * a_33 * a_44
            - a_00 * a_11 * a_22 * a_34 * a_43
            - a_00 * a_11 * a_23 * a_32 * a_44
            + a_00 * a_11 * a_23 * a_34 * a_42
            + a_00 * a_11 * a_24 * a_32 * a_43
            - a_00 * a_11 * a_24 * a_33 * a_42
            - a_00 * a_12 * a_21 * a_33 * a_44
            + a_00 * a_12 * a_21 * a_34 * a_43
            + a_00 * a_12 * a_23 * a_31 * a_44
            - a_00 * a_12 * a_23 * a_34 * a_41
            - a_00 * a_12 * a_24 * a_31 * a_43
            + a_00 * a_12 * a_24 * a_33 * a_41
            + a_00 * a_13 * a_21 * a_32 * a_44
            - a_00 * a_13 * a_21 * a_34 * a_42
            - a_00 * a_13 * a_22 * a_31 * a_44
            + a_00 * a_13 * a_22 * a_34 * a_41
            + a_00 * a_13 * a_24 * a_31 * a_42
            - a_00 * a_13 * a_24 * a_32 * a_41
            - a_00 * a_14 * a_21 * a_32 * a_43
            + a_00 * a_14 * a_21 * a_33 * a_42
            + a_00 * a_14 * a_22 * a_31 * a_43
            - a_00 * a_14 * a_22 * a_33 * a_41
            - a_00 * a_14 * a_23 * a_31 * a_42
            + a_00 * a_14 * a_23 * a_32 * a_41
            - a_01 * a_10 * a_22 * a_33 * a_44
            + a_01 * a_10 * a_22 * a_34 * a_43
            + a_01 * a_10 * a_23 * a_32 * a_44
            - a_01 * a_10 * a_23 * a_34 * a_42
            - a_01 * a_10 * a_24 * a_32 * a_43
            + a_01 * a_10 * a_24 * a_33 * a_42
            + a_01 * a_12 * a_20 * a_33 * a_44
            - a_01 * a_12 * a_20 * a_34 * a_43
            - a_01 * a_12 * a_23 * a_30 * a_44
            + a_01 * a_12 * a_23 * a_34 * a_40
            + a_01 * a_12 * a_24 * a_30 * a_43
            - a_01 * a_12 * a_24 * a_33 * a_40
            - a_01 * a_13 * a_20 * a_32 * a_44
            + a_01 * a_13 * a_20 * a_34 * a_42
            + a_01 * a_13 * a_22 * a_30 * a_44
            - a_01 * a_13 * a_22 * a_34 * a_40
            - a_01 * a_13 * a_24 * a_30 * a_42
            + a_01 * a_13 * a_24 * a_32 * a_40
            + a_01 * a_14 * a_20 * a_32 * a_43
            - a_01 * a_14 * a_20 * a_33 * a_42
            - a_01 * a_14 * a_22 * a_30 * a_43
            + a_01 * a_14 * a_22 * a_33 * a_40
            + a_01 * a_14 * a_23 * a_30 * a_42
            - a_01 * a_14 * a_23 * a_32 * a_40
            + a_02 * a_10 * a_21 * a_33 * a_44
            - a_02 * a_10 * a_21 * a_34 * a_43
            - a_02 * a_10 * a_23 * a_31 * a_44
            + a_02 * a_10 * a_23 * a_34 * a_41
            + a_02 * a_10 * a_24 * a_31 * a_43
            - a_02 * a_10 * a_24 * a_33 * a_41
            - a_02 * a_11 * a_20 * a_33 * a_44
            + a_02 * a_11 * a_20 * a_34 * a_43
            + a_02 * a_11 * a_23 * a_30 * a_44
            - a_02 * a_11 * a_23 * a_34 * a_40
            - a_02 * a_11 * a_24 * a_30 * a_43
            + a_02 * a_11 * a_24 * a_33 * a_40
            + a_02 * a_13 * a_20 * a_31 * a_44
            - a_02 * a_13 * a_20 * a_34 * a_41
            - a_02 * a_13 * a_21 * a_30 * a_44
            + a_02 * a_13 * a_21 * a_34 * a_40
            + a_02 * a_13 * a_24 * a_30 * a_41
            - a_02 * a_13 * a_24 * a_31 * a_40
            - a_02 * a_14 * a_20 * a_31 * a_43
            + a_02 * a_14 * a_20 * a_33 * a_41
            + a_02 * a_14 * a_21 * a_30 * a_43
            - a_02 * a_14 * a_21 * a_33 * a_40
            - a_02 * a_14 * a_23 * a_30 * a_41
            + a_02 * a_14 * a_23 * a_31 * a_40
            - a_03 * a_10 * a_21 * a_32 * a_44
            + a_03 * a_10 * a_21 * a_34 * a_42
            + a_03 * a_10 * a_22 * a_31 * a_44
            - a_03 * a_10 * a_22 * a_34 * a_41
            - a_03 * a_10 * a_24 * a_31 * a_42
            + a_03 * a_10 * a_24 * a_32 * a_41
            + a_03 * a_11 * a_20 * a_32 * a_44
            - a_03 * a_11 * a_20 * a_34 * a_42
            - a_03 * a_11 * a_22 * a_30 * a_44
            + a_03 * a_11 * a_22 * a_34 * a_40
            + a_03 * a_11 * a_24 * a_30 * a_42
            - a_03 * a_11 * a_24 * a_32 * a_40
            - a_03 * a_12 * a_20 * a_31 * a_44
            + a_03 * a_12 * a_20 * a_34 * a_41
            + a_03 * a_12 * a_21 * a_30 * a_44
            - a_03 * a_12 * a_21 * a_34 * a_40
            - a_03 * a_12 * a_24 * a_30 * a_41
            + a_03 * a_12 * a_24 * a_31 * a_40
            + a_03 * a_14 * a_20 * a_31 * a_42
            - a_03 * a_14 * a_20 * a_32 * a_41
            - a_03 * a_14 * a_21 * a_30 * a_42
            + a_03 * a_14 * a_21 * a_32 * a_40
            + a_03 * a_14 * a_22 * a_30 * a_41
            - a_03 * a_14 * a_22 * a_31 * a_40
            + a_04 * a_10 * a_21 * a_32 * a_43
            - a_04 * a_10 * a_21 * a_33 * a_42
            - a_04 * a_10 * a_22 * a_31 * a_43
            + a_04 * a_10 * a_22 * a_33 * a_41
            + a_04 * a_10 * a_23 * a_31 * a_42
            - a_04 * a_10 * a_23 * a_32 * a_41
            - a_04 * a_11 * a_20 * a_32 * a_43
            + a_04 * a_11 * a_20 * a_33 * a_42
            + a_04 * a_11 * a_22 * a_30 * a_43
            - a_04 * a_11 * a_22 * a_33 * a_40
            - a_04 * a_11 * a_23 * a_30 * a_42
            + a_04 * a_11 * a_23 * a_32 * a_40
            + a_04 * a_12 * a_20 * a_31 * a_43
            - a_04 * a_12 * a_20 * a_33 * a_41
            - a_04 * a_12 * a_21 * a_30 * a_43
            + a_04 * a_12 * a_21 * a_33 * a_40
            + a_04 * a_12 * a_23 * a_30 * a_41
            - a_04 * a_12 * a_23 * a_31 * a_40
            - a_04 * a_13 * a_20 * a_31 * a_42
            + a_04 * a_13 * a_20 * a_32 * a_41
            + a_04 * a_13 * a_21 * a_30 * a_42
            - a_04 * a_13 * a_21 * a_32 * a_40
            - a_04 * a_13 * a_22 * a_30 * a_41
            + a_04 * a_13 * a_22 * a_31 * a_40
    }
    ///
    /// Calculate the inverse of the Matrix6x6 via tha Adjoint Matrix:
    /// A^(-1) = 1/det Adj
    /// where Adj = Cofactor.Transpose()
    /// Cofactor = (-1)^(i+j) M(i, j).det()
    ///
    fn inverse(&self) -> Option<Self> {
        let one = T::one();
        let det = self.det();
        if det.abs() > T::epsilon() {
            let mut cofactors: M55<T> = M55::zeros();
            for i in 0..self.rows() {
                for j in 0..self.cols() {
                    let sign = if (i + j) % 2 == 0 {one} else {-one};
                    // transpose in place
                    cofactors[(j, i)] = sign * self.get_submatrix((i, j)).det();
                }
            }
            Some(cofactors * (T::one() / det))
        } else {
            None
        }
    }

    /// Calculate de QR factorization of the M55 via gram-schmidt
    /// orthogonalization process
    fn qr(&self) -> Option<(Self, Self)> {
        let det = self.det();
        if det.abs() > T::epsilon() {
            let cols = self.get_cols();
            let mut q: [V5<T>; 5] = *M55::zeros().get_cols();
            for i in 0..q.len() {
                let mut q_tilde = cols[i];
                for k in 0..i {
                    q_tilde -= q[k] * project_x_over_y(&*cols[i], &*q[k]);
                }
                normalize(&mut *q_tilde);
                q[i] = q_tilde;
            }
            let basis = V5::new([q[0], q[1], q[2], q[3], q[4]]);
            let q     = M55::new_from_vecs(basis);
            let r     = q.transpose() * (*self);
            Some((q, r))
        } else {
            None
        }
    }
}

impl<T> M55<T> {
    pub fn new(data_input: [[T; 5]; 5]) -> Self {
        Self(data_input)
    }
    pub fn rows(&self) -> usize {
        self.0.len()
    }
    pub fn cols(&self) -> usize {
        self.rows()
    }
}

impl<T: Num + Copy> Add for M55<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let a_00 = self[(0, 0)];
        let a_01 = self[(0, 1)];
        let a_02 = self[(0, 2)];
        let a_03 = self[(0, 3)];
        let a_04 = self[(0, 4)];
        let a_10 = self[(1, 0)];
        let a_11 = self[(1, 1)];
        let a_12 = self[(1, 2)];
        let a_13 = self[(1, 3)];
        let a_14 = self[(1, 4)];
        let a_20 = self[(2, 0)];
        let a_21 = self[(2, 1)];
        let a_22 = self[(2, 2)];
        let a_23 = self[(2, 3)];
        let a_24 = self[(2, 4)];
        let a_30 = self[(3, 0)];
        let a_31 = self[(3, 1)];
        let a_32 = self[(3, 2)];
        let a_33 = self[(3, 3)];
        let a_34 = self[(3, 4)];
        let a_40 = self[(4, 0)];
        let a_41 = self[(4, 1)];
        let a_42 = self[(4, 2)];
        let a_43 = self[(4, 3)];
        let a_44 = self[(4, 4)];

        let b_00 = rhs[(0, 0)];
        let b_01 = rhs[(0, 1)];
        let b_02 = rhs[(0, 2)];
        let b_03 = rhs[(0, 3)];
        let b_04 = rhs[(0, 4)];
        let b_10 = rhs[(1, 0)];
        let b_11 = rhs[(1, 1)];
        let b_12 = rhs[(1, 2)];
        let b_13 = rhs[(1, 3)];
        let b_14 = rhs[(1, 4)];
        let b_20 = rhs[(2, 0)];
        let b_21 = rhs[(2, 1)];
        let b_22 = rhs[(2, 2)];
        let b_23 = rhs[(2, 3)];
        let b_24 = rhs[(2, 4)];
        let b_30 = rhs[(3, 0)];
        let b_31 = rhs[(3, 1)];
        let b_32 = rhs[(3, 2)];
        let b_33 = rhs[(3, 3)];
        let b_34 = rhs[(3, 4)];
        let b_40 = rhs[(4, 0)];
        let b_41 = rhs[(4, 1)];
        let b_42 = rhs[(4, 2)];
        let b_43 = rhs[(4, 3)];
        let b_44 = rhs[(4, 4)];

        M55::new([
            [
                a_00 + b_00,
                a_01 + b_01,
                a_02 + b_02,
                a_03 + b_03,
                a_04 + b_04,
            ],
            [
                a_10 + b_10,
                a_11 + b_11,
                a_12 + b_12,
                a_13 + b_13,
                a_14 + b_14,
            ],
            [
                a_20 + b_20,
                a_21 + b_21,
                a_22 + b_22,
                a_23 + b_23,
                a_24 + b_24,
            ],
            [
                a_30 + b_30,
                a_31 + b_31,
                a_32 + b_32,
                a_33 + b_33,
                a_34 + b_34,
            ],
            [
                a_40 + b_40,
                a_41 + b_41,
                a_42 + b_42,
                a_43 + b_43,
                a_44 + b_44,
            ],
        ])
    }
}

// M55 - M55
impl<T: Num + Copy> Sub for M55<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let a_00 = self[(0, 0)];
        let a_01 = self[(0, 1)];
        let a_02 = self[(0, 2)];
        let a_03 = self[(0, 3)];
        let a_04 = self[(0, 4)];
        let a_10 = self[(1, 0)];
        let a_11 = self[(1, 1)];
        let a_12 = self[(1, 2)];
        let a_13 = self[(1, 3)];
        let a_14 = self[(1, 4)];
        let a_20 = self[(2, 0)];
        let a_21 = self[(2, 1)];
        let a_22 = self[(2, 2)];
        let a_23 = self[(2, 3)];
        let a_24 = self[(2, 4)];
        let a_30 = self[(3, 0)];
        let a_31 = self[(3, 1)];
        let a_32 = self[(3, 2)];
        let a_33 = self[(3, 3)];
        let a_34 = self[(3, 4)];
        let a_40 = self[(4, 0)];
        let a_41 = self[(4, 1)];
        let a_42 = self[(4, 2)];
        let a_43 = self[(4, 3)];
        let a_44 = self[(4, 4)];

        let b_00 = rhs[(0, 0)];
        let b_01 = rhs[(0, 1)];
        let b_02 = rhs[(0, 2)];
        let b_03 = rhs[(0, 3)];
        let b_04 = rhs[(0, 4)];
        let b_10 = rhs[(1, 0)];
        let b_11 = rhs[(1, 1)];
        let b_12 = rhs[(1, 2)];
        let b_13 = rhs[(1, 3)];
        let b_14 = rhs[(1, 4)];
        let b_20 = rhs[(2, 0)];
        let b_21 = rhs[(2, 1)];
        let b_22 = rhs[(2, 2)];
        let b_23 = rhs[(2, 3)];
        let b_24 = rhs[(2, 4)];
        let b_30 = rhs[(3, 0)];
        let b_31 = rhs[(3, 1)];
        let b_32 = rhs[(3, 2)];
        let b_33 = rhs[(3, 3)];
        let b_34 = rhs[(3, 4)];
        let b_40 = rhs[(4, 0)];
        let b_41 = rhs[(4, 1)];
        let b_42 = rhs[(4, 2)];
        let b_43 = rhs[(4, 3)];
        let b_44 = rhs[(4, 4)];

        M55::new([
            [
                a_00 - b_00,
                a_01 - b_01,
                a_02 - b_02,
                a_03 - b_03,
                a_04 - b_04,
            ],
            [
                a_10 - b_10,
                a_11 - b_11,
                a_12 - b_12,
                a_13 - b_13,
                a_14 - b_14,
            ],
            [
                a_20 - b_20,
                a_21 - b_21,
                a_22 - b_22,
                a_23 - b_23,
                a_24 - b_24,
            ],
            [
                a_30 - b_30,
                a_31 - b_31,
                a_32 - b_32,
                a_33 - b_33,
                a_34 - b_34,
            ],
            [
                a_40 - b_40,
                a_41 - b_41,
                a_42 - b_42,
                a_43 - b_43,
                a_44 - b_44,
            ],
        ])
    }
}

// M55 * constant
impl<T: Num + Copy> Mul<T> for M55<T> {
    type Output = M55<T>;

    fn mul(self, rhs: T) -> M55<T> {
        let a_00 = self[(0, 0)] * rhs;
        let a_01 = self[(0, 1)] * rhs;
        let a_02 = self[(0, 2)] * rhs;
        let a_03 = self[(0, 3)] * rhs;
        let a_04 = self[(0, 4)] * rhs;
        let a_10 = self[(1, 0)] * rhs;
        let a_11 = self[(1, 1)] * rhs;
        let a_12 = self[(1, 2)] * rhs;
        let a_13 = self[(1, 3)] * rhs;
        let a_14 = self[(1, 4)] * rhs;
        let a_20 = self[(2, 0)] * rhs;
        let a_21 = self[(2, 1)] * rhs;
        let a_22 = self[(2, 2)] * rhs;
        let a_23 = self[(2, 3)] * rhs;
        let a_24 = self[(2, 4)] * rhs;
        let a_30 = self[(3, 0)] * rhs;
        let a_31 = self[(3, 1)] * rhs;
        let a_32 = self[(3, 2)] * rhs;
        let a_33 = self[(3, 3)] * rhs;
        let a_34 = self[(3, 4)] * rhs;
        let a_40 = self[(4, 0)] * rhs;
        let a_41 = self[(4, 1)] * rhs;
        let a_42 = self[(4, 2)] * rhs;
        let a_43 = self[(4, 3)] * rhs;
        let a_44 = self[(4, 4)] * rhs;

        M55::new([
            [a_00, a_01, a_02, a_03, a_04],
            [a_10, a_11, a_12, a_13, a_14],
            [a_20, a_21, a_22, a_23, a_24],
            [a_30, a_31, a_32, a_33, a_34],
            [a_40, a_41, a_42, a_43, a_44],
        ])
    }
}


// M55 * V5
impl<T: Num + Copy> Mul<V5<T>> for M55<T> {
    type Output = V5<T>;

    fn mul(self, rhs: V5<T>) -> V5<T> {

        let a_00 = self[(0, 0)];
        let a_01 = self[(0, 1)];
        let a_02 = self[(0, 2)];
        let a_03 = self[(0, 3)];
        let a_04 = self[(0, 4)];
        let a_10 = self[(1, 0)];
        let a_11 = self[(1, 1)];
        let a_12 = self[(1, 2)];
        let a_13 = self[(1, 3)];
        let a_14 = self[(1, 4)];
        let a_20 = self[(2, 0)];
        let a_21 = self[(2, 1)];
        let a_22 = self[(2, 2)];
        let a_23 = self[(2, 3)];
        let a_24 = self[(2, 4)];
        let a_30 = self[(3, 0)];
        let a_31 = self[(3, 1)];
        let a_32 = self[(3, 2)];
        let a_33 = self[(3, 3)];
        let a_34 = self[(3, 4)];
        let a_40 = self[(4, 0)];
        let a_41 = self[(4, 1)];
        let a_42 = self[(4, 2)];
        let a_43 = self[(4, 3)];
        let a_44 = self[(4, 4)];

        let a = rhs[0];
        let b = rhs[1];
        let c = rhs[2];
        let d = rhs[3];
        let e = rhs[4];

        let v0 = a_00 * a + a_01 * b + a_02 * c + a_03 * d + a_04 * e;
        let v1 = a_10 * a + a_11 * b + a_12 * c + a_13 * d + a_14 * e;
        let v2 = a_20 * a + a_21 * b + a_22 * c + a_23 * d + a_24 * e;
        let v3 = a_30 * a + a_31 * b + a_32 * c + a_33 * d + a_34 * e;
        let v4 = a_40 * a + a_41 * b + a_42 * c + a_43 * d + a_44 * e;

        V5::new([v0, v1, v2, v3, v4])
    }

}

impl<T: Num + Copy> Mul for M55<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let a_00 = self[(0, 0)];
        let a_01 = self[(0, 1)];
        let a_02 = self[(0, 2)];
        let a_03 = self[(0, 3)];
        let a_04 = self[(0, 4)];
        let a_10 = self[(1, 0)];
        let a_11 = self[(1, 1)];
        let a_12 = self[(1, 2)];
        let a_13 = self[(1, 3)];
        let a_14 = self[(1, 4)];
        let a_20 = self[(2, 0)];
        let a_21 = self[(2, 1)];
        let a_22 = self[(2, 2)];
        let a_23 = self[(2, 3)];
        let a_24 = self[(2, 4)];
        let a_30 = self[(3, 0)];
        let a_31 = self[(3, 1)];
        let a_32 = self[(3, 2)];
        let a_33 = self[(3, 3)];
        let a_34 = self[(3, 4)];
        let a_40 = self[(4, 0)];
        let a_41 = self[(4, 1)];
        let a_42 = self[(4, 2)];
        let a_43 = self[(4, 3)];
        let a_44 = self[(4, 4)];

        let b_00 = rhs[(0, 0)];
        let b_01 = rhs[(0, 1)];
        let b_02 = rhs[(0, 2)];
        let b_03 = rhs[(0, 3)];
        let b_04 = rhs[(0, 4)];
        let b_10 = rhs[(1, 0)];
        let b_11 = rhs[(1, 1)];
        let b_12 = rhs[(1, 2)];
        let b_13 = rhs[(1, 3)];
        let b_14 = rhs[(1, 4)];
        let b_20 = rhs[(2, 0)];
        let b_21 = rhs[(2, 1)];
        let b_22 = rhs[(2, 2)];
        let b_23 = rhs[(2, 3)];
        let b_24 = rhs[(2, 4)];
        let b_30 = rhs[(3, 0)];
        let b_31 = rhs[(3, 1)];
        let b_32 = rhs[(3, 2)];
        let b_33 = rhs[(3, 3)];
        let b_34 = rhs[(3, 4)];
        let b_40 = rhs[(4, 0)];
        let b_41 = rhs[(4, 1)];
        let b_42 = rhs[(4, 2)];
        let b_43 = rhs[(4, 3)];
        let b_44 = rhs[(4, 4)];

        M55::new([
            [
                a_00 * b_00 + a_01 * b_10 + a_02 * b_20 + a_03 * b_30 + a_04 * b_40,
                a_00 * b_01 + a_01 * b_11 + a_02 * b_21 + a_03 * b_31 + a_04 * b_41,
                a_00 * b_02 + a_01 * b_12 + a_02 * b_22 + a_03 * b_32 + a_04 * b_42,
                a_00 * b_03 + a_01 * b_13 + a_02 * b_23 + a_03 * b_33 + a_04 * b_43,
                a_00 * b_04 + a_01 * b_14 + a_02 * b_24 + a_03 * b_34 + a_04 * b_44,
            ],
            [
                a_10 * b_00 + a_11 * b_10 + a_12 * b_20 + a_13 * b_30 + a_14 * b_40,
                a_10 * b_01 + a_11 * b_11 + a_12 * b_21 + a_13 * b_31 + a_14 * b_41,
                a_10 * b_02 + a_11 * b_12 + a_12 * b_22 + a_13 * b_32 + a_14 * b_42,
                a_10 * b_03 + a_11 * b_13 + a_12 * b_23 + a_13 * b_33 + a_14 * b_43,
                a_10 * b_04 + a_11 * b_14 + a_12 * b_24 + a_13 * b_34 + a_14 * b_44,
            ],
            [
                a_20 * b_00 + a_21 * b_10 + a_22 * b_20 + a_23 * b_30 + a_24 * b_40,
                a_20 * b_01 + a_21 * b_11 + a_22 * b_21 + a_23 * b_31 + a_24 * b_41,
                a_20 * b_02 + a_21 * b_12 + a_22 * b_22 + a_23 * b_32 + a_24 * b_42,
                a_20 * b_03 + a_21 * b_13 + a_22 * b_23 + a_23 * b_33 + a_24 * b_43,
                a_20 * b_04 + a_21 * b_14 + a_22 * b_24 + a_23 * b_34 + a_24 * b_44,
            ],
            [
                a_30 * b_00 + a_31 * b_10 + a_32 * b_20 + a_33 * b_30 + a_34 * b_40,
                a_30 * b_01 + a_31 * b_11 + a_32 * b_21 + a_33 * b_31 + a_34 * b_41,
                a_30 * b_02 + a_31 * b_12 + a_32 * b_22 + a_33 * b_32 + a_34 * b_42,
                a_30 * b_03 + a_31 * b_13 + a_32 * b_23 + a_33 * b_33 + a_34 * b_43,
                a_30 * b_04 + a_31 * b_14 + a_32 * b_24 + a_33 * b_34 + a_34 * b_44,
            ],
            [
                a_40 * b_00 + a_41 * b_10 + a_42 * b_20 + a_43 * b_30 + a_44 * b_40,
                a_40 * b_01 + a_41 * b_11 + a_42 * b_21 + a_43 * b_31 + a_44 * b_41,
                a_40 * b_02 + a_41 * b_12 + a_42 * b_22 + a_43 * b_32 + a_44 * b_42,
                a_40 * b_03 + a_41 * b_13 + a_42 * b_23 + a_43 * b_33 + a_44 * b_43,
                a_40 * b_04 + a_41 * b_14 + a_42 * b_24 + a_43 * b_34 + a_44 * b_44,
            ],
        ])
    }
}

impl<T: Num + Copy> Zero for M55<T> {
    fn zero() -> M55<T> {
        M55::new([[T::zero(); 5]; 5])
    }

    fn is_zero(&self) -> bool {
        *self == M55::zero()
    }
}

impl<T: Num + Copy> One for M55<T> {
    /// Create an identity matrix
    fn one() -> M55<T> {
        let one = T::one();
        let zero = T::zero();
        M55::new([
            [one, zero, zero, zero, zero],
            [zero, one, zero, zero, zero],
            [zero, zero, one, zero, zero],
            [zero, zero, zero, one, zero],
            [zero, zero, zero, zero, one],
        ])
    }
}

// NOTE(elsuizo:2020-04-26): poniendo ese Trait anda el norm2 funcional
impl<T: Num + Copy> M55<T> {
    /// contruct identity matrix
    pub fn identity() -> M55<T> {
        <M55<T> as One>::one()
    }

    /// construct the matrix with all zeros
    pub fn zeros() -> M55<T> {
        <M55<T> as Zero>::zero()
    }

    /// transform the matrix to a flatten vector
    pub fn as_vec(&self) -> [T; 25] {
        let mut result = [T::zero(); 25];
        for (index, element) in self.iter().flatten().enumerate() {
            result[index] = *element;
        }
        result
    }

    pub fn new_from_vecs(cols: V5<V5<T>>) -> Self {
        let mut result = Self::zeros();

        for i in 0..result.cols() {
            result[(i, 0)] = cols[0][i];
            result[(i, 1)] = cols[1][i];
            result[(i, 2)] = cols[2][i];
            result[(i, 3)] = cols[3][i];
            result[(i, 4)] = cols[4][i];
        }
        result
    }

    /// get the diagonal of the matrix
    pub fn get_diagonal(&self) -> V5<T> {
        let mut result = V5::zeros();
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

    /// get the a submatrix from discard row `i` and column `j`
    fn get_submatrix(&self, selected: (usize, usize)) -> M44<T> {
        let mut values: [T; 16] = [T::zero(); 16];
        let mut result: M44<T> = M44::zeros();
        let mut index: usize = 0;
        for i in 0..self.rows() {
            for j in 0..self.cols() {
                if !(i == selected.0 || j == selected.1) {
                    values[index] = self[(i, j)];
                    index += 1;
                }
            }
        }
        index = 0;
        for r in 0..result.rows() {
            for c in 0..result.cols() {
                result[(r, c)] = values[index];
                index += 1;
            }
        }
        return result;
    }

    pub fn get_rows(self) -> V5<V5<T>> {
        let mut r0 = V5::zeros();
        let mut r1 = V5::zeros();
        let mut r2 = V5::zeros();
        let mut r3 = V5::zeros();
        let mut r4 = V5::zeros();

        for j in 0..self.rows() {
            r0[j] = self[(0, j)];
            r1[j] = self[(1, j)];
            r2[j] = self[(2, j)];
            r3[j] = self[(3, j)];
            r4[j] = self[(4, j)];
        }
        V5::new([r0, r1, r2, r3, r4])
    }

    pub fn get_cols(self) -> V5<V5<T>> {
        let mut c0 = V5::zeros();
        let mut c1 = V5::zeros();
        let mut c2 = V5::zeros();
        let mut c3 = V5::zeros();
        let mut c4 = V5::zeros();

        for i in 0..self.cols() {
            c0[i] = self[(i, 0)];
            c1[i] = self[(i, 1)];
            c2[i] = self[(i, 2)];
            c3[i] = self[(i, 3)];
            c4[i] = self[(i, 4)];
        }
        V5::new([c0, c1, c2, c3, c4])
    }

    pub fn get_upper_triagular(&self) -> [T; 10] {
        let zero = T::zero();
        let mut result: [T; 10] = [zero; 10];
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

    pub fn get_lower_triangular(&self) -> [T; 10] {
        let zero = T::zero();
        let mut result: [T; 10] = [zero; 10];
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

}

impl<T> Deref for M55<T> {
    type Target = [[T; 5]; 5];
    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> DerefMut for M55<T> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T> From<[[T; 5]; 5]> for M55<T> {
    fn from(data: [[T; 5]; 5]) -> M55<T> {
        M55(data)
    }
}

impl<T> Index<(usize, usize)> for M55<T> {
    type Output = T;
    fn index(&self, index: (usize, usize)) -> &T {
        &self.0[index.0][index.1]
    }
}

impl<T> IndexMut<(usize, usize)> for M55<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut T {
        &mut self.0[index.0][index.1]
    }
}

//-------------------------------------------------------------------------
//                        Display for M55
//-------------------------------------------------------------------------
impl<T: Num + fmt::Display> fmt::Display for M55<T> {
    fn fmt(&self, dest: &mut fmt::Formatter) -> fmt::Result {
        println!("");
        write!(
            dest,
            "|{0:<7.2} {1:^7.2} {2:^7.2} {3:^7.2} {4:>7.2}|\n",
            self[(0, 0)],
            self[(0, 1)],
            self[(0, 2)],
            self[(0, 3)],
            self[(0, 4)]
        )?;
        write!(
            dest,
            "|{0:<7.2} {1:^7.2} {2:^7.2} {3:^7.2} {4:>7.2}|\n",
            self[(1, 0)],
            self[(1, 1)],
            self[(1, 2)],
            self[(1, 3)],
            self[(1, 4)]
        )?;
        write!(
            dest,
            "|{0:<7.2} {1:^7.2} {2:^7.2} {3:^7.2} {4:>7.2}|\n",
            self[(2, 0)],
            self[(2, 1)],
            self[(2, 2)],
            self[(2, 3)],
            self[(2, 4)]
        )?;
        write!(
            dest,
            "|{0:<7.2} {1:^7.2} {2:^7.2} {3:^7.2} {4:>7.2}|\n",
            self[(3, 0)],
            self[(3, 1)],
            self[(3, 2)],
            self[(3, 3)],
            self[(3, 4)]
        )?;
        write!(
            dest,
            "|{0:<7.2} {1:^7.2} {2:^7.2} {3:^7.2} {4:>7.2}|\n",
            self[(4, 0)],
            self[(4, 1)],
            self[(4, 2)],
            self[(4, 3)],
            self[(4, 4)]
        )
    }
}

//-------------------------------------------------------------------------
//                        macros
//-------------------------------------------------------------------------
#[macro_export]
macro_rules! m55_new {
    ($($first_row:expr),*;
     $($second_row:expr),*;
     $($third_row:expr),*;
     $($fourth_row:expr),*;
     $($fifth_row:expr),*
     ) => {
        M55::new([[$($first_row),*],
                 [$($second_row),*],
                 [$($third_row),*],
                 [$($fourth_row),*],
                 [$($fifth_row),*]])
    }
}

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------

#[cfg(test)]
mod test_matrix5x5 {

    use crate::matrix5x5::M55;
    use crate::traits::LinearAlgebra;
    use crate::utils::nearly_equal;
    use crate::utils::compare_vecs;
    use crate::vector5::V5;

    const EPS: f32 = 1e-6;

    #[test]
    fn matrix5x5_det_test() {

        use super::test_matrix5x5::EPS;

        let m = M55::new([
            [10.0, 1.0, 7.0, 1.0, 5.0],
            [2.0, 4.0, 8.0, 3.0, 2.0],
            [5.0, 1.0, 2.0, 9.0, 10.0],
            [6.0, 9.0, 9.0, 7.0, 3.0],
            [1.0, 8.0, 8.0, 10.0, 5.0],
        ]);
        let result = m.det();
        let expected = 49.99999999999798;
        assert!(nearly_equal(result, expected, EPS));
    }

    #[test]
    fn matrix5x5_sum_test() {

        use super::test_matrix5x5::EPS;

        let m = M55::new([
            [10.0, 1.0, 7.0, 1.0, 5.0],
            [2.0, 4.0, 8.0, 3.0, 2.0],
            [5.0, 1.0, 2.0, 9.0, 10.0],
            [6.0, 9.0, 9.0, 7.0, 3.0],
            [1.0, 8.0, 8.0, 10.0, 5.0],
        ]);

        let expected = M55::new([
            [20.0, 2.0, 14.0, 2.0, 10.0],
            [4.0, 8.0, 16.0, 6.0, 4.0],
            [10.0, 2.0, 4.0, 18.0, 20.0],
            [12.0, 18.0, 18.0, 14.0, 6.0],
            [2.0, 16.0, 16.0, 20.0, 10.0],
        ]);
        let result = m + m;

        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn sub_test() {
        use super::test_matrix5x5::EPS;

        let m1 = m55_new!(10.0, 1.0, 7.0,  1.0,  5.0;
                           2.0, 4.0, 8.0,  3.0,  2.0;
                           5.0, 1.0, 2.0,  9.0, 10.0;
                           6.0, 9.0, 9.0,  7.0,  3.0;
                           1.0, 8.0, 8.0, 10.0,  5.0);

        let m2 = m55_new!(10.0, 1.0, 7.0,  1.0,  5.0;
                           2.0, 4.0, 8.0,  3.0,  2.0;
                           5.0, 1.0, 2.0,  9.0, 10.0;
                           6.0, 9.0, 9.0,  7.0,  3.0;
                           1.0, 8.0, 8.0, 10.0,  5.0);

        let result = m1 - m2;
        let expected = m55_new!( 0.0, 0.0, 0.0,  0.0,  0.0;
                                 0.0, 0.0, 0.0,  0.0,  0.0;
                                 0.0, 0.0, 0.0,  0.0,  0.0;
                                 0.0, 0.0, 0.0,  0.0,  0.0;
                                 0.0, 0.0, 0.0,  0.0,  0.0);

        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn matrix5x5_product_test() {

        use super::test_matrix5x5::EPS;

        let m = M55::new([
            [10.0, 1.0, 7.0, 1.0, 5.0],
            [2.0, 4.0, 8.0, 3.0, 2.0],
            [5.0, 1.0, 2.0, 9.0, 10.0],
            [6.0, 9.0, 9.0, 7.0, 3.0],
            [1.0, 8.0, 8.0, 10.0, 5.0],
        ]);
        let result = m * m;
        let expected = M55::new([
            [148.0, 70.0, 141.0, 133.0, 150.0],
            [88.0, 69.0, 105.0, 127.0, 117.0],
            [126.0, 172.0, 208.0, 189.0, 124.0],
            [168.0, 138.0, 219.0, 193.0, 174.0],
            [131.0, 171.0, 217.0, 217.0, 156.0],
        ]);

        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }
    #[test]
    fn matrix5x5_norm2_test() {

        use super::test_matrix5x5::EPS;

        let m = M55::new([
            [10.0, 1.0, 7.0, 1.0, 5.0],
            [2.0, 4.0, 8.0, 3.0, 2.0],
            [5.0, 1.0, 2.0, 9.0, 10.0],
            [6.0, 9.0, 9.0, 7.0, 3.0],
            [1.0, 8.0, 8.0, 10.0, 5.0],
        ]);

        let result = m.norm2();
        let expected = 31.52776554086889;
        assert!(nearly_equal(result, expected, EPS));
    }

    #[test]
    fn matrix5x5_const_product_test() {

        use super::test_matrix5x5::EPS;

        let m = M55::new([
            [10.0, 1.0, 7.0, 1.0, 5.0],
            [2.0, 4.0, 8.0, 3.0, 2.0],
            [5.0, 1.0, 2.0, 9.0, 10.0],
            [6.0, 9.0, 9.0, 7.0, 3.0],
            [1.0, 8.0, 8.0, 10.0, 5.0],
        ]);

        let result = m * 0.5;
        let expected = M55::new([
            [5.0, 0.5, 3.5, 0.5, 2.5],
            [1.0, 2.0, 4.0, 1.5, 1.0],
            [2.5, 0.5, 1.0, 4.5, 5.0],
            [3.0, 4.5, 4.5, 3.5, 1.5],
            [0.5, 4.0, 4.0, 5.0, 2.5],
        ]);
        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }
    #[test]
    fn matrix5x5_inv_test() {

        use super::test_matrix5x5::EPS;

        let m = M55::new([
            [10.0, 1.0, 7.0, 1.0, 5.0],
            [2.0, 4.0, 8.0, 3.0, 2.0],
            [5.0, 1.0, 2.0, 9.0, 10.0],
            [6.0, 9.0, 9.0, 7.0, 3.0],
            [1.0, 8.0, 8.0, 10.0, 5.0],
        ]);
        let expected = M55::new([
            [-11.98, 15.64, 9.32, 10.34, -19.12],
            [33.62, -44.16, -26.08, -28.46, 53.28],
            [-9.36, 12.48, 7.24, 7.88, -14.84],
            [-37.2, 48.6, 28.8, 31.6, -58.8],
            [37.98, -49.64, -29.32, -32.34, 60.12],
        ]);

        if let Some(result) = m.inverse() {
            assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
        }
    }

    #[test]
    fn get_cols_test() {
        let m1 = m55_new!(10.0, 1.0, 7.0,  1.0,  5.0;
                           2.0, 4.0, 8.0,  3.0,  2.0;
                           5.0, 1.0, 2.0,  9.0, 10.0;
                           6.0, 9.0, 9.0,  7.0,  3.0;
                           1.0, 8.0, 8.0, 10.0,  5.0);

        let result = m1.get_cols();

        let expected0 = V5::new([10.0, 2.0, 5.0, 6.0, 1.0]);
        let expected1 = V5::new([1.0, 4.0, 1.0, 9.0, 8.0]);
        let expected2 = V5::new([7.0, 8.0, 2.0, 9.0, 8.0]);
        let expected3 = V5::new([1.0, 3.0, 9.0, 7.0, 10.0]);
        let expected4 = V5::new([5.0, 2.0, 10.0, 3.0, 5.0]);

        let expected = V5::new([expected0, expected1, expected2, expected3, expected4]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn mul_vector_rhs() {
        let m = m55_new!(10.0, 1.0, 7.0,  1.0,  5.0;
                          2.0, 4.0, 8.0,  3.0,  2.0;
                          5.0, 1.0, 2.0,  9.0, 10.0;
                          6.0, 9.0, 9.0,  7.0,  3.0;
                          1.0, 8.0, 8.0, 10.0,  5.0);

        let v = V5::new([0.0, 1.0, 2.0, 3.0, 4.0]);
        let result = m * v;
        let expected = V5::new([38.0, 37.0, 72.0, 60.0, 74.0]);

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
        let m1 = m55_new!(10.0, 1.0, 7.0,  1.0,  5.0;
                           2.0, 4.0, 8.0,  3.0,  2.0;
                           5.0, 1.0, 2.0,  9.0, 10.0;
                           6.0, 9.0, 9.0,  7.0,  3.0;
                           1.0, 8.0, 8.0, 10.0,  5.0);

        let result = m1.transpose().get_rows();

        let expected0 = V5::new([10.0, 2.0, 5.0, 6.0, 1.0]);
        let expected1 = V5::new([1.0, 4.0, 1.0, 9.0, 8.0]);
        let expected2 = V5::new([7.0, 8.0, 2.0, 9.0, 8.0]);
        let expected3 = V5::new([1.0, 3.0, 9.0, 7.0, 10.0]);
        let expected4 = V5::new([5.0, 2.0, 10.0, 3.0, 5.0]);

        let expected = V5::new([expected0, expected1, expected2, expected3, expected4]);
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
        let expected = m55_new!(10.0, 1.0, 7.0,  1.0,  5.0;
                                 2.0, 4.0, 8.0,  3.0,  2.0;
                                 5.0, 1.0, 2.0,  9.0, 10.0;
                                 6.0, 9.0, 9.0,  7.0,  3.0;
                                 1.0, 8.0, 8.0, 10.0,  5.0);

        let cols = expected.get_cols();

        let result = M55::new_from_vecs(cols);

        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn qr_test() {
        let expected = m55_new!(10.0, 1.0, 7.0,  1.0,  5.0;
                                 2.0, 1.0, 8.0,  3.0,  2.0;
                                 5.0, 1.0, 1.0,  9.0, 10.0;
                                 6.0, 9.0, 9.0,  1.0,  3.0;
                                 1.0, 8.0, 8.0, 10.0,  5.0);
        if let Some((q, r)) = expected.qr() {
            let result = q * r;
            assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
            assert!(nearly_equal(q.det().abs(), 1.0, EPS));
        }
    }

    #[test]
    fn get_diagonal() {
        let m = m55_new!(10.0, 1.0, 7.0,  1.0,  5.0;
                          2.0, 1.0, 8.0,  3.0,  2.0;
                          5.0, 1.0, 1.0,  9.0, 10.0;
                          6.0, 9.0, 9.0,  1.0,  3.0;
                          1.0, 8.0, 8.0, 10.0,  5.0);
        let result = m.get_diagonal();
        let expected = V5::new([10.0, 1.0, 1.0, 1.0, 5.0]);
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
        let m = m55_new!(10.0, 1.0, 7.0,  1.0,  5.0;
                          2.0, 1.0, 8.0,  3.0,  2.0;
                          5.0, 1.0, 1.0,  9.0, 10.0;
                          6.0, 9.0, 9.0,  1.0,  3.0;
                          1.0, 8.0, 8.0, 10.0,  5.0);
        let result = m.get_upper_triagular();
        let expected = [1.0, 7.0, 1.0, 5.0, 8.0, 3.0, 2.0, 9.0, 10.0, 3.0];
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
        let m = m55_new!(10.0, 1.0, 7.0,  1.0,  5.0;
                          2.0, 1.0, 8.0,  3.0,  2.0;
                          5.0, 1.0, 1.0,  9.0, 10.0;
                          6.0, 9.0, 9.0,  1.0,  3.0;
                          1.0, 8.0, 8.0, 10.0,  5.0);
        let result = m.get_lower_triangular();
        let expected = [2.0, 5.0, 1.0, 6.0, 9.0, 9.0, 1.0, 8.0, 8.0, 10.0];
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

}
