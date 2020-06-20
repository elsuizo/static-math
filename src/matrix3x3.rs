//-------------------------------------------------------------------------
// @file matrix3x3.rs
//
// @date 06/02/20 18:41:39
// @author Martin Noblia
// @email mnoblia@disroot.org
//
// @brief
//
// @detail
//
//  Licence:
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
//--------------------------------------------------------------------------
use std::fmt;
use std::ops::{Add, Mul};
use std::ops::{Deref, DerefMut, Index, IndexMut};

use crate::traits::LinearAlgebra;
use num::{Float, One, Zero, Num};

use crate::vector3::*;
//-------------------------------------------------------------------------
//                        code
//-------------------------------------------------------------------------

/// A static matrix of 3x3 shape
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct M33<T>([[T; 3]; 3]);

impl<T: Float> LinearAlgebra<T> for M33<T> {
    fn rows(&self) -> usize {
        self.0.len()
    }

    fn cols(&self) -> usize {
        self.rows()
    }

    fn transpose(&self) -> M33<T> {
        M33::new([
            [self[(0, 0)], self[(1, 0)], self[(2, 0)]],
            [self[(0, 1)], self[(1, 1)], self[(2, 1)]],
            [self[(0, 2)], self[(1, 2)], self[(2, 2)]],
        ])
    }

    fn trace(&self) -> T {
        return self[(0, 0)] + self[(1, 1)] + self[(2, 2)];
    }

    fn norm2(&self) -> T {
        T::sqrt(
            self[(0, 0)] * self[(0, 0)]
                + self[(1, 0)] * self[(1, 0)]
                + self[(2, 0)] * self[(2, 0)]
                + self[(0, 1)] * self[(0, 1)]
                + self[(1, 1)] * self[(1, 1)]
                + self[(2, 1)] * self[(2, 1)]
                + self[(0, 2)] * self[(0, 2)]
                + self[(1, 2)] * self[(1, 2)]
                + self[(2, 2)] * self[(2, 2)],
        )
    }

    fn det(&self) -> T {
        self[(0, 0)] * (self[(1, 1)] * self[(2, 2)] - self[(2, 1)] * self[(1, 2)])
            - self[(0, 1)] * (self[(1, 0)] * self[(2, 2)] - self[(1, 2)] * self[(2, 0)])
            + self[(0, 2)] * (self[(1, 0)] * self[(2, 1)] - self[(1, 1)] * self[(2, 0)])
    }

    // TODO(elsuizo:2020-06-02): use here utils::nearly_equal()
    fn inverse(&self) -> Option<Self> {
        let det = self.det();
        if det.abs() > T::epsilon() {
            let invdet = T::one() / det;
            let mut res = M33::zero();
            res[(0, 0)] = (self[(1, 1)] * self[(2, 2)] - self[(2, 1)] * self[(1, 2)]) * invdet;
            res[(0, 1)] = (self[(0, 2)] * self[(2, 1)] - self[(0, 1)] * self[(2, 2)]) * invdet;
            res[(0, 2)] = (self[(0, 1)] * self[(1, 2)] - self[(0, 2)] * self[(1, 1)]) * invdet;
            res[(1, 0)] = (self[(1, 2)] * self[(2, 0)] - self[(1, 0)] * self[(2, 2)]) * invdet;
            res[(1, 1)] = (self[(0, 0)] * self[(2, 2)] - self[(0, 2)] * self[(2, 0)]) * invdet;
            res[(1, 2)] = (self[(1, 0)] * self[(0, 2)] - self[(0, 0)] * self[(1, 2)]) * invdet;
            res[(2, 0)] = (self[(1, 0)] * self[(2, 1)] - self[(2, 0)] * self[(1, 1)]) * invdet;
            res[(2, 1)] = (self[(2, 0)] * self[(0, 1)] - self[(0, 0)] * self[(2, 1)]) * invdet;
            res[(2, 2)] = (self[(0, 0)] * self[(1, 1)] - self[(1, 0)] * self[(0, 1)]) * invdet;
            Some(res)
        } else {
            None
        }
    }
}

impl<T> M33<T> {
    pub fn new(data_input: [[T; 3]; 3]) -> M33<T> {
        M33(data_input)
    }
    pub fn rows(&self) -> usize {
        self.0.len()
    }
    pub fn cols(&self) -> usize {
        self.rows()
    }
}

impl<T: Num + Copy> M33<T> {
    /// contruct identity matrix
    pub fn identity() -> M33<T> {
        <M33<T> as One>::one()
    }

    /// construct the matrix with all zeros
    pub fn zeros() -> M33<T> {
        <M33<T> as Zero>::zero()
    }

    /// transform the matrix to a flatten vector
    pub fn as_vec(&self) -> [T; 9] {
        let mut result = [T::zero(); 9];
        for i in 0..self.rows() {
            for j in 0..self.cols() {
                result[i] = self[(i, j)];
            }
        }
        result
    }
}

impl<T: Num + Copy> M33<T> {
    pub fn get_cols(self) -> V3<V3<T>> {
        let mut c0 = V3::zeros();
        let mut c1 = V3::zeros();
        let mut c2 = V3::zeros();

        for j in 0..self.cols() {
            c0[j] = self[(0, j)];
            c1[j] = self[(1, j)];
            c2[j] = self[(2, j)];
        }
        V3::new([c0, c1, c2])
    }

    pub fn get_rows(self) -> V3<V3<T>> {
        let mut r0 = V3::zeros();
        let mut r1 = V3::zeros();
        let mut r2 = V3::zeros();

        for i in 0..self.rows() {
            r0[i] = self[(i, 0)];
            r1[i] = self[(i, 1)];
            r2[i] = self[(i, 2)];
        }
        V3::new([r0, r1, r2])
    }
}

impl<T: Num + Copy> Add for M33<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        M33::new([
            [
                self[(0, 0)] + rhs[(0, 0)],
                self[(0, 1)] + rhs[(0, 1)],
                self[(0, 2)] + rhs[(0, 2)],
            ],
            [
                self[(1, 0)] + rhs[(1, 0)],
                self[(1, 1)] + rhs[(1, 1)],
                self[(1, 2)] + rhs[(1, 2)],
            ],
            [
                self[(2, 0)] + rhs[(2, 0)],
                self[(2, 1)] + rhs[(2, 1)],
                self[(2, 2)] + rhs[(2, 2)],
            ],
        ])
    }
}

// M3 * V3
impl<T: Num + Copy> Mul<V3<T>> for M33<T> {
    type Output = V3<T>;

    fn mul(self, rhs: V3<T>) -> V3<T> {

        let a_00 = self[(0, 0)];
        let a_01 = self[(0, 1)];
        let a_02 = self[(0, 2)];
        let a_10 = self[(1, 0)];
        let a_11 = self[(1, 1)];
        let a_12 = self[(1, 2)];
        let a_20 = self[(2, 0)];
        let a_21 = self[(2, 1)];
        let a_22 = self[(2, 2)];

        let v0 = rhs[0];
        let v1 = rhs[1];
        let v2 = rhs[2];
        V3::new([a_00 * v0 + a_01 * v1 + a_02 * v2,
                 a_10 * v0 + a_11 * v1 + a_12 * v2,
                 a_20 * v0 + a_21 * v1 + a_22 * v2])
    }
}

// M3 * constant
impl<T: Num + Copy> Mul<T> for M33<T> {
    type Output = M33<T>;

    fn mul(self, rhs: T) -> Self::Output {
        let a_00 = self[(0, 0)] * rhs;
        let a_01 = self[(0, 1)] * rhs;
        let a_02 = self[(0, 2)] * rhs;
        let a_10 = self[(1, 0)] * rhs;
        let a_11 = self[(1, 1)] * rhs;
        let a_12 = self[(1, 2)] * rhs;
        let a_20 = self[(2, 0)] * rhs;
        let a_21 = self[(2, 1)] * rhs;
        let a_22 = self[(2, 2)] * rhs;

        M33::new([[a_00, a_01, a_02], [a_10, a_11, a_12], [a_20, a_21, a_22]])
    }
}

// M3 * M3
impl<T: Num + Copy> Mul for M33<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let m00 =
            self[(0, 0)] * rhs[(0, 0)] + self[(0, 1)] * rhs[(1, 0)] + self[(0, 2)] * rhs[(2, 0)];
        let m01 =
            self[(0, 0)] * rhs[(0, 1)] + self[(0, 1)] * rhs[(1, 1)] + self[(0, 2)] * rhs[(2, 1)];
        let m02 =
            self[(0, 0)] * rhs[(0, 2)] + self[(0, 1)] * rhs[(1, 2)] + self[(0, 2)] * rhs[(2, 2)];

        let m10 =
            self[(1, 0)] * rhs[(0, 0)] + self[(1, 1)] * rhs[(1, 0)] + self[(1, 2)] * rhs[(2, 0)];
        let m11 =
            self[(1, 0)] * rhs[(0, 1)] + self[(1, 1)] * rhs[(1, 1)] + self[(1, 2)] * rhs[(2, 1)];
        let m12 =
            self[(1, 0)] * rhs[(0, 2)] + self[(1, 1)] * rhs[(1, 2)] + self[(1, 2)] * rhs[(2, 2)];

        let m20 =
            self[(2, 0)] * rhs[(0, 0)] + self[(2, 1)] * rhs[(1, 0)] + self[(2, 2)] * rhs[(2, 0)];
        let m21 =
            self[(2, 0)] * rhs[(0, 1)] + self[(2, 1)] * rhs[(1, 1)] + self[(2, 2)] * rhs[(2, 1)];
        let m22 =
            self[(2, 0)] * rhs[(0, 2)] + self[(2, 1)] * rhs[(1, 2)] + self[(2, 2)] * rhs[(2, 2)];

        M33::new([[m00, m01, m02], [m10, m11, m12], [m20, m21, m22]])
    }
}

impl<T: Num + Copy> Zero for M33<T> {
    fn zero() -> M33<T> {
        M33::new([[T::zero(); 3]; 3])
    }

    fn is_zero(&self) -> bool {
        *self == M33::zero()
    }
}

impl<T: Num + Copy> One for M33<T> {
    /// Create an identity matrix
    fn one() -> M33<T> {
        let one = T::one();
        let zero = T::zero();
        M33::new([[one, zero, zero], [zero, one, zero], [zero, zero, one]])
    }
}
//
impl<T> Deref for M33<T> {
    type Target = [[T; 3]; 3];
    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> DerefMut for M33<T> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T> From<[[T; 3]; 3]> for M33<T> {
    fn from(data: [[T; 3]; 3]) -> M33<T> {
        M33(data)
    }
}

impl<T> Index<(usize, usize)> for M33<T> {
    type Output = T;
    fn index(&self, index: (usize, usize)) -> &T {
        &self.0[index.0][index.1]
    }
}

impl<T> IndexMut<(usize, usize)> for M33<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut T {
        &mut self.0[index.0][index.1]
    }
}

//-------------------------------------------------------------------------
//                        Display for M33
//-------------------------------------------------------------------------
impl<T: Num + fmt::Display> fmt::Display for M33<T> {
    fn fmt(&self, dest: &mut fmt::Formatter) -> fmt::Result {
        println!("");
        write!(
            dest,
            "|{0:<7.2} {1:^7.2} {2:>7.2}|\n",
            self[(0, 0)],
            self[(0, 1)],
            self[(0, 2)]
        )?;
        write!(
            dest,
            "|{0:<7.2} {1:^7.2} {2:>7.2}|\n",
            self[(1, 0)],
            self[(1, 1)],
            self[(1, 2)]
        )?;
        write!(
            dest,
            "|{0:<7.2} {1:^7.2} {2:>7.2}|\n",
            self[(2, 0)],
            self[(2, 1)],
            self[(2, 2)]
        )
    }
}

//-------------------------------------------------------------------------
//                        macros
//-------------------------------------------------------------------------
#[macro_export]
macro_rules! m33_new {
    ($($first_row:expr),*;
     $($second_row:expr),*;
     $($third_row:expr),*
     ) => {
        M33::new([[$($first_row),*], [$($second_row),*], [$($third_row),*]])
    }
}

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------

#[cfg(test)]
mod test_matrix3x3 {
    use crate::traits::LinearAlgebra;
    use crate::matrix3x3::M33;
    use crate::utils::nearly_equal;
    use crate::utils::compare_vecs;
    use crate::vector3::*;

    const EPS: f32 = 1e-8;

    #[test]
    fn create_matrix() {
        let matrix = M33::new([[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]]);
        assert_eq!(matrix[(0, 0)], 0.0);
        assert_eq!(matrix[(0, 1)], 1.0);
        assert_eq!(matrix[(0, 2)], 2.0);
        assert_eq!(matrix[(1, 0)], 3.0);
        assert_eq!(matrix[(1, 1)], 4.0);
        assert_eq!(matrix[(1, 2)], 5.0);
        assert_eq!(matrix[(2, 0)], 6.0);
        assert_eq!(matrix[(2, 1)], 7.0);
        assert_eq!(matrix[(2, 2)], 8.0);
    }

    #[test]
    fn trace_test() {
        let matrix = M33::new([[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]]);
        assert_eq!(matrix.trace(), 12.0);
    }

    #[test]
    fn add_matrix() {
        use super::test_matrix3x3::EPS;
        let m1 = M33::new([[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]]);

        let m2 = M33::new([[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]]);

        let expected = M33::new([[0.0, 2.0, 4.0], [6.0, 8.0, 10.0], [12.0, 14.0, 16.0]]);
        let result = m1 + m2;
        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn test_identity_creation() {
        use super::test_matrix3x3::EPS;
        let identity: M33<f32> = M33::identity();

        let expected = M33::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        assert!(compare_vecs(&identity.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn test_zeros_creation() {
        use super::test_matrix3x3::EPS;
        let zero: M33<f32> = M33::zeros();

        let expected = M33::new([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
        assert!(compare_vecs(&zero.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn test_trace() {
        let m: M33<f32> = M33::identity();
        assert_eq!(m.trace(), 3.0);
    }

    #[test]
    fn test_norm2() {
        use super::test_matrix3x3::EPS;
        let m = M33::new([[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]]);
        assert!(nearly_equal(m.norm2(), 14.2828568570857, EPS));
    }

    #[test]
    fn determinant_test() {
        let m = M33::new([[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 9.0]]);
        let expected = -3.0;
        let result = m.det();

        assert!(nearly_equal(result, expected, EPS));
    }

    #[test]
    fn inverse_test() {
        use super::test_matrix3x3::EPS;
        let m = M33::new([[1.0, 0.0, 3.0], [2.0, 1.0, 6.0], [1.0, 0.0, 9.0]]);
        // NOTE(elsuizo:2019-09-25): hay que buscar una que tenga una inversa mas facil jasjdfjas
        let expected = M33::new([
            [1.5, 0.0, -0.5],
            [-2.0, 1.0, 0.0],
            [-0.16666666666666666, 0.0, 0.16666666666666666],
        ]);

        if let Some(result) = m.inverse() {
            assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
        }
    }

    #[test]
    fn test_mul_m33_v3() {
        let m = M33::new([[1.0,  1.0,  1.0],
                          [3.0,  2.0,  1.0],
                          [7.0,  3.0,  3.0],]);
        let v = V3::new([2.0, 7.0, 6.0]);
        let result = m * v;
        let expected = V3::new([15.0, 26.0, 53.0]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }
}
