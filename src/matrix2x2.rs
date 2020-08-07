//-------------------------------------------------------------------------
// @file matrix2x2.rs
//
// @date 06/01/20 22:14:25
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
// imports
use std::fmt;
use std::ops::{Add, Mul, Sub};
use std::ops::{Deref, DerefMut, Index, IndexMut};

use num::{Float, Num, Zero, One};
use crate::vector2::*;
use crate::traits::LinearAlgebra;
use crate::slices_methods::*;

// code

/// A static Matrix of 2x2 shape
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct M22<T>([[T; 2]; 2]);

impl<T> M22<T> {

    pub fn new(data_input: [[T; 2]; 2]) -> Self {
        Self(data_input)
    }

    pub fn create(a: T, b: T, c: T, d: T) -> Self {
        Self::new([[a, b], [c, d]])
    }

    pub fn rows(&self) -> usize {
        self.0.len()
    }

    pub fn cols(&self) -> usize {
        self.rows()
    }
}

impl<T: Float + std::iter::Sum> LinearAlgebra<T> for M22<T> {
    fn rows(&self) -> usize {
        self.0.len()
    }

    fn cols(&self) -> usize {
        self.rows()
    }

    fn det(&self) -> T {
        let a = self[(0, 0)];
        let b = self[(0, 1)];
        let c = self[(1, 0)];
        let d = self[(1, 1)];
        (a * d) - (c * b)
    }

    fn transpose(&self) -> M22<T> {
        let a = self[(0, 0)];
        let b = self[(0, 1)];
        let c = self[(1, 0)];
        let d = self[(1, 1)];
        M22::new([[a, c], [b, d]])
    }

    fn trace(&self) -> T {
        self[(0, 0)] + self[(1, 1)]
    }

    fn norm2(&self) -> T {
        let a = self[(0, 0)];
        let b = self[(0, 1)];
        let c = self[(1, 0)];
        let d = self[(1, 1)];
        T::sqrt(a * a + b * b + c * c + d * d)
    }

    fn inverse(&self) -> Option<Self> {
        let a = self[(0, 0)];
        let b = self[(0, 1)];
        let c = self[(1, 0)];
        let d = self[(1, 1)];
        let det = self.det();
        if det.abs() > T::epsilon() {
            Some(M22::new([[d / det, -b / det], [-c / det, a / det]]))
        } else {
            None
        }
    }

    /// Calculate de QR factorization of the M22 via gram-schmidt
    /// orthogonalization process
    fn qr(&self) -> Option<(Self, Self)> {
        let det = self.det();
        if det.abs() > T::epsilon() {
            let cols = self.get_cols();
            let mut q: [V2<T>; 2] = *M22::zeros().get_cols();
            for i in 0..q.len() {
                let mut q_tilde = cols[i];
                for k in 0..i {
                    q_tilde -= q[k] * project_x_over_y(&*cols[i], &*q[k]);
                }
                normalize(&mut *q_tilde);
                q[i] = q_tilde;
            }
            let basis = V2::new([q[0], q[1]]);
            let q     = M22::new_from_vecs(basis);
            let r     = q.transpose() * (*self);
            Some((q, r))
        } else {
            None
        }
    }

}

impl<T: Num + Copy> M22<T> {
    /// contruct identity matrix
    pub fn identity() -> M22<T> {
        <M22<T> as One>::one()
    }

    /// construct the matrix with all zeros
    pub fn zeros() -> M22<T> {
        <M22<T> as Zero>::zero()
    }

    /// transform the matrix to a flatten vector
    pub fn as_vec(&self) -> [T; 4] {
        let mut result = [T::zero(); 4];
        for i in 0..self.rows() {
            for j in 0..self.cols() {
                result[i] = self[(i, j)];
            }
        }
        result
    }

    /// construct the matrix from columns-vectors
    pub fn new_from_vecs(cols: V2<V2<T>>) -> Self {
        let mut result = Self::zeros();

        for i in 0..result.cols() {
            result[(i, 0)] = cols[0][i];
            result[(i, 1)] = cols[1][i];
        }
        result
    }

    /// get the diagonal of the matrix
    pub fn get_diagonal(&self) -> V2<T> {
        let mut result = V2::zeros();
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
}

// M22 * V2
impl<T: Num + Copy> Mul<V2<T>> for M22<T> {
    type Output = V2<T>;

    fn mul(self, rhs: V2<T>) -> V2<T> {
        let a1 = self[(0, 0)];
        let b1 = self[(0, 1)];
        let c1 = self[(1, 0)];
        let d1 = self[(1, 1)];

        let v1 = rhs[0];
        let v2 = rhs[1];
        V2::new([a1 * v1 + b1 * v2, c1 * v1 + d1 * v2])
    }
}

// M22 + M22
impl<T: Num + Copy> Add for M22<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let a1 = self[(0, 0)];
        let b1 = self[(0, 1)];
        let c1 = self[(1, 0)];
        let d1 = self[(1, 1)];

        let a2 = rhs[(0, 0)];
        let b2 = rhs[(0, 1)];
        let c2 = rhs[(1, 0)];
        let d2 = rhs[(1, 1)];
        M22::new([[a1 + a2, b1 + b2], [c1 + c2, d1 + d2]])
    }
}

// M22 - M22
impl<T: Num + Copy> Sub for M22<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let a1 = self[(0, 0)];
        let b1 = self[(0, 1)];
        let c1 = self[(1, 0)];
        let d1 = self[(1, 1)];

        let a2 = rhs[(0, 0)];
        let b2 = rhs[(0, 1)];
        let c2 = rhs[(1, 0)];
        let d2 = rhs[(1, 1)];
        M22::new([[a1 - a2, b1 - b2], [c1 - c2, d1 - d2]])
    }
}

impl<T: Num + Copy> M22<T> {
    /// get the rows of the matrix as a vectors
    pub fn get_rows(self) -> V2<V2<T>> {
        let mut r0 = V2::zeros();
        let mut r1 = V2::zeros();

        for j in 0..self.rows() {
            r0[j] = self[(0, j)];
            r1[j] = self[(1, j)]
        }

        V2::new([r0, r1])
    }

    /// get the columns of the matrix as a vectors
    pub fn get_cols(self) -> V2<V2<T>> {
        let mut c0 = V2::zeros();
        let mut c1 = V2::zeros();

        for i in 0..self.cols() {
            c0[i] = self[(i, 0)];
            c1[i] = self[(i, 1)]
        }

        V2::new([c0, c1])
    }
}

// NOTE(elsuizo:2020-06-10): maybe an error here is better
impl<T: Float + std::iter::Sum> M22<T> {
    /// calculate the real eigen values for the matrix
    pub fn real_eigenvals(&self) -> Option<V2<T>> {
        let tau = self.trace();
        let delta = self.det();
        let tau_2 = tau * tau;
        let four = T::from(4)?;
        let discr = tau_2 - four * delta;
        if discr < T::zero() {
            None
        } else {
            let two = T::from(2)?;
            let lambda2 = (tau - T::sqrt(discr)) / two;
            let lambda1 = (tau + T::sqrt(discr)) / two;
            Some(V2::new([lambda1, lambda2]))
        }
    }
}

// FIXME(elsuizo:2020-06-19): this is a hack
// f32 * M22<f32>
impl Mul<M22<f32>> for f32 {
    type Output = M22<f32>;

    fn mul(self, rhs: M22<f32>) -> M22<f32> {
        let a_00 = rhs[(0, 0)] * self;
        let a_01 = rhs[(0, 1)] * self;
        let a_10 = rhs[(1, 0)] * self;
        let a_11 = rhs[(1, 1)] * self;

        M22::new([[a_00, a_01], [a_10, a_11]])
    }
}

// M22 * constant
impl<T: Num + Copy> Mul<T> for M22<T> {
    type Output = M22<T>;

    fn mul(self, rhs: T) -> M22<T> {
        let a_00 = self[(0, 0)] * rhs;
        let a_01 = self[(0, 1)] * rhs;
        let a_10 = self[(1, 0)] * rhs;
        let a_11 = self[(1, 1)] * rhs;

        M22::new([[a_00, a_01], [a_10, a_11]])
    }
}

// M22 * M22
impl<T: Num + Copy> Mul for M22<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let a1 = self[(0, 0)];
        let b1 = self[(0, 1)];
        let c1 = self[(1, 0)];
        let d1 = self[(1, 1)];

        let a2 = rhs[(0, 0)];
        let b2 = rhs[(0, 1)];
        let c2 = rhs[(1, 0)];
        let d2 = rhs[(1, 1)];

        let m00 = a1 * a2 + b1 * c2;
        let m01 = a1 * b2 + b1 * d2;

        let m10 = c1 * a2 + d1 * c2;
        let m11 = c1 * b2 + d1 * d2;
        M22::new([[m00, m01], [m10, m11]])
    }
}

impl<T: Num + Copy> Zero for M22<T> {
    fn zero() -> M22<T> {
        M22::new([[T::zero(); 2]; 2])
    }

    fn is_zero(&self) -> bool {
        *self == M22::zero()
    }
}

impl<T: Num + Copy> One for M22<T> {
    /// Create an identity matrix
    fn one() -> M22<T> {
        let one = T::one();
        let zero = T::zero();
        M22::new([[one, zero], [zero, one]])
    }
}

impl<T> Deref for M22<T> {
    type Target = [[T; 2]; 2];
    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> DerefMut for M22<T> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T> From<[[T; 2]; 2]> for M22<T> {
    fn from(data: [[T; 2]; 2]) -> M22<T> {
        M22(data)
    }
}

impl<T> Index<(usize, usize)> for M22<T> {
    type Output = T;
    fn index(&self, index: (usize, usize)) -> &T {
        &self.0[index.0][index.1]
    }
}

impl<T> IndexMut<(usize, usize)> for M22<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut T {
        &mut self.0[index.0][index.1]
    }
}

//-------------------------------------------------------------------------
//                        macros
//-------------------------------------------------------------------------
#[macro_export]
macro_rules! m22_new {
    ($($first_row:expr),* ; $($second_row:expr),*) => {
        M22::new([[$($first_row),*], [$($second_row),*]])
    }
}

//-------------------------------------------------------------------------
//                        Display for M22
//-------------------------------------------------------------------------
impl<T: Num + fmt::Display> fmt::Display for M22<T> {
    fn fmt(&self, dest: &mut fmt::Formatter) -> fmt::Result {
        println!("");
        write!(dest, "|{0:^3.2} {1:^3.2}|\n", self[(0, 0)], self[(0, 1)])?;
        write!(dest, "|{0:^3.2} {1:^3.2}|\n", self[(1, 0)], self[(1, 1)])
    }
}

//-------------------------------------------------------------------------
//                        testing
//-------------------------------------------------------------------------

#[cfg(test)]
mod test_matrix2x2 {
    use crate::traits::LinearAlgebra;
    use crate::matrix2x2::M22;
    use crate::utils::{compare_vecs, nearly_equal};
    use crate::vector2::V2;

    const EPS: f32 = 1e-7;

    #[test]
    fn create_m22_floats() {
        let matrix = M22::new([[0.0, 1.0], [2.0, 3.0]]);
        assert_eq!(matrix[(0, 0)], 0.0);
        assert_eq!(matrix[(0, 1)], 1.0);
        assert_eq!(matrix[(1, 0)], 2.0);
        assert_eq!(matrix[(1, 1)], 3.0);
    }

    #[test]
    fn create_m22_test() {
        let m = m22_new!(0.0, 1.0;
                         2.0, 3.0);

        assert_eq!(m[(0, 0)], 0.0);
        assert_eq!(m[(0, 1)], 1.0);
        assert_eq!(m[(1, 0)], 2.0);
        assert_eq!(m[(1, 1)], 3.0);
    }

    #[test]
    fn create_m22_ints() {
        let m = M22::new([[0, 1], [2, 3]]);
        assert_eq!(m[(0, 0)], 0);
        assert_eq!(m[(0, 1)], 1);
        assert_eq!(m[(1, 0)], 2);
        assert_eq!(m[(1, 1)], 3);
    }

    #[test]
    fn create_identity_floats() {
        let expected = M22::new([[1.0, 0.0], [0.0, 1.0]]);
        let result: M22<f64> = M22::identity();
        assert_eq!(result.as_vec(), expected.as_vec());
    }

    #[test]
    fn create_identity_ints() {
        let expected = M22::new([[1, 0], [0, 1]]);
        let result: M22<i32> = M22::identity();
        assert_eq!(result.as_vec(), expected.as_vec());
    }

    #[test]
    fn add_m22_floats() {
        let m1 = M22::new([[1.0, 2.0], [3.0, 4.0]]);
        let m2 = M22::new([[5.0, 6.0], [7.0, 8.0]]);
        let expected = M22::new([[6.0, 8.0], [10.0, 12.0]]);
        let result = m1 + m2;
        assert_eq!(result.as_vec(), expected.as_vec());
    }

    #[test]
    fn sub_test() {
        let m1 = m22_new!(1.0, 2.0;
                          3.0, 4.0);
        let m2 = m22_new!(5.0, 6.0;
                          7.0, 8.0);
        let expected = m22_new!( -4.0,  -4.0;
                                 -4.0,  -4.0);
        let result = m1 - m2;
        assert_eq!(result.as_vec(), expected.as_vec());
    }

    #[test]
    fn add_m22_ints() {
        let m1 = M22::new([[1, 2], [3, 4]]);
        let m2 = M22::new([[5, 6], [7, 8]]);
        let expected = M22::new([[6, 8], [10, 12]]);
        let result = m1 + m2;
        assert_eq!(result.as_vec(), expected.as_vec());
    }

    // TODO(elsuizo:2020-06-02): no se como hacer para que ande con Ints y Floats
    // #[test]
    // #[ignore]
    // fn test_determinant() {
    //     let m1 = M22::new([[1, 2], [3, 4]]);
    //     let result = m1.det();
    //     let expected = -2;
    //     assert_eq!(result, expected);
    // }

    #[test]
    fn product_with_vector2_rhs_test() {
        let m1 = M22::new([[1.0, 2.0], [3.0, 4.0]]);
        let v = V2::new([1.0, 2.0]);

        let result = m1 * v;
        let expected = V2::new([5.0, 11.0]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn product_with_matrix2x2_rhs_test() {
        let v = V2::new([1.0, 2.0]);
        let m1 = M22::new([[1.0, 2.0], [3.0, 4.0]]);
        let result = v * m1;
        let expected = V2::new([7.0, 10.0]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn inverse_test() {
        // NOTE(elsuizo:2020-06-02): no se si conviene asi o poner el numero
        // directamente
        use super::test_matrix2x2::EPS;
        let m1 = M22::new([[1.0, 2.0], [3.0, 4.0]]);
        let expected = M22::new([[-2.0, 1.0], [1.5, -0.5]]);
        if let Some(result) = m1.inverse() {
            assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
        }
    }

    #[test]
    fn get_columns_test() {
        let m1 = m22_new!(1.0, 2.0;
                          3.0, 4.0);
        let result = m1.get_cols();

        let expected1 = V2::new([1.0, 3.0]);
        let expected2 = V2::new([2.0, 4.0]);
        let expected = V2::new([expected1, expected2]);
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
        let m1 = m22_new!(1.0, 2.0;
                          3.0, 4.0);
        let result = m1.get_rows();

        let expected1 = V2::new([1.0, 2.0]);
        let expected2 = V2::new([3.0, 4.0]);
        let expected = V2::new([expected1, expected2]);
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
        let expected = m22_new!(1.0, 2.0;
                                3.0, 4.0);

        let cols = expected.get_cols();

        let result = M22::new_from_vecs(cols);

        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn qr_test() {
        let expected = m22_new!(10.0, 2.0;
                                3.0, -4.0);
        if let Some((q, r)) = expected.qr() {
            let result = q * r;
            assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
            assert!(nearly_equal(q.det().abs(), 1.0, EPS));
        }
    }

    #[test]
    fn get_diagonal() {
        let m = m22_new!(10.0, 2.0;
                         3.0, -4.0);
        let result = m.get_diagonal();
        let expected = V2::new([10.0, -4.0]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }
}
