//-------------------------------------------------------------------------
// @file vector5.rs
//
// @date 06/02/20 20:49:59
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
use std::fmt;
use num::{Float, Zero, Num, Signed};
use std::ops::{Deref, DerefMut};

use std::ops::{Add, Sub, Div, Mul, SubAssign, AddAssign, Neg};

use crate::slices_methods::{norm_inf, norm_l};
use crate::errors::VectorErrors;
use crate::matrix5x5::M55;
//-------------------------------------------------------------------------
//                        code
//-------------------------------------------------------------------------
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct V5<T>([T; 5]);

impl<T> V5<T> {
    /// create a new V5 from a static array
    pub fn new(input: [T; 5]) -> Self {
        Self(input)
    }

    /// create a new V5 from numbers
    pub fn new_from(a: T, b: T, c: T, d: T, e: T) -> Self {
        Self::new([a, b, c, d, e])
    }
}

impl<T: Num + Copy> V5<T> {
    /// create a V5 with all elements zeros
    pub fn zeros() -> Self {
        <V5<T> as Zero>::zero()
    }

    /// create a V5 with all elements one
    pub fn ones() -> Self {
        let one = T::one();
        Self::new([one, one, one, one, one])
    }
}

impl<T: Num + Copy + std::cmp::PartialOrd> V5<T> {
    pub fn norm_inf(&self) -> T {
        norm_inf(&**self)
    }
}

impl<T: Num + Copy + Signed + std::iter::Sum> V5<T> {
    pub fn norm_l(&self) -> T {
        norm_l(&**self)
    }
}

impl<T: Num + Copy + Signed> Neg for V5<T> {
    type Output = Self;

    fn neg(self) -> Self {
        let a = self[0];
        let b = self[1];
        let c = self[2];
        let d = self[3];
        let e = self[4];
        V5::new([-a, -b, -c, -d, -e])
    }
}
impl<T: Float> V5<T> {
    /// calculate the euclidean norm of the V5
    pub fn norm2(&self) -> T {
        let a = self[0];
        let b = self[1];
        let c = self[2];
        let d = self[3];
        let e = self[4];
        T::sqrt(a * a + b * b + c * c + d * d + e * e)
    }

    /// normalize the current V5 and return a new one
    pub fn normalize(&mut self) -> Result<Self, VectorErrors> {
        let n = self.norm2();
        if n != T::zero() {
            // this method return a new fresh and clean vector :)
            let mut result = Self::zeros();
            for i in 0..self.len() {
                result[i] = self[i] / n;
            }
            Ok(result)
        } else {
            Err(VectorErrors::Norm2IsZero)
        }
    }
}

// V5 * V5
impl<T: Num + Copy> Mul for V5<T> {
    type Output = T;

    fn mul(self, rhs: Self) -> T {
        let a0 = self[0];
        let a1 = self[1];
        let a2 = self[2];
        let a3 = self[3];
        let a4 = self[4];

        let b0 = rhs[0];
        let b1 = rhs[1];
        let b2 = rhs[2];
        let b3 = rhs[3];
        let b4 = rhs[4];

        a0 * b0 + a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4
    }
}

// V5 * constant
impl<T: Num + Copy> Mul<T> for V5<T> {
    type Output = V5<T>;

    fn mul(self, rhs: T) -> V5<T> {
        let a0 = self[0] * rhs;
        let a1 = self[1] * rhs;
        let a2 = self[2] * rhs;
        let a3 = self[3] * rhs;
        let a4 = self[4] * rhs;
        V5::new([a0, a1, a2, a3, a4])
    }
}

// V5 / const
impl<T: Num + Copy> Div<T> for V5<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        let a0 = self[0] / rhs;
        let a1 = self[1] / rhs;
        let a2 = self[2] / rhs;
        let a3 = self[3] / rhs;
        let a4 = self[4] / rhs;
        V5::new([a0, a1, a2, a3, a4])
    }
}
// f32 * V5<f32>
impl Mul<V5<f32>> for f32 {
    type Output = V5<f32>;

    fn mul(self, rhs: V5<f32>) -> V5<f32> {
        let a0 = self * rhs[0];
        let a1 = self * rhs[1];
        let a2 = self * rhs[2];
        let a3 = self * rhs[3];
        let a4 = self * rhs[4];
        V5::new([a0, a1, a2, a3, a4])
    }
}

// V5 * M55
impl<T: Num + Copy> Mul<M55<T>> for V5<T> {
    type Output = V5<T>;

    fn mul(self, rhs: M55<T>) -> V5<T> {
        let v0 = self[0];
        let v1 = self[1];
        let v2 = self[2];
        let v3 = self[3];
        let v4 = self[4];

        let a_00 = rhs[(0, 0)];
        let a_01 = rhs[(0, 1)];
        let a_02 = rhs[(0, 2)];
        let a_03 = rhs[(0, 3)];
        let a_04 = rhs[(0, 4)];
        let a_10 = rhs[(1, 0)];
        let a_11 = rhs[(1, 1)];
        let a_12 = rhs[(1, 2)];
        let a_13 = rhs[(1, 3)];
        let a_14 = rhs[(1, 4)];
        let a_20 = rhs[(2, 0)];
        let a_21 = rhs[(2, 1)];
        let a_22 = rhs[(2, 2)];
        let a_23 = rhs[(2, 3)];
        let a_24 = rhs[(2, 4)];
        let a_30 = rhs[(3, 0)];
        let a_31 = rhs[(3, 1)];
        let a_32 = rhs[(3, 2)];
        let a_33 = rhs[(3, 3)];
        let a_34 = rhs[(3, 4)];
        let a_40 = rhs[(4, 0)];
        let a_41 = rhs[(4, 1)];
        let a_42 = rhs[(4, 2)];
        let a_43 = rhs[(4, 3)];
        let a_44 = rhs[(4, 4)];

        V5::new([
            a_00 * v0 + a_10 * v1 + a_20 * v2 + a_30 * v3 + a_40 * v4,
            a_01 * v0 + a_11 * v1 + a_21 * v2 + a_31 * v3 + a_41 * v4,
            a_02 * v0 + a_12 * v1 + a_22 * v2 + a_32 * v3 + a_42 * v4,
            a_03 * v0 + a_13 * v1 + a_23 * v2 + a_33 * v3 + a_43 * v4,
            a_04 * v0 + a_14 * v1 + a_24 * v2 + a_34 * v3 + a_44 * v4,
        ])
    }
}

// V5 - V5
impl<T: Num + Copy> Sub for V5<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let v0 = self[0];
        let v1 = self[1];
        let v2 = self[2];
        let v3 = self[3];
        let v4 = self[4];

        let a0 = rhs[0];
        let a1 = rhs[1];
        let a2 = rhs[2];
        let a3 = rhs[3];
        let a4 = rhs[4];

        V5::new([v0 - a0, v1 - a1, v2 - a2, v3 - a3, v4 - a4])
    }
}

// V5 -= V5
impl<T: Num + Copy> SubAssign for V5<T> {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other
    }
}

// V5 + V5
impl<T: Num + Copy> Add for V5<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let v0 = self[0];
        let v1 = self[1];
        let v2 = self[2];
        let v3 = self[3];
        let v4 = self[4];

        let a0 = rhs[0];
        let a1 = rhs[1];
        let a2 = rhs[2];
        let a3 = rhs[3];
        let a4 = rhs[4];

        V5::new([v0 + a0, v1 + a1, v2 + a2, v3 + a3, v4 + a4])
    }
}

// V5 += V5
impl<T: Num + Copy> AddAssign for V5<T> {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other
    }
}

impl<T: Num + Copy> Zero for V5<T> {
    fn zero() -> V5<T> {
        V5::new([T::zero(); 5])
    }

    fn is_zero(&self) -> bool {
        *self == V5::zero()
    }
}

impl<T> Deref for V5<T> {
    type Target = [T; 5];
    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> DerefMut for V5<T> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T> From<[T; 5]> for V5<T> {
    fn from(data: [T; 5]) -> V5<T> {
        V5(data)
    }
}

//-------------------------------------------------------------------------
//                        Display impl
//-------------------------------------------------------------------------
impl<T: Num + fmt::Display> fmt::Display for V5<T> {
    fn fmt(&self, dest: &mut fmt::Formatter) -> fmt::Result {
        write!(dest, "[{0:^3.2} {1:^3.2} {2:^3.2} {3:^3.2} {4:^3.2}]\n", self[0], self[1], self[2], self[3], self[4])
    }
}

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod vector5_test {

    use crate::matrix5x5::M55;
    use crate::vector5::V5;

    #[test]
    fn vector5_creation_test() {
        let v = V5::new([1.0, 2.0, 3.0, 4.0, 5.0]);
        assert_eq!(v[0], 1.0);
        assert_eq!(v[1], 2.0);
        assert_eq!(v[2], 3.0);
        assert_eq!(v[3], 4.0);
        assert_eq!(v[4], 5.0);
    }

    #[test]
    fn vector5_zeros_test() {
        let result: V5<f32> = V5::zeros();
        let expected = V5::new([0.0, 0.0, 0.0, 0.0, 0.0]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn vector5_add_test() {
        let v1 = V5::new([1.0, 2.0, 3.0, 4.0, 5.0]);
        let v2 = V5::new([6.0, 7.0, 8.0, 9.0, 10.0]);
        let result = v1 + v2;
        let expected = V5::new([7.0, 9.0, 11.0, 13.0, 15.0]);
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
        let v1 = V5::new([1.0, 2.0, 3.0, 4.0, 5.0]);
        let v2 = V5::new([6.0, 7.0, 8.0, 9.0, 10.0]);
        let result = v1 - v2;
        let expected = V5::new([-5.0, -5.0, -5.0, -5.0, -5.0]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn vector5_mul_test() {
        let v1 = V5::new([1.0, 2.0, 3.0, 4.0, 5.0]);
        let v2 = V5::new([6.0, 7.0, 8.0, 9.0, 10.0]);
        let result = v1 * v2;
        let expected = 130.0;
        assert_eq!(result, expected);
    }

    #[test]
    fn vector5_norm_test() {
        let v1 = V5::new([1.0, 2.0, 3.0, 4.0, 5.0]);
        let result = v1.norm2();
        let expected = 7.416198487095663;
        assert_eq!(result, expected);
    }

    #[test]
    fn mul_const_rhs() {
        let v = V5::new([1.0, 2.0, 3.0, 4.0, 5.0]);
        let result = 2.0 * v;
        let expected = V5::new([2.0, 4.0, 6.0, 8.0, 10.0]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn mul_const() {
        let v = V5::new([1.0, 2.0, 3.0, 4.0, 5.0]);
        let result = v * 2.0;
        let expected = V5::new([2.0, 4.0, 6.0, 8.0, 10.0]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn vector5_mul_matrix5x5_test() {
        let v1 = V5::new([1.0, 2.0, 3.0, 4.0, 5.0]);
        let m = M55::new([
            [10.0, 1.0, 7.0, 1.0, 5.0],
            [2.0, 4.0, 8.0, 3.0, 2.0],
            [5.0, 1.0, 2.0, 9.0, 10.0],
            [6.0, 9.0, 9.0, 7.0, 3.0],
            [1.0, 8.0, 8.0, 10.0, 5.0],
        ]);
        let result = v1 * m;
        let expected = V5::new([58.0, 88.0, 105.0, 112.0, 76.0]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn normalize_test() {
        let result = V5::new([1.0, 1.0, 1.0, 1.0, 1.0]).normalize().unwrap();
        let expected = V5::new([0.4472135954999579, 0.4472135954999579, 0.4472135954999579, 0.4472135954999579, 0.4472135954999579]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn sub_assigment_test() {
        let mut result = V5::new([1.0, 2.0, 3.0, 4.0, 5.0]);
        let v = V5::new([0.0, 1.0, 2.0, 3.0, 4.0]);
        let expected = V5::new([1.0, 1.0, 1.0, 1.0, 1.0]);
        result -= v;
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn add_assigment_test() {
        let mut result = V5::new_from(1.0, 2.0, 3.0, 4.0, 5.0);
        let v = V5::new_from(0.0, 1.0, 2.0, 3.0, 4.0);
        let expected = V5::new_from(1.0, 3.0, 5.0, 7.0, 9.0);
        result += v;
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn norm_inf_test() {
        let v = V5::new_from(1, 2, 30, 91, 10);
        let result = v.norm_inf();
        let expected = 91;
        assert_eq!(result, expected);
    }

    #[test]
    fn norm_l_test() {
        let v = V5::new_from(-1, 1, -1, 1, -1);
        let result = v.norm_l();
        let expected = 5;
        assert_eq!(result, expected);
    }
}

