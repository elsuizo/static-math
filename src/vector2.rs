//-------------------------------------------------------------------------
// @file vector2.rs
//
// @date 06/01/20 22:37:02
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
use num::{Float, Zero, Num, Signed};
use core::ops::{Deref, DerefMut};
use core::ops::{Add, Sub, Div, Mul, SubAssign, AddAssign, Neg};

use crate::errors::VectorErrors;
use crate::slices_methods::{norm_inf, norm_l};
use crate::matrix2x2::*;

//-------------------------------------------------------------------------
//                        code
//-------------------------------------------------------------------------
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct V2<T>([T; 2]);

// NOTE(elsuizo:2020-06-27): whit this definition you could create a new vector
// of any types
impl<T> V2<T> {
    /// create a new V2 from a static array
    pub fn new(input: [T; 2]) -> Self {
        Self(input)
    }

    /// create a new V2 from numbers
    pub fn new_from(a: T, b: T) -> Self {
        Self::new([a, b])
    }
}

impl<T: Num + Copy> V2<T> {
    /// create a V2 with all elements zero
    pub fn zeros() -> Self {
        <Self as Zero>::zero()
    }

    /// create a V2 with all elements one
    pub fn ones() -> Self {
        let one = T::one();
        Self::new([one, one])
    }

}

impl<T: Num + Copy + core::cmp::PartialOrd> V2<T> {
    pub fn norm_inf(&self) -> T {
        norm_inf(&**self)
    }
}

impl<T: Num + Copy + Signed + core::iter::Sum> V2<T> {
    pub fn norm_l(&self) -> T {
        norm_l(&**self)
    }
}

impl<T: Num + Copy + Signed> Neg for V2<T> {
    type Output = Self;

    fn neg(self) -> Self {
        let a = self[0];
        let b = self[1];
        V2::new([-a, -b])
    }
}

// NOTE(elsuizo:2020-06-27): this impl must be whith the Float trait because
// the use of the sqrt method
impl<T: Float> V2<T> {

    /// calculate the euclidean norm of the V2
    pub fn norm2(&self) -> T {
        let a = self[0];
        let b = self[1];

        T::sqrt(a * a + b * b)
    }

    /// normalize the current V2 and return a new one
    pub fn normalize(&self) -> Result<Self, VectorErrors> {
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

impl<T: Num + Copy> Mul for V2<T> {
    type Output = T;

    fn mul(self, rhs: Self) -> T {
        let a1 = self[0];
        let a2 = self[1];

        let b1 = rhs[0];
        let b2 = rhs[1];

        a1 * b1 + a2 * b2
    }
}

// TODO(elsuizo:2020-05-01): missing constant * V2
// V2 * constant
impl<T: Num + Copy> Mul<T> for V2<T> {
    type Output = V2<T>;

    fn mul(self, rhs: T) -> V2<T> {
        let a0 = self[0] * rhs;
        let a1 = self[1] * rhs;
        V2::new([a0, a1])
    }
}

// V2 / const
impl<T: Num + Copy> Div<T> for V2<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        let a0 = self[0] / rhs;
        let a1 = self[1] / rhs;
        V2::new([a0, a1])
    }
}

// FIXME(elsuizo:2020-06-19): this code dont compile i dont known why
//impl<T: Num + Copy> Mul<V2<T>> for T {
//    type Output = V2<T>;
//
//    fn mul(self, rhs: V2<T>) -> V2<T> {
//        let a0 = self * rhs[0];
//        let a1 = self * rhs[1];
//        V2::new([a0, a1])
//    }
//}

impl Mul<V2<f32>> for f32 {
    type Output = V2<f32>;

    fn mul(self, rhs: V2<f32>) -> V2<f32> {
        let a0 = self * rhs[0];
        let a1 = self * rhs[1];
        V2::new([a0, a1])
    }
}

// V2 * M22
impl<T: Num + Copy> Mul<M22<T>> for V2<T> {
    type Output = V2<T>;

    fn mul(self, rhs: M22<T>) -> V2<T> {
        let a1 = rhs[(0, 0)];
        let b1 = rhs[(0, 1)];
        let c1 = rhs[(1, 0)];
        let d1 = rhs[(1, 1)];

        let v1 = self[0];
        let v2 = self[1];
        V2::new([v1 * a1 + v2 * c1, v1 * b1 + v2 * d1])
    }
}

// V2 - V2
impl<T: Num + Copy> Sub for V2<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let v0 = self[0];
        let v1 = self[1];

        let a0 = rhs[0];
        let a1 = rhs[1];

        V2::new([v0 - a0, v1 - a1])
    }
}

// V2 -= V2
impl<T: Num + Copy> SubAssign for V2<T> {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other
    }
}

// V2 + V2
impl<T: Num + Copy> Add for V2<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let a1 = self[0];
        let a2 = self[1];

        let b1 = rhs[0];
        let b2 = rhs[1];

        V2::new([a1 + b1, a2 + b2])
    }
}

// V2 += V2
impl<T: Num + Copy> AddAssign for V2<T> {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other
    }
}

// Zero trait
impl<T: Num + Copy> Zero for V2<T> {
    fn zero() -> V2<T> {
        V2::new([T::zero(); 2])
    }

    fn is_zero(&self) -> bool {
        *self == V2::zero()
    }
}

impl<T> Deref for V2<T> {
    type Target = [T; 2];
    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> DerefMut for V2<T> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

//-------------------------------------------------------------------------
//                        Display impl
//-------------------------------------------------------------------------
impl<T: Num + fmt::Display> fmt::Display for V2<T> {
    fn fmt(&self, dest: &mut fmt::Formatter) -> fmt::Result {
        writeln!(dest, "[{0:^3.2} {1:^3.2}]", self[0], self[1])
    }
}

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod vector2_test {

    use crate::vector2::V2;

    #[test]
    fn create_vector2_test() {
        let v = V2::new([1.0, 2.0]);
        assert_eq!(v[0], 1.0);
        assert_eq!(v[1], 2.0);
    }

    #[test]
    fn zero_vector2_test() {
        let result: V2<f64> = V2::zeros();
        let expected = V2::new([0.0, 0.0]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn product_test() {
        let v1 = V2::new([1.0, 2.0]);
        let v2 = V2::new([3.0, 4.0]);
        let result = v1 * v2;
        let expected = 11.0;
        assert_eq!(result, expected);
    }

    #[test]
    fn add_test() {
        let v1 = V2::new([1.0, 2.0]);
        let v2 = V2::new([3.0, 4.0]);
        let result = v1 + v2;
        let expected = V2::new([4.0, 6.0]);
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
        let v1 = V2::new([1.0, 2.0]);
        let v2 = V2::new([2.0, 3.0]);
        let expected = V2::new([-1.0, -1.0]);
        let result = v1 - v2;

        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn mul_const_rhs() {
        let v1 = V2::new_from(1.0, 2.0);
        let result = 2.0 * v1;
        let expected = V2::new_from(2.0, 4.0);
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
        let v1 = V2::new_from(1.0, 2.0);
        let result = v1 * 10.0;
        let expected = V2::new_from(10.0, 20.0);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn norm2_test() {
        let v1 = V2::new_from(1.0, 2.0);
        let expected = 2.23606797749979;
        let result = v1.norm2();
        assert_eq!(result, expected);
    }

    #[test]
    fn normalize_test() {
        let result = V2::new([1.0, 1.0]).normalize().unwrap();
        let expected = V2::new([0.7071067811865475, 0.7071067811865475]);
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
        let mut result = V2::new_from(1.0, 3.0);
        let v2 = V2::new_from(3.0, 3.0);
        let expected = V2::new_from(-2.0, 0.0);
        result -= v2;
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
        let mut result = V2::new_from(1.0, 3.0);
        let v2 = V2::new_from(3.0, 3.0);
        let expected = V2::new_from(4.0, 6.0);
        result += v2;
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
        let v = V2::new_from(1, 10);
        let result = v.norm_inf();
        let expected = 10;
        assert_eq!(result, expected);
    }

    #[test]
    fn norm_l_test() {
        let v = V2::new_from(-1, 1);
        let result = v.norm_l();
        let expected = 2;
        assert_eq!(result, expected);
    }
}
