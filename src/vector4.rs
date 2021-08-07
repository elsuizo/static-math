//-------------------------------------------------------------------------
// @file vector4.rs
//
// @date 06/02/20 20:08:45
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
use core::ops::BitAnd;
use crate::slices_methods::{norm_inf, norm_l};
use crate::errors::VectorErrors;
use crate::matrix4x4::M44;

//-------------------------------------------------------------------------
//                        code
//-------------------------------------------------------------------------
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct V4<T>([T; 4]);

impl<T> V4<T> {
    /// create a new V4 from a static array
    pub const fn new(input: [T; 4]) -> Self {
        Self(input)
    }

    /// create a new V4 from raw numbers
    pub const fn new_from(a: T, b: T, c: T, d: T) -> Self {
        Self::new([a, b, c, d])
    }
}

impl<T: Num + Copy> V4<T> {
    /// create a V4 with all elements zeros
    pub fn zeros() -> V4<T> {
        <V4<T> as Zero>::zero()
    }

    /// create a V4 with all elements one
    pub fn ones() -> Self {
        let one = T::one();
        Self::new([one, one, one, one])
    }

}

impl BitAnd for V4<bool> {
    type Output = Self;
    fn bitand(self, rhs: Self) -> Self::Output {
        Self::new_from(self[0] & rhs[0], self[1] & rhs[1], self[2] & rhs[2], self[3] & rhs[3])
    }
}

// norm inf
impl<T: Num + Copy + core::cmp::PartialOrd> V4<T> {
    pub fn norm_inf(&self) -> T {
        norm_inf(&**self)
    }
}

// norm l
impl<T: Num + Copy + Signed + core::iter::Sum> V4<T> {
    pub fn norm_l(&self) -> T {
        norm_l(&**self)
    }
}

// -V4
impl<T: Num + Copy + Signed> Neg for V4<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        V4::new_from(-self[0], -self[1], -self[2], -self[3])
    }
}

impl<T: Float> V4<T> {
    /// calculate the euclidean norm of the V4
    #[inline]
    pub fn norm2(&self) -> T {
        T::sqrt(self[0] * self[0] + self[1] * self[1] + self[2] * self[2] + self[3] * self[3])
    }

    /// normalize the current vector and return a new one
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

// V4 * V4(dot product)
impl<T: Num + Copy> Mul for V4<T> {
    type Output = T;

    #[inline]
    fn mul(self, rhs: Self) -> T {
        self[0] * rhs[0] + self[1] * rhs[1] + self[2] * rhs[2] + self[3] * rhs[3]
    }
}

// V4 * const
impl<T: Num + Copy> Mul<T> for V4<T> {
    type Output = V4<T>;

    #[inline(always)]
    fn mul(self, rhs: T) -> V4<T> {
        Self::new_from(self[0] * rhs, self[1] * rhs, self[2] * rhs, self[3] * rhs)
    }
}

// NOTE(elsuizo:2021-08-03): this should panics if rhs == zero
// V4 / const
impl<T: Num + Copy> Div<T> for V4<T> {
    type Output = Self;

    #[inline(always)]
    fn div(self, rhs: T) -> Self::Output {
        Self::new_from(self[0] / rhs, self[1] / rhs, self[2] / rhs, self[3] / rhs)
    }
}

// f32 * V4<f32>
impl Mul<V4<f32>> for f32 {
    type Output = V4<f32>;

    #[inline]
    fn mul(self, rhs: V4<f32>) -> V4<f32> {
        V4::new_from(self * rhs[0], self * rhs[1], self * rhs[2], self * rhs[3])
    }
}

// V4 * M44
impl<T: Num + Copy> Mul<M44<T>> for V4<T> {
    type Output = V4<T>;

    #[inline]
    fn mul(self, rhs: M44<T>) -> V4<T> {
        Self::new_from(
            self[0] * rhs[(0, 0)] + self[1] * rhs[(1, 0)] + self[2] * rhs[(2, 0)] + self[3] * rhs[(3, 0)],
            self[0] * rhs[(0, 1)] + self[1] * rhs[(1, 1)] + self[2] * rhs[(2, 1)] + self[3] * rhs[(3, 1)],
            self[0] * rhs[(0, 2)] + self[1] * rhs[(1, 2)] + self[2] * rhs[(2, 2)] + self[3] * rhs[(3, 2)],
            self[0] * rhs[(0, 3)] + self[1] * rhs[(1, 3)] + self[2] * rhs[(2, 3)] + self[3] * rhs[(3, 3)],
        )
    }
}

// V4 - V4
impl<T: Num + Copy> Sub for V4<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        V4::new_from(self[0] - rhs[0], self[1] - rhs[1], self[2] - rhs[2], self[3] - rhs[3])
    }
}

// V4 -= V4
impl<T: Num + Copy> SubAssign for V4<T> {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other
    }
}

// V4 + V4
impl<T: Num + Copy> Add for V4<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        V4::new_from(self[0] + rhs[0], self[1] + rhs[1], self[2] + rhs[2], self[3] + rhs[3])
    }
}

// V4 += V4
impl<T: Num + Copy> AddAssign for V4<T> {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        *self = *self + other
    }
}

impl<T: Num + Copy> Zero for V4<T> {
    #[inline]
    fn zero() -> V4<T> {
        V4::new([T::zero(); 4])
    }

    fn is_zero(&self) -> bool {
        *self == V4::zero()
    }
}

impl<T> Deref for V4<T> {
    type Target = [T; 4];
    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> DerefMut for V4<T> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T> From<[T; 4]> for V4<T> {
    fn from(data: [T; 4]) -> V4<T> {
        V4(data)
    }
}


//-------------------------------------------------------------------------
//                        Display impl
//-------------------------------------------------------------------------
impl<T: Num + fmt::Display> fmt::Display for V4<T> {
    fn fmt(&self, dest: &mut fmt::Formatter) -> fmt::Result {
        writeln!(dest, "[{0:^3.2} {1:^3.2} {2:^3.2} {3:^3.2}]", self[0], self[1], self[2], self[3])
    }
}

//-------------------------------------------------------------------------
//                        constants
//-------------------------------------------------------------------------
/// constant `V4` zeros
pub const V4_ZEROS: V4<f32> = V4::new_from(0.0, 0.0, 0.0, 0.0);

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod vector4_test {
    use crate::matrix4x4::M44;
    use crate::vector4::V4;

    #[test]
    fn vector4_creation_test() {
        let v = V4::new([1, 1, 1, 1]);
        assert_eq!(v[0], 1);
        assert_eq!(v[1], 1);
        assert_eq!(v[2], 1);
        assert_eq!(v[3], 1);
    }

    #[test]
    fn vector4_zeros_test() {
        let result: V4<f32> = V4::zeros();
        let expected = V4::new([0.0, 0.0, 0.0, 0.0]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn vector4_add_test() {
        let v1 = V4::new([1.0, 2.0, 3.0, 4.0]);
        let v2 = V4::new([5.0, 6.0, 7.0, 8.0]);
        let result = v1 + v2;
        let expected = V4::new([6.0, 8.0, 10.0, 12.0]);
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
        let v1 = V4::new([1.0, 2.0, 3.0, 4.0]);
        let v2 = V4::new([5.0, 6.0, 7.0, 8.0]);
        let result = v1 - v2;
        let expected = V4::new([-4.0, -4.0, -4.0, -4.0]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn vector4_product_test() {
        let v1 = V4::new([1.0, 2.0, 3.0, 4.0]);
        let v2 = V4::new([5.0, 6.0, 7.0, 8.0]);
        let result = v1 * v2;
        let expected = 70.0;
        assert_eq!(result, expected);
    }

    #[test]
    fn vector4_norm_test() {
        let v1 = V4::new([1.0, 2.0, 3.0, 4.0]);
        let result = v1.norm2();
        let expected = 5.477225575051661;
        assert_eq!(result, expected);
    }

    #[test]
    fn mul_const_rhs() {
        let v = V4::new([1.0, 2.0, 3.0, 4.0]);
        let result = 2.0 * v;
        let expected = V4::new([2.0, 4.0, 6.0, 8.0]);
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
        let v = V4::new([1.0, 2.0, 3.0, 4.0]);
        let result = v * 2.0;
        let expected = V4::new([2.0, 4.0, 6.0, 8.0]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn vector4_mul_matrix4x4_test() {
        let m = M44::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0, 16.0],
        ]);

        let v1 = V4::new([1.0, 2.0, 3.0, 4.0]);
        let result = v1 * m;
        let expected = V4::new([90.0, 100.0, 110.0, 120.0]);
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
        let result = V4::new([1.0, 1.0, 1.0, 1.0]).normalize().unwrap();
        let expected = V4::new([0.5, 0.5, 0.5, 0.5]);
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
        let mut result = V4::new([1.0, 2.0, 3.0, 4.0]);
        let v = V4::new([5.0, 6.0, 7.0, 8.0]);
        let expected = V4::new([-4.0, -4.0, -4.0, -4.0]);
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
        let mut result = V4::new_from(1.0, 2.0, 3.0, 4.0);
        let v = V4::new_from(5.0, 6.0, 7.0, 8.0);
        let expected = V4::new_from(6.0, 8.0, 10.0, 12.0);
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
        let v = V4::new_from(10, 100, -9, 0);
        let result = v.norm_inf();
        let expected = 100;
        assert_eq!(result, expected);
    }

    #[test]
    fn norm_l_test() {
        let v = V4::new_from(1, -1, 1, -1);
        let result = v.norm_l();
        let expected = 4;
        assert_eq!(result, expected);
    }

    #[test]
    fn product_lhs_f32_test() {
        let v = V4::new_from(1.0f32, 2.0, 3.0, 4.0);
        let result = 2.0f32 * v;
        let expected = V4::new_from(2.0, 4.0, 6.0, 8.0);
        assert_eq!(result, expected);
    }
}
