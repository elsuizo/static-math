//-------------------------------------------------------------------------
// @file vector3.rs
//
// @date 06/02/20 18:50:41
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
use core::ops::{Deref, DerefMut};
use num::{Float, Num, Signed, Zero};

use crate::slices_methods::{norm_inf, norm_l};
use core::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

use crate::errors::VectorErrors;
use crate::matrix3x3::M33;
//-------------------------------------------------------------------------
//                        code
//-------------------------------------------------------------------------
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct V3<T>([T; 3]);

impl<T> V3<T> {
    /// create a new V3 from a static array
    pub const fn new(input: [T; 3]) -> Self {
        Self(input)
    }

    /// create a new V3 from raw numbers
    pub const fn new_from(a: T, b: T, c: T) -> Self {
        Self::new([a, b, c])
    }
}

impl<T: Num + Copy> V3<T> {
    /// create a V3 with all elements zero
    #[inline(always)]
    pub fn zeros() -> Self {
        <Self as Zero>::zero()
    }

    /// create a V3 with all elements one
    pub fn ones() -> Self {
        let one = T::one();
        Self::new_from(one, one, one)
    }

    /// calculate the cross product
    #[inline]
    pub fn cross(&self, rhs: Self) -> Self {
        Self::new_from(
            self[1] * rhs[2] - self[2] * rhs[1],
            self[2] * rhs[0] - self[0] * rhs[2],
            self[0] * rhs[1] - self[1] * rhs[0],
        )
    }

    /// create a unitary x axis vector
    pub fn x_axis() -> Self {
        let one = T::one();
        let zero = T::zero();
        Self::new_from(one, zero, zero)
    }

    /// create a unitary y axis vector
    pub fn y_axis() -> Self {
        let one = T::one();
        let zero = T::zero();
        Self::new_from(zero, one, zero)
    }

    /// create a unitary z axis vector
    pub fn z_axis() -> Self {
        let one = T::one();
        let zero = T::zero();
        Self::new_from(zero, zero, one)
    }
}

impl<T: Num + Copy + core::cmp::PartialOrd> V3<T> {
    pub fn norm_inf(&self) -> T {
        norm_inf(&**self)
    }
}

impl<T: Num + Copy + Signed + core::iter::Sum> V3<T> {
    pub fn norm_l(&self) -> T {
        norm_l(&**self)
    }
}

impl<T: Num + Copy + Signed> Neg for V3<T> {
    type Output = Self;

    #[inline(always)]
    fn neg(self) -> Self {
        Self::new_from(-self[0], -self[1], -self[2])
    }
}

impl<T: Float> V3<T> {
    /// calculate the euclidean norm of the V3
    #[inline(always)]
    pub fn norm2(&self) -> T {
        T::sqrt(self[0] * self[0] + self[1] * self[1] + self[2] * self[2])
    }

    /// normalize the current V3
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

//-------------------------------------------------------------------------
//                        maths basic operations
//-------------------------------------------------------------------------

// V3 * const
impl<T: Num + Copy> Mul<T> for V3<T> {
    type Output = Self;

    #[inline(always)]
    fn mul(self, rhs: T) -> Self::Output {
        Self::new_from(self[0] * rhs, self[1] * rhs, self[2] * rhs)
    }
}

// V3 / const
impl<T: Num + Copy> Div<T> for V3<T> {
    type Output = Self;

    #[inline(always)]
    fn div(self, rhs: T) -> Self::Output {
        Self::new_from(self[0] / rhs, self[1] / rhs, self[2] / rhs)
    }
}

// FIXME(elsuizo:2020-06-19): this is a hack
impl Mul<V3<f32>> for f32 {
    type Output = V3<f32>;

    #[inline(always)]
    fn mul(self, rhs: V3<f32>) -> Self::Output {
        V3::new_from(rhs[0] * self, rhs[1] * self, rhs[2] * self)
    }
}

// f64 * V3<64>
impl Mul<V3<f64>> for f64 {
    type Output = V3<f64>;

    #[inline(always)]
    fn mul(self, rhs: V3<f64>) -> Self::Output {
        V3::new_from(rhs[0] * self, rhs[1] * self, rhs[2] * self)
    }
}

// V3 * V3
impl<T: Num + Copy> Mul for V3<T> {
    type Output = T;

    #[inline(always)]
    fn mul(self, rhs: Self) -> T {
        self[0] * rhs[0] + self[1] * rhs[1] + self[2] * rhs[2]
    }
}

// V3 * M33
impl<T: Num + Copy> Mul<M33<T>> for V3<T> {
    type Output = V3<T>;

    #[inline]
    fn mul(self, rhs: M33<T>) -> V3<T> {
        Self::new_from(
            rhs[(0, 0)] * self[0] + rhs[(0, 1)] * self[1] + rhs[(0, 2)] * self[2],
            rhs[(1, 0)] * self[0] + rhs[(1, 1)] * self[1] + rhs[(1, 2)] * self[2],
            rhs[(2, 0)] * self[0] + rhs[(2, 1)] * self[1] + rhs[(2, 2)] * self[2])
    }
}

// V3 - V3
impl<T: Num + Copy> Sub for V3<T> {
    type Output = Self;

    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        Self::new_from(self[0] - rhs[0], self[1] - rhs[1], self[2] - rhs[2])
    }
}

// V3 -= V3
impl<T: Num + Copy> SubAssign for V3<T> {
    #[inline(always)]
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other
    }
}

// V3 + V3
impl<T: Num + Copy> Add for V3<T> {
    type Output = Self;

    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        Self::new_from(self[0] + rhs[0], self[1] + rhs[1], self[2] + rhs[2])
    }
}

// V3 += V3
impl<T: Num + Copy> AddAssign for V3<T> {
    #[inline(always)]
    fn add_assign(&mut self, other: Self) {
        *self = *self + other
    }
}

// impl the Zero trait
impl<T: Num + Copy> Zero for V3<T> {
    #[inline(always)]
    fn zero() -> V3<T> {
        Self::new_from(T::zero(), T::zero(), T::zero())
    }

    fn is_zero(&self) -> bool {
        *self == V3::zero()
    }
}

impl<T> Deref for V3<T> {
    type Target = [T; 3];
    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> DerefMut for V3<T> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T> From<[T; 3]> for V3<T> {
    fn from(data: [T; 3]) -> V3<T> {
        V3(data)
    }
}

//-------------------------------------------------------------------------
//                        Display impl
//-------------------------------------------------------------------------
impl<T: Num + fmt::Display> fmt::Display for V3<T> {
    fn fmt(&self, dest: &mut fmt::Formatter) -> fmt::Result {
        writeln!(
            dest,
            "[{0:^3.2} {1:^3.2} {2:^3.2}]",
            self[0], self[1], self[2]
        )
    }
}

//-------------------------------------------------------------------------
//                        constants
//-------------------------------------------------------------------------
/// constant `x` axis
pub const X_AXIS: V3<f32>   = V3::new_from(1.0, 0.0, 0.0);
/// constant `y` axis
pub const Y_AXIS: V3<f32>   = V3::new_from(0.0, 1.0, 0.0);
/// constant `z` axis
pub const Z_AXIS: V3<f32>   = V3::new_from(0.0, 0.0, 1.0);
/// constant `V3` zeros
pub const V3_ZEROS: V3<f32> = V3::new_from(0.0, 0.0, 0.0);
//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod vector3_test {

    use crate::vector3::V3;

    #[test]
    fn create_vector3_test() {
        let v = V3::new([1.0, 1.0, 1.0]);
        assert_eq!(v[0], 1.0);
        assert_eq!(v[1], 1.0);
        assert_eq!(v[2], 1.0);
    }

    #[test]
    fn zero_vector3_test() {
        let result: V3<f64> = V3::zeros();
        let expected = V3::new([0.0, 0.0, 0.0]);
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
        let v1 = V3::new([1.0, 2.0, 3.0]);
        let v2 = V3::new([4.0, 5.0, 6.0]);
        let result = v1 * v2;
        let expected = 32.0;
        assert_eq!(result, expected);
    }

    #[test]
    fn add_test() {
        let v1 = V3::new([1.0, 2.0, 3.0]);
        let v2 = V3::new([4.0, 5.0, 6.0]);
        let result = v1 + v2;
        let expected = V3::new([5.0, 7.0, 9.0]);
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
        let v1 = V3::new([1.0, 2.0, 3.0]);
        let expected = 3.7416573867739413;
        let result = v1.norm2();
        assert_eq!(result, expected);
    }

    #[test]
    fn mul_const_rhs() {
        let v = V3::new([1.0, 2.0, 3.0]);
        let result: V3<f64> = 2.0 * v;
        let expected = V3::new([2.0, 4.0, 6.0]);
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
        let v = V3::new([1.0, 2.0, 3.0]);
        let result = v * 2.0;
        let expected = V3::new([2.0, 4.0, 6.0]);
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
        let v1 = V3::new([1.0, 1.0, 1.0]);
        let v2 = V3::new([2.0, 3.0, 4.0]);
        let result = v1 - v2;
        let expected = V3::new([-1.0, -2.0, -3.0]);
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
        let result = V3::new([1.0, 1.0, 1.0]).normalize().unwrap();
        let expected = V3::new([0.5773502691896258, 0.5773502691896258, 0.5773502691896258]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn cross_test() {
        let x = V3::new([1.0, 0.0, 0.0]);
        let y = V3::new([0.0, 1.0, 0.0]);

        let result = x.cross(y);
        // z
        let expected = V3::new([0.0, 0.0, 1.0]);
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
        let mut result = V3::new([1.0, 2.0, 3.0]);
        let v = V3::new([4.0, 5.0, 6.0]);
        let expected = V3::new([-3.0, -3.0, -3.0]);
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
        let mut result = V3::new_from(1.0, 2.0, 3.0);
        let v = V3::new_from(4.0, 5.0, 6.0);
        let expected = V3::new_from(5.0, 7.0, 9.0);
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
        let v = V3::new_from(1, -10, 73);
        let result = v.norm_inf();
        let expected = 73;
        assert_eq!(result, expected);
    }

    #[test]
    fn norm_l_test() {
        let v = V3::new_from(1, -1, 1);
        let result = v.norm_l();
        let expected = 3;
        assert_eq!(result, expected);
    }
}
