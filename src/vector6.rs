//-------------------------------------------------------------------------
// @file vector6.rs
//
// @date 06/02/20 21:07:05
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

use crate::slices_methods::{norm_inf, norm_l};
use crate::errors::VectorErrors;
use crate::matrix6x6::M66;

//-------------------------------------------------------------------------
//                        code
//-------------------------------------------------------------------------
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct V6<T>([T; 6]);

impl<T> V6<T> {
    /// create a new V6 from a static array
    pub const fn new(input: [T; 6]) -> Self {
        V6(input)
    }

    /// create a new V6 from raw numbers
    pub const fn new_from(num1: T, num2: T, num3: T, num4: T, num5: T, num6: T) -> Self {
        Self::new([num1, num2, num3, num4, num5, num6])
    }
}

impl<T: Num + Copy> V6<T> {
    /// create a V6 with all elements zeros
    pub fn zeros() -> Self {
        <V6<T> as Zero>::zero()
    }

    /// create a V6 with all elements one
    pub fn ones() -> Self {
        let one = T::one();
        Self::new([one, one, one, one, one, one])
    }
}

impl<T: Num + Copy + core::cmp::PartialOrd> V6<T> {
    pub fn norm_inf(&self) -> T {
        norm_inf(&**self)
    }
}

impl<T: Num + Copy + Signed + core::iter::Sum> V6<T> {
    pub fn norm_l(&self) -> T {
        norm_l(&**self)
    }
}

impl<T: Num + Copy + Signed> Neg for V6<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::new_from(-self[0], -self[1], -self[2], -self[3], -self[4], -self[5])
    }
}

impl<T: Float> V6<T> {
    /// calculate the euclidean norm of the V6
    #[inline]
    pub fn norm2(&self) -> T {
        T::sqrt(self[0] * self[0] + self[1] * self[1] + self[2] * self[2]
                + self[3] * self[3] + self[4] * self[4] + self[5] * self[5])
    }

    /// normalize the current V6 and return a new one
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

// V6 * constant
impl<T: Num + Copy> Mul<T> for V6<T> {
    type Output = V6<T>;

    #[inline]
    fn mul(self, rhs: T) -> V6<T> {
        Self::new_from(self[0] * rhs, self[1] * rhs, self[2] * rhs,
                       self[3] * rhs, self[4] * rhs, self[5] * rhs)
    }
}

// V6 / const
impl<T: Num + Copy> Div<T> for V6<T> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        Self::new_from(self[0] / rhs, self[1] / rhs, self[2] / rhs,
                       self[3] / rhs, self[4] / rhs, self[5] / rhs)
    }
}

// f32 * V6<f32>
impl Mul<V6<f32>> for f32 {
    type Output = V6<f32>;

    #[inline]
    fn mul(self, rhs: V6<f32>) -> V6<f32> {
        V6::new_from(self * rhs[0], self * rhs[1], self * rhs[2],
                     self * rhs[3], self * rhs[4], self * rhs[5])
    }
}

// V6 * V6(dot product)
impl<T: Num + Copy> Mul for V6<T> {
    type Output = T;

    #[inline]
    fn mul(self, rhs: Self) -> T {
        self[0] * rhs[0] + self[1] * rhs[1] + self[2] * rhs[2] + self[3] * rhs[3]
        + self[4] * rhs[4] + self[5] * rhs[5]
    }
}

// V6 * M66
impl<T: Num + Copy> Mul<M66<T>> for V6<T> {
    type Output = V6<T>;

    fn mul(self, rhs: M66<T>) -> V6<T> {
        let a_00 = rhs[(0, 0)];
        let a_01 = rhs[(0, 1)];
        let a_02 = rhs[(0, 2)];
        let a_03 = rhs[(0, 3)];
        let a_04 = rhs[(0, 4)];
        let a_05 = rhs[(0, 5)];
        let a_10 = rhs[(1, 0)];
        let a_11 = rhs[(1, 1)];
        let a_12 = rhs[(1, 2)];
        let a_13 = rhs[(1, 3)];
        let a_14 = rhs[(1, 4)];
        let a_15 = rhs[(1, 5)];
        let a_20 = rhs[(2, 0)];
        let a_21 = rhs[(2, 1)];
        let a_22 = rhs[(2, 2)];
        let a_23 = rhs[(2, 3)];
        let a_24 = rhs[(2, 4)];
        let a_25 = rhs[(2, 5)];
        let a_30 = rhs[(3, 0)];
        let a_31 = rhs[(3, 1)];
        let a_32 = rhs[(3, 2)];
        let a_33 = rhs[(3, 3)];
        let a_34 = rhs[(3, 4)];
        let a_35 = rhs[(3, 5)];
        let a_40 = rhs[(4, 0)];
        let a_41 = rhs[(4, 1)];
        let a_42 = rhs[(4, 2)];
        let a_43 = rhs[(4, 3)];
        let a_44 = rhs[(4, 4)];
        let a_45 = rhs[(4, 5)];
        let a_50 = rhs[(5, 0)];
        let a_51 = rhs[(5, 1)];
        let a_52 = rhs[(5, 2)];
        let a_53 = rhs[(5, 3)];
        let a_54 = rhs[(5, 4)];
        let a_55 = rhs[(5, 5)];

        let v0 = self[0];
        let v1 = self[1];
        let v2 = self[2];
        let v3 = self[3];
        let v4 = self[4];
        let v5 = self[5];

        V6::new([
            a_00 * v0 + a_10 * v1 + a_20 * v2 + a_30 * v3 + a_40 * v4 + a_50 * v5,
            a_01 * v0 + a_11 * v1 + a_21 * v2 + a_31 * v3 + a_41 * v4 + a_51 * v5,
            a_02 * v0 + a_12 * v1 + a_22 * v2 + a_32 * v3 + a_42 * v4 + a_52 * v5,
            a_03 * v0 + a_13 * v1 + a_23 * v2 + a_33 * v3 + a_43 * v4 + a_53 * v5,
            a_04 * v0 + a_14 * v1 + a_24 * v2 + a_34 * v3 + a_44 * v4 + a_54 * v5,
            a_05 * v0 + a_15 * v1 + a_25 * v2 + a_35 * v3 + a_45 * v4 + a_55 * v5,
        ])
    }
}

// V6 - V6
impl<T: Num + Copy> Sub for V6<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new_from(self[0] - rhs[0], self[1] - rhs[1], self[2] - rhs[2],
                       self[3] - rhs[3], self[4] - rhs[4], self[5] - rhs[5])
    }
}

// V6 -= V6
impl<T: Num + Copy> SubAssign for V6<T> {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other
    }
}

// V6 + V6
impl<T: Num + Copy> Add for V6<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new_from(self[0] + rhs[0], self[1] + rhs[1], self[2] + rhs[2],
                       self[3] + rhs[3], self[4] + rhs[4], self[5] + rhs[5])
    }
}

// V6 += V6
impl<T: Num + Copy> AddAssign for V6<T> {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        *self = *self + other
    }
}

// impl Zero trait
impl<T: Num + Copy> Zero for V6<T> {
    fn zero() -> V6<T> {
        V6::new([T::zero(); 6])
    }

    fn is_zero(&self) -> bool {
        *self == V6::zero()
    }
}

impl<T> Deref for V6<T> {
    type Target = [T; 6];
    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> DerefMut for V6<T> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T> From<[T; 6]> for V6<T> {
    fn from(data: [T; 6]) -> V6<T> {
        V6(data)
    }
}

//-------------------------------------------------------------------------
//                        Display impl
//-------------------------------------------------------------------------
impl<T: Num + fmt::Display> fmt::Display for V6<T> {
    fn fmt(&self, dest: &mut fmt::Formatter) -> fmt::Result {
        writeln!(dest, "[{0:^3.2} {1:^3.2} {2:^3.2} {3:^3.2} {4:^3.2} {5:^3.2}]",
                 self[0], self[1], self[2], self[3], self[4], self[5])
    }
}

//-------------------------------------------------------------------------
//                        constants
//-------------------------------------------------------------------------
/// constant `V6` zeros
pub const V6_ZEROS: V6<f32> = V6::new_from(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod vector6_tests {

    use crate::matrix6x6::M66;
    use crate::vector6::V6;

    #[test]
    fn vector6_creation_test() {
        let v = V6::new([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(v[0], 1.0);
        assert_eq!(v[1], 2.0);
        assert_eq!(v[2], 3.0);
        assert_eq!(v[3], 4.0);
        assert_eq!(v[4], 5.0);
        assert_eq!(v[5], 6.0);
    }

    #[test]
    fn vector6_add_test() {
        let v = V6::new([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let result = v + v;
        let expected = V6::new([2.0, 4.0, 6.0, 8.0, 10.0, 12.0]);
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
        let v = V6::new([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let result = v - v;
        let expected = V6::new([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
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
        let v = V6::new([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let result = 2.0 * v;
        let expected = V6::new([2.0, 4.0, 6.0, 8.0, 10.0, 12.0]);
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
        let v = V6::new([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let result = v * 2.0;
        let expected = V6::new([2.0, 4.0, 6.0, 8.0, 10.0, 12.0]);
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
        let v = V6::new([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let result = v * v;
        let expected = 91.0;
        assert_eq!(result, expected);
    }
    #[test]
    fn product_matrix6x6_test() {
        let v = V6::new([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);

        let m = M66::new([
            [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            [6.0, 7.0, 8.0, 9.0, 10.0, 11.0],
            [12.0, 13.0, 14.0, 15.0, 16.0, 17.0],
            [18.0, 19.0, 20.0, 21.0, 22.0, 23.0],
            [24.0, 25.0, 26.0, 27.0, 28.0, 29.0],
            [30.0, 31.0, 32.0, 33.0, 34.0, 35.0],
        ]);
        let result = v * m;
        let expected = V6::new([420.0, 441.0, 462.0, 483.0, 504.0, 525.0]);
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
        let result = V6::new([1.0, 1.0, 1.0, 1.0, 1.0, 1.0]).normalize().unwrap();
        let expected = V6::new([0.4082482904638631, 0.4082482904638631, 0.4082482904638631, 0.4082482904638631, 0.4082482904638631, 0.4082482904638631]);
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
        let mut result = V6::new([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let v = V6::new([0.0, 1.0, 2.0, 3.0, 4.0, 5.0]);
        let expected = V6::new([1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
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
        let mut result = V6::new_from(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
        let v = V6::new_from(0.0, 1.0, 2.0, 3.0, 4.0, 5.0);
        let expected = V6::new_from(1.0, 3.0, 5.0, 7.0, 9.0, 11.0);
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
        let v = V6::new_from(1, 10, -10, 10, 100, 3);
        let result = v.norm_inf();
        let expected = 100;
        assert_eq!(result, expected);
    }

    #[test]
    fn norm_l_test() {
        let v = V6::new_from(1, -1, 1, -1, 1, -1);
        let result = v.norm_l();
        let expected = 6;
        assert_eq!(result, expected);
    }
}
