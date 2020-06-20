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
use num::{Float, Zero, Num};
use std::ops::{Deref, DerefMut};

use std::ops::{Add, Mul, Sub};

use crate::errors::VectorErrors;
use crate::matrix5x5::M55;
//-------------------------------------------------------------------------
//                        code
//-------------------------------------------------------------------------
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct V5<T>([T; 5]);

impl<T> V5<T> {
    pub fn new(input: [T; 5]) -> Self {
        Self(input)
    }
}

impl<T: Num + Copy> V5<T> {
    pub fn zeros() -> Self {
        <V5<T> as Zero>::zero()
    }
}

impl<T: Float> V5<T> {
    pub fn norm2(&self) -> T {
        let a = self[0];
        let b = self[1];
        let c = self[2];
        let d = self[3];
        let e = self[4];
        T::sqrt(a * a + b * b + c * c + d * d + e * e)
    }

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
        println!("");
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
}

