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
use crate::matrix4x4::M44;

//-------------------------------------------------------------------------
//                        code
//-------------------------------------------------------------------------
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct V4<T>([T; 4]);

impl<T> V4<T> {
    pub fn new(input: [T; 4]) -> V4<T> {
        V4(input)
    }
}

impl<T: Num + Copy> V4<T> {
    pub fn zeros() -> V4<T> {
        <V4<T> as Zero>::zero()
    }

}

impl<T: Float> V4<T> {
    pub fn norm2(&self) -> T {
        let a = self[0];
        let b = self[1];
        let c = self[2];
        let d = self[3];
        T::sqrt(a * a + b * b + c * c + d * d)
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

impl<T: Num + Copy> Mul for V4<T> {
    type Output = T;

    fn mul(self, rhs: Self) -> T {
        let a1 = self[0];
        let a2 = self[1];
        let a3 = self[2];
        let a4 = self[3];

        let b1 = rhs[0];
        let b2 = rhs[1];
        let b3 = rhs[2];
        let b4 = rhs[3];

        a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4
    }
}

/// V4 * constant
impl<T: Num + Copy> Mul<T> for V4<T> {
    type Output = V4<T>;

    fn mul(self, rhs: T) -> V4<T> {
        let a0 = self[0] * rhs;
        let a1 = self[1] * rhs;
        let a2 = self[2] * rhs;
        let a3 = self[3] * rhs;
        V4::new([a0, a1, a2, a3])
    }
}

/// f32 * V4<f32>
impl Mul<V4<f32>> for f32 {
    type Output = V4<f32>;

    fn mul(self, rhs: V4<f32>) -> V4<f32> {
        let a0 = self * rhs[0];
        let a1 = self * rhs[1];
        let a2 = self * rhs[2];
        let a3 = self * rhs[3];
        V4::new([a0, a1, a2, a3])
    }
}

impl<T: Num + Copy> Mul<M44<T>> for V4<T> {
    type Output = V4<T>;

    fn mul(self, rhs: M44<T>) -> V4<T> {
        let a11 = rhs[(0, 0)];
        let a12 = rhs[(0, 1)];
        let a13 = rhs[(0, 2)];
        let a14 = rhs[(0, 3)];
        let a21 = rhs[(1, 0)];
        let a22 = rhs[(1, 1)];
        let a23 = rhs[(1, 2)];
        let a24 = rhs[(1, 3)];
        let a31 = rhs[(2, 0)];
        let a32 = rhs[(2, 1)];
        let a33 = rhs[(2, 2)];
        let a34 = rhs[(2, 3)];
        let a41 = rhs[(3, 0)];
        let a42 = rhs[(3, 1)];
        let a43 = rhs[(3, 2)];
        let a44 = rhs[(3, 3)];

        let v1 = self[0];
        let v2 = self[1];
        let v3 = self[2];
        let v4 = self[3];

        V4::new([
            v1 * a11 + v2 * a21 + v3 * a31 + v4 * a41,
            v1 * a12 + v2 * a22 + v3 * a32 + v4 * a42,
            v1 * a13 + v2 * a23 + v3 * a33 + v4 * a43,
            v1 * a14 + v2 * a24 + v3 * a34 + v4 * a44,
        ])
    }
}

// V4 - V4
impl<T: Num + Copy> Sub for V4<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let v0 = self[0];
        let v1 = self[1];
        let v2 = self[2];
        let v3 = self[3];

        let a0 = rhs[0];
        let a1 = rhs[1];
        let a2 = rhs[2];
        let a3 = rhs[3];

        V4::new([v0 - a0, v1 - a1, v2 - a2, v3 - a3])
    }
}

impl<T: Num + Copy> Add for V4<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let v1 = self[0];
        let v2 = self[1];
        let v3 = self[2];
        let v4 = self[3];

        let a1 = rhs[0];
        let a2 = rhs[1];
        let a3 = rhs[2];
        let a4 = rhs[3];

        V4::new([v1 + a1, v2 + a2, v3 + a3, v4 + a4])
    }
}

impl<T: Num + Copy> Zero for V4<T> {
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
        println!("");
        write!(dest, "[{0:^3.2} {1:^3.2} {2:^3.2} {3:^3.2}]\n", self[0], self[1], self[2], self[3])
    }
}

// TODO(elsuizo:2020-06-02): impl fmt for this type
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
}
