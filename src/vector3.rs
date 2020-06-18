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

use std::ops::{Add, Mul};
// use std::ops::{AddAssign, DivAssign, MulAssign, SubAssign};

use crate::errors::VectorErrors;
use crate::matrix3x3::M33;
//-------------------------------------------------------------------------
//                        code
//-------------------------------------------------------------------------
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct V3<T>([T; 3]);

impl<T> V3<T> {
    pub fn new(input: [T; 3]) -> V3<T> {
        V3(input)
    }
}

impl<T: Num + Copy> V3<T> {
    pub fn zeros() -> V3<T> {
        <V3<T> as Zero>::zero()
    }

}

impl<T: Float> V3<T> {
    pub fn norm2(&self) -> T {
        let a = self[0];
        let b = self[1];
        let c = self[2];
        T::sqrt(a * a + b * b + c * c)
    }
}

impl<T: Float> V3<T> {
    pub fn normalize(&mut self) -> Result<Self, VectorErrors> {
        let n = self.norm2();
        if n != T::zero() {
            for i in 0..3 {
                self[i] = self[i] / n;
            }
            Ok(*self)
        } else {
            Err(VectorErrors::Norm2IsZero)
        }
    }
}

impl<T: Num + Copy> Mul<T> for V3<T> {
    type Output = V3<T>;

    fn mul(self, rhs: T) -> V3<T> {
        let a0 = self[0] * rhs;
        let a1 = self[1] * rhs;
        let a2 = self[2] * rhs;
        V3::new([a0, a1, a2])
    }
}

impl<T: Num + Copy> Mul for V3<T> {
    type Output = T;

    fn mul(self, rhs: Self) -> T {
        let a1 = self[0];
        let a2 = self[1];
        let a3 = self[2];

        let b1 = rhs[0];
        let b2 = rhs[1];
        let b3 = rhs[2];

        a1 * b1 + a2 * b2 + a3 * b3
    }
}

impl<T: Num + Copy> Mul<M33<T>> for V3<T> {
    type Output = V3<T>;

    fn mul(self, rhs: M33<T>) -> V3<T> {
        let a11 = rhs[(0, 0)];
        let a12 = rhs[(0, 1)];
        let a13 = rhs[(0, 2)];
        let a21 = rhs[(1, 0)];
        let a22 = rhs[(1, 1)];
        let a23 = rhs[(1, 2)];
        let a31 = rhs[(2, 0)];
        let a32 = rhs[(2, 1)];
        let a33 = rhs[(2, 2)];

        let v1 = self[0];
        let v2 = self[1];
        let v3 = self[2];

        V3::new([
            a11 * v1 + a12 * v2 + a13 * v3,
            a21 * v1 + a22 * v2 + a23 * v3,
            a31 * v1 + a32 * v2 + a33 * v3,
        ])
    }
}

impl<T: Num + Copy> Add for V3<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let v1 = self[0];
        let v2 = self[1];
        let v3 = self[2];

        let a1 = rhs[0];
        let a2 = rhs[1];
        let a3 = rhs[2];

        V3::new([v1 + a1, v2 + a2, v3 + a3])
    }
}

impl<T: Num + Copy> Zero for V3<T> {
    fn zero() -> V3<T> {
        V3::new([T::zero(); 3])
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
        println!("");
        write!(dest, "[{0:^3.2} {1:^3.2} {2:^3.2}]\n", self[0], self[1], self[2])
    }
}

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
// TODO(elsuizo:2020-06-02): faltan tests
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
}
