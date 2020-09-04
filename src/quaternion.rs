//-------------------------------------------------------------------------
// @file quaternions.rs
//
// @date 08/29/20 20:26:13
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
use std::ops::{Mul, Neg};
use num::{Num, Signed};

use crate::vector3::*;

/// Quaternion type
#[derive(Copy, Debug, Clone)]
pub struct Quaternion<T> {
    q0: T,   // scalar part
    q: V3<T> // imaginary values
}

impl<T> Quaternion<T> {

    pub fn new(q0: T, q: V3<T>) -> Self {
        Self {q0, q}
    }

    pub fn new_from(q0: T, q1: T, q2: T, q3: T) -> Self {
        Self {q0, q: V3::new([q1, q2, q3])}
    }
}

impl<T: Num + Copy> Quaternion<T> {
    pub fn dot(&self, rhs: Self) -> T {
        self.q0 * rhs.q0 + self.q * rhs.q
    }
}

impl<T: Num + Copy> Mul for Quaternion<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let q0 = self.q0 * rhs.q0  - self.q * rhs.q;
        let q = rhs.q * self.q0 + self.q * rhs.q0 + self.q.cross(rhs.q);
        Self {q0, q}
    }

}

impl<T: Num + Copy + Signed> Neg for Quaternion<T> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self{q0: self.q0, q: -self.q}
    }
}
//-------------------------------------------------------------------------
//                        testing
//-------------------------------------------------------------------------
#[cfg(test)]
mod test_quaternion {
    use crate::vector3::V3;
    use crate::quaternion::Quaternion;
    // use crate::utils::{compare_vecs, nearly_equal};

    #[test]
    fn quaternion_creation_test() {
        let q = Quaternion::new(0, V3::ones());

        let expected = V3::new([1, 1, 1]);
        assert_eq!(q.q0, 0);
        assert_eq!(
            &q.q[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &q.q[..],
            &expected[..]
        );
    }

    #[test]
    fn quaternion_product_test() {
        let a = Quaternion::new(1, V3::ones());
        let b = Quaternion::new(1, V3::ones());
        let result = a * b;

        assert_eq!(result.q0, -2);
        assert_eq!(result.q[0], 2);
        assert_eq!(result.q[1], 2);
        assert_eq!(result.q[2], 2);

        let q1 = Quaternion::new(1, V3::ones());
        let q2 = -q1;

        let result = q1 * q2;
        let expected = Quaternion::new(q1.dot(q1), V3::zeros());

        assert_eq!(result.q0, expected.q0);
        assert_eq!(result.q[0], expected.q[0]);
        assert_eq!(result.q[1], expected.q[1]);
        assert_eq!(result.q[2], expected.q[2]);
    }

    #[test]
    fn quaternion_conj() {
        let a = Quaternion::new(1, V3::ones());
        let result = -a;
        assert_eq!(result.q0, 1);
        assert_eq!(result.q[0], -1);
        assert_eq!(result.q[1], -1);
        assert_eq!(result.q[2], -1);


        let a_float = Quaternion::new(1.0, V3::ones());
        let result_float = -a_float;
        assert_eq!(result_float.q0, 1.0);
        assert_eq!(result_float.q[0], -1.0);
        assert_eq!(result_float.q[1], -1.0);
        assert_eq!(result_float.q[2], -1.0);
    }
}
