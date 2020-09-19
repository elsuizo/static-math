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
use std::ops::{Mul, Add, Sub, Neg, Div};
use num::{Num, Float, Signed};

use crate::vector3::*;
use crate::matrix3x3::M33;

/// Quaternion type
#[derive(Copy, Debug, Clone)]
pub struct Quaternion<T> {
    /// Scalar part
    pub q0: T,
    /// Imaginary part
    pub q: V3<T>,
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
    /// dot product
    pub fn dot(&self, rhs: Self) -> T {
        self.q0 * rhs.q0 + self.q * rhs.q
    }

    /// get the real part
    pub fn real(&self) -> T {
        self.q0
    }

    /// get the imaginary part
    pub fn imag(&self) -> V3<T> {
        self.q
    }

}

// q + q
impl<T: Num + Copy> Add for Quaternion<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self{q0: self.q0 + rhs.q0, q: self.q + rhs.q}
    }
}

// q - q
impl<T: Num + Copy> Sub for Quaternion<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self{q0: self.q0 - rhs.q0, q: self.q - rhs.q}
    }
}

// q / const
impl<T: Num + Copy> Div<T> for Quaternion<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        let q0 = self.q0 / rhs;
        let q  = self.q / rhs;
        Self{q0, q}
    }
}

// q * q
impl<T: Num + Copy> Mul for Quaternion<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let q0 = self.q0 * rhs.q0  - self.q * rhs.q;
        let q = rhs.q * self.q0 + self.q * rhs.q0 + self.q.cross(rhs.q);
        Self{q0, q}
    }
}

// NOTE(elsuizo:2020-09-10): this implementation comes from this nice simplification
// https://fgiesen.wordpress.com/2019/02/09/rotating-a-single-vector-using-a-quaternion/
// from: Fabian “ryg” Giesen
impl<T: Num + Copy + Signed> Mul<V3<T>> for Quaternion<T> {
    type Output = V3<T>;
    fn mul(self, rhs: V3<T>) -> Self::Output {
        let one = T::one();
        let two = one + one;
        let t = (self.q * two).cross(rhs);
        rhs + t * self.q0 + self.q.cross(t)
    }
}

impl<T: Num + Copy + Signed> Neg for Quaternion<T> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self{q0: -self.q0, q: -self.q}
    }
}

impl<T: Num + Copy + Signed> Quaternion<T> {
    pub fn conj(&self) -> Self {
        Self{q0: self.q0, q: -self.q}
    }
}

// TODO(elsuizo:2020-09-09): maybe here is better a Error
impl<T: Float> Quaternion<T> {
    /// normalize the Quaternion
    pub fn normalize(&self) -> Option<Self> {
        let norm_sqr = self.dot(*self);
        if norm_sqr > T::epsilon() {
            Some(*self / self.dot(*self).sqrt())
        } else {
            None
        }
    }

    /// generate a Quaternion that represents a rotation of a angle `theta`
    /// around the axis(normalized) `v`
    pub fn rotation(theta: T, v: V3<T>) -> Self {
        let two = T::from(2u8).unwrap();
        let n = v.normalize().expect("the input has to be a non zero vector");
        Self::new((theta.to_radians() / two).cos(), n * (theta.to_radians() / two).sin())
    }

    /// generate a Quaternion that represents a rotation of a angle `theta`
    /// around the axis(normalized) `v`, the angle `theta` is encoded in the
    /// norm of the vector `v`
    pub fn rotation_norm_encoded(v: V3<T>) -> Self {
        let one = T::one();
        let two = T::from(2.0).unwrap();
        let theta = v.norm2();
        if theta > T::epsilon() {
            let s = T::sin(theta.to_radians() / two) / theta;
            Self{q0: T::cos(theta.to_radians() / two), q: v * s}
        } else {
            Self{q0: one, q: V3::zeros()}
        }
    }

    pub fn to_dcm(&self) -> M33<T> {
        let (q0, q) = (self.real(), self.imag());
        let q0_s = q0 * q0;
        let (q1, q2, q3) = (q[0], q[1], q[2]);
        let q1_s = q1 * q1;
        let q2_s = q2 * q2;
        let q3_s = q3 * q3;
        let two = T::one() + T::one();

        m33_new!(q0_s + q1_s - q2_s - q3_s, two*q1*q2 - two*q0*q3, two*q1*q3 + two*q0*q2;
                 two*q1*q2 + two*q0*q3, q0_s - q1_s + q2_s - q3_s, two*q2*q3 - two*q0*q1;
                 two*q1*q3 - two*q0*q2, two*q2*q3 + two*q0*q1, q0_s - q1_s - q2_s + q3_s)
    }

    pub fn get_angle(&self) -> T {
        let two = T::from(2.0).unwrap();
        let n = self.q.norm2();

        two * T::atan2(n, self.q0)
    }

    pub fn get_axis(&self) -> Option<V3<T>> {
        let qn = self.normalize()?;
        let s = T::cos(qn.get_angle() / T::from(2.0)?);
        if s.abs() > T::epsilon() {
            Some(qn.q / s)
        } else {
            None
        }
    }
}

//-------------------------------------------------------------------------
//                        testing
//-------------------------------------------------------------------------
#[cfg(test)]
mod test_quaternion {
    use crate::vector3::V3;
    use crate::quaternion::Quaternion;
    use crate::utils::{compare_floats};

    const EPS: f32 = 1e-6;

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
        let q2 = q1.conj();

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
        let result = a.conj();
        assert_eq!(result.q0, 1);
        assert_eq!(result.q[0], -1);
        assert_eq!(result.q[1], -1);
        assert_eq!(result.q[2], -1);


        let a_float = Quaternion::new(1.0, V3::ones());
        let result_float = a_float.conj();
        assert_eq!(result_float.q0, 1.0);
        assert_eq!(result_float.q[0], -1.0);
        assert_eq!(result_float.q[1], -1.0);
        assert_eq!(result_float.q[2], -1.0);
    }

    #[test]
    fn rotate_vec() {
        let q1 = Quaternion::rotation(90.0, V3::new_from(0.0, 0.0, 1.0));
        let x = V3::new_from(1.0, 0.0, 0.0);
        // rotate x around z 90 degrees
        let result = q1 * x;
        let expected = V3::new_from(0.0, -1.0, 0.0);
        assert!(compare_floats(result[0], expected[0], EPS));
        assert!(compare_floats(result[1], expected[1], EPS));
        assert!(compare_floats(result[2], expected[2], EPS));
    }

    #[test]
    fn rotate_vec_composition_360() {
        let q1 = Quaternion::rotation(90.0, V3::new_from(0.0, 0.0, 1.0));
        let x = V3::new_from(1.0, 0.0, 0.0);
        // rotate x around z (90 * 4 = 360) degrees
        let result = q1 * q1 * q1 * q1 * x;
        assert!(compare_floats(result[0], x[0], EPS));
        assert!(compare_floats(result[1], x[1], EPS));
        assert!(compare_floats(result[2], x[2], EPS));
    }

    #[test]
    fn rotate_vec_angle_encode() {
        let q = Quaternion::rotation_norm_encoded(V3::new_from(0.0, 0.0, 90.0));
        let x = V3::x_axis();
        let result = q * x;
        let expected = V3::new_from(0.0, -1.0, 0.0);
        assert!(compare_floats(result[0], expected[0], EPS));
        assert!(compare_floats(result[1], expected[1], EPS));
        assert!(compare_floats(result[2], expected[2], EPS));
    }

    #[test]
    fn convert_dcm_test() {
        let q = Quaternion::rotation_norm_encoded(V3::new_from(0.0, 0.0, 90.0));
        let x = V3::x_axis();
        // rotate the x around z axis 360 degrees
        let expected = q * q * q * q * x;
        // convert the quaternion to a rotation matrix
        let m = q.to_dcm();
        // rotate the x around z axis 360 degrees with the rotation matrix
        let result = m * m * m * m * x;

        assert!(compare_floats(result[0], expected[0], EPS));
        assert!(compare_floats(result[1], expected[1], EPS));
        assert!(compare_floats(result[2], expected[2], EPS));
    }
}
