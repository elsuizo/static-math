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
use num::{Num, Float, Signed, Zero, One};
use num::traits::FloatConst;
use crate::vector3::*;
use crate::matrix3x3::M33;

/// Quaternion type
#[derive(Copy, Debug, Clone, PartialEq)]
pub struct Quaternion<T> {
    /// Scalar part
    pub q0: T,
    /// Imaginary part
    pub q: V3<T>,
    /// flag to signaling if the Quaternion is normalized
    normalized: bool
}

impl<T> Quaternion<T> {

    /// construct a new Quaternion from a number(real part) and a vector(imag part)
    pub const fn new(q0: T, q: V3<T>) -> Self {
        Self {q0, q, normalized: false}
    }

    /// construct a new Quaternion from four numbers
    pub fn new_from(q0: T, q1: T, q2: T, q3: T) -> Self {
        Self {q0, q: V3::new([q1, q2, q3]), normalized: false}
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

    /// construct a unit Quaternion
    pub fn one() -> Quaternion<T> {
        <Quaternion<T> as One>::one()
    }

    /// construct a zero Quaternion
    pub fn zero() -> Self {
        <Quaternion<T> as Zero>::zero()
    }

    /// construct a pure "real" Quaternion
    pub fn new_real(q0: T) -> Self {
        Self {q0, q: V3::zeros(), normalized: true}
    }

    /// construct a pure "imaginary" Quaternion
    pub fn new_imag(q: V3<T>) -> Self {
        Self {q0: T::zero(), q, normalized: false}
    }

    /// calculate the abs2 of the Quaternion
    pub fn abs2(&self) -> T {
        self.q0 * self.q0 + self.q[0] * self.q[0] + self.q[1] * self.q[1] + self.q[2] * self.q[2]
    }

}

impl<T: Num + Copy> Zero for Quaternion<T> {
    fn zero() -> Self {
        Self::new(T::zero(), V3::zeros())
    }

    fn is_zero(&self) -> bool {
        *self == Quaternion::zero()
    }
}

impl<T: Num + Copy> One for Quaternion<T> {
    /// Create an identity Quaternion
    fn one() -> Self {
        let one = T::one();
        Self::new(one, V3::zeros())
    }
}

// q + q
impl<T: Num + Copy> Add for Quaternion<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.q0 + rhs.q0, self.q + rhs.q)
    }
}

// q - q
impl<T: Num + Copy> Sub for Quaternion<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.q0 - rhs.q0, self.q - rhs.q)
    }
}

// q / const
impl<T: Num + Copy> Div<T> for Quaternion<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        let q0 = self.q0 / rhs;
        let q  = self.q / rhs;
        Self::new(q0, q)
    }
}

// q / q
impl<T: Float + Signed> Div for Quaternion<T> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inv().expect("the input has to be a non zero vector")
    }
}

// q * q
impl<T: Num + Copy> Mul for Quaternion<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let q0 = self.q0 * rhs.q0  - self.q * rhs.q;
        let q = rhs.q * self.q0 + self.q * rhs.q0 + self.q.cross(rhs.q);
        Self::new(q0, q)
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

impl<T: Num + Copy> Mul<T> for Quaternion<T> {
    type Output = Quaternion<T>;
    fn mul(self, rhs: T) -> Self::Output {
        Self {q0: self.q0 * rhs, q: self.q * rhs, normalized: false}
    }
}

impl<T: Num + Copy + Signed> Neg for Quaternion<T> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.q0, -self.q)
    }
}

impl<T: Num + Copy + Signed> Quaternion<T> {
    pub fn conj(&self) -> Self {
        Self::new(self.q0, -self.q)
    }
}

impl<T: Float + Signed> Quaternion<T> {
    pub fn inv(&self) -> Option<Self> {
        if !self.normalized {
            let norm_sqr = self.dot(*self);
            if norm_sqr > T::epsilon() {
                Some(self.conj() / self.abs2())
            } else {
                None
            }
        } else {
            Some(self.conj())
        }
    }
}

// TODO(elsuizo:2021-04-14): a warning from here but i dont known why???
impl<T: Float + FloatConst + Signed> Quaternion<T> {
    /// normalize the Quaternion only if necesary
    pub fn normalize(&self) -> Option<Self> {
        if self.normalized {
            Some(*self)
        } else {
            let norm_sqr = self.dot(*self);
            if norm_sqr > T::epsilon() {
                let mut result = Self::zero();
                result = *self / self.dot(*self).sqrt();
                result.normalized = true;
                Some(result)
            } else {
                None
            }
        }
    }

    /// get the norm of the "imaginary" part
    pub fn abs_imag(&self) -> T {
        self.imag().norm2()
    }

    pub fn norm2(&self) -> T {
        self.dot(*self).sqrt()
    }

    /// generate a Quaternion that represents a rotation of a angle `theta`
    /// around the axis(normalized) `v`
    pub fn rotation(theta: T, v: V3<T>) -> Self {
        let two = T::from(2u8).unwrap();
        let n   = v.normalize().expect("the input has to be a non zero vector");
        let q0  = (theta / two).cos();
        let q   = n * (theta / two).sin();
        Self{q0, q, normalized: true}
    }

    /// generate a Quaternion that represents a rotation of a angle `theta`
    /// around the axis(normalized) `v`, the angle `theta` is encoded in the
    /// norm of the vector `v`
    pub fn rotation_norm_encoded(v: V3<T>) -> Self {
        let one = T::one();
        let two = T::from(2.0).unwrap();
        let theta = v.norm2();
        if theta > T::epsilon() {
            let s = T::sin(theta / two) / theta;
            Self::new(T::cos(theta / two), v * s)
        } else {
            Self::new(one, V3::zeros())
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

    // NOTE(elsuizo:2021-04-14): this implementation avoid the use of sqrt(which is expensive
    // computationally)
    // TODO(elsuizo:2021-04-14): maybe here we could acept a tuple (T, T, T)
    /// create a quaternion that represents the rotation from a Euler angles
    /// with the roll-pitch-yay convention
    pub fn from_euler_angles(yay: T, pitch: T, roll: T) -> Self {
        let two = T::one() + T::one();
        let c1 = T::cos(yay / two);
        let s1 = T::sin(yay / two);
        let c2 = T::cos(pitch / two);
        let s2 = T::sin(pitch / two);
        let c3 = T::cos(roll / two);
        let s3 = T::sin(roll / two);
        let c1c2 = c1 * c2;
        let s1s2 = s1 * s2;

        let q0 = c1c2 * c3 - s1s2 * s3;
        let q1 = c1c2 * s3 + s1s2 * c3;
        let q2 = s1 * c2 * c3 + c1 * s2 * s3;
        let q3 = c1 * s2 * c3 - s1 * c2 * s3;

        Self::new_from(q0, q1, q2, q3)
    }

    /// get the angle of representation from this Quaternion
    pub fn get_angle(&self) -> T {
        let two = T::from(2.0).unwrap();
        let n = self.q.norm2();

        two * T::atan2(n, self.q0)
    }

    // TODO(elsuizo:2021-04-15): this is with `sin` or `cos` ???
    /// get the axis of rotation from which this Quaternion represent
    pub fn get_axis(&self) -> Option<V3<T>> {
        let qn = self.normalize()?;
        let s = T::sin(qn.get_angle() / T::from(2.0)?);
        (s.abs() > T::epsilon()).then(|| qn.q / s)
    }

    /// combine the two previous methods: `get_axis` and `get_angle`
    pub fn axis_angle(&self) -> (T, Option<V3<T>>) {
        (self.get_angle(), self.get_axis())
    }

    // NOTE(elsuizo:2021-04-14): this implementation is a translation from this
    // webpage from Martin John Baker.
    // http://www.euclideanspace.com/maths/geometry/rotations/conversions/eulerToQuaternion/index.htm
    //
    pub fn to_euler_angles(&self) -> (T, T, T) {
        let one = T::one();
        let two = one + one;
        let q0_2 = self.q0 * self.q0;
        let q1_2 = self.q[0] * self.q[0];
        let q2_2 = self.q[1] * self.q[1];
        let q3_2 = self.q[2] * self.q[2];
        // if the Quaternion is normalized this is one, otherwise is a correction factor
        let unit = q0_2 + q1_2 + q2_2 + q3_2;
        let test = self.q[0] * self.q[1] + self.q[2] * self.q0;
        // TODO(elsuizo:2021-04-14): maybe this has to be 0.499999...
        let sing_number = one / two;
        if test - sing_number * unit > T::epsilon() {
            let yay   = two * T::atan2(self.q[0], self.q0);
            let pitch = T::FRAC_PI_2();
            let roll  = T::zero();
            return (yay, pitch, roll);
        }
        if test + sing_number * unit < T::epsilon() {
            let yay   = -two * T::atan2(self.q[0], self.q0);
            let pitch = -one * T::FRAC_PI_2();
            let roll  = T::zero();
            return (yay, pitch, roll);
        }
        let aux1  = two * self.q[1] * self.q0 - two * self.q[0] * self.q[2];
        let aux2  = q1_2 - q2_2 - q3_2 + q0_2;
        let yay   = T::atan2(aux1, aux2);
        let pitch = T::asin(two * test / unit);
        let aux3  = two * self.q[0] * self.q0 - two * self.q[1] * self.q[2];
        let aux4  = -q1_2 + q2_2 - q3_2 + q0_2;
        let roll  = T::atan2(aux3, aux4);
        (yay, pitch, roll)
    }

    fn normalize_q(&self) -> Self {
        let a = self.dot(*self);
        if a > T::epsilon() {
            let mut result = Self::zero();
            result = *self / a.sqrt();
            result.normalized = true;
            result
        } else {
            Self {q0: T::zero(), q: V3::x_axis(), normalized: true}
        }
    }

    fn normalize_a(&self) -> (Self, T) {
        if self.normalized {
            return (*self, T::one())
        }
        let a = self.norm2();
        let mut result = *self / a;
        result.normalized = true;
        return (result, a)
    }

    pub fn argq(&self) -> Self {
        let result = Quaternion::new(T::zero(), self.q);
        result.normalize_q()
    }

    /// exponential function apply to the current Quaternion
    pub fn exp(&self) -> Self {
        let real = self.real();
        let real_exp = T::exp(real);
        let mut scale = real_exp;
        let imag_norm = self.abs_imag();

        if imag_norm > T::epsilon() {
            scale = scale * (T::sin(imag_norm) / imag_norm);
        }

        Self {q0: real_exp * T::cos(imag_norm), q: self.q * scale, normalized: self.norm2() < T::epsilon()}
    }

    /// natural logaritmic function apply to the current Quaternion
    pub fn ln(&self) -> Self {
        let (q_norm, a) = self.normalize_a();
        let real = q_norm.real();
        let mut imag_norm = q_norm.abs_imag();
        let arg_angle = T::atan2(imag_norm, real);
        if imag_norm > T::epsilon() {
            imag_norm = arg_angle / imag_norm;
            Self {q0: T::ln(a), q: q_norm.q * imag_norm, normalized: false}
        } else {
            Self {q0: T::ln(a), q: V3::new_from(arg_angle, T::zero(), T::zero()), normalized: false}
        }
    }

    /// sin function apply to the current Quaternion
    pub fn sin(&self) -> Self {
        let one = T::one();
        let two = one + one;
        let l = self.argq();
        ((*self * l).exp() - (*self * -l).exp())/ (l * two)
    }

    /// cos function apply to the current Quaternion
    pub fn cos(&self) -> Self {
        let one = T::one();
        let two = one + one;
        let l = self.argq();
        ((*self * l).exp() + (*self * -l).exp()) / two
    }

    /// sqrt function apply to the current Quaternion
    pub fn sqrt(&self) -> Self {
        let one = T::one();
        let two = one + one;
        (self.ln() * (one / two)).exp()
    }

    /// power the current Quaternion to the rhs argument
    pub fn pow(&self, rhs: Self) -> Self {
        (rhs * self.ln()).exp()
    }
}
//-------------------------------------------------------------------------
//                        tests
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

    // NOTE(elsuizo:2021-04-14): we assume all the values of the angles in radians!!!
    #[test]
    fn rotate_vec() {
        let q1 = Quaternion::rotation(90.0f32.to_radians(), V3::new_from(0.0, 0.0, 1.0));
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
        let q1 = Quaternion::rotation(90.0f32.to_radians(), V3::new_from(0.0, 0.0, 1.0));
        let x = V3::new_from(1.0, 0.0, 0.0);
        // rotate x around z (90 * 4 = 360) degrees
        let result = q1 * q1 * q1 * q1 * x;
        assert!(compare_floats(result[0], x[0], EPS));
        assert!(compare_floats(result[1], x[1], EPS));
        assert!(compare_floats(result[2], x[2], EPS));
    }

    #[test]
    fn rotate_vec_angle_encode() {
        let q = Quaternion::rotation_norm_encoded(V3::new_from(0.0, 0.0, 90.0f32.to_radians()));
        let x = V3::x_axis();
        let result = q * x;
        let expected = V3::new_from(0.0, -1.0, 0.0);
        assert!(compare_floats(result[0], expected[0], EPS));
        assert!(compare_floats(result[1], expected[1], EPS));
        assert!(compare_floats(result[2], expected[2], EPS));
    }

    #[test]
    fn convert_dcm_test() {
        let q = Quaternion::rotation_norm_encoded(V3::new_from(0.0, 0.0, 90.0f32.to_radians()));
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

    #[test]
    fn inverse_test() {
        let q = Quaternion::new_from(1.0, 1.0, 1.0, 10.0);
        if let Some(inv) = q.inv() {
            let result = q * inv;
            let expected = Quaternion::one();
            assert!(compare_floats(result.q0, expected.q0, EPS));
            assert!(compare_floats(result.q[0], expected.q[0], EPS));
            assert!(compare_floats(result.q[1], expected.q[1], EPS));
            assert!(compare_floats(result.q[2], expected.q[2], EPS));
        }
    }

    #[test]
    fn euler_and_quaternions() {
        let expected = (0.1, 0.2, 0.3);
        let q = Quaternion::from_euler_angles(expected.0, expected.1, expected.2);
        let result = q.to_euler_angles();
        assert!(compare_floats(result.0, expected.0, EPS));
        assert!(compare_floats(result.1, expected.1, EPS));
        assert!(compare_floats(result.2, expected.2, EPS));
    }
}
