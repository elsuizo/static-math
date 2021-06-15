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
use core::fmt;
use core::ops::{Mul, Add, Sub, Neg, Div};
use num::{Num, Float, Signed, Zero, One};
use num::traits::FloatConst;
use crate::vector3::*;
use crate::vector4::V4;
use crate::matrix3x3::M33;
use crate::transformations::rotation_to_euler;
use crate::utils::nearly_zero;
use crate::traits::LinearAlgebra;

/// Quaternion type
#[derive(Copy, Debug, Clone, PartialEq)]
pub struct Quaternion<T> {
    /// Scalar part
    q0: T,
    /// Imaginary part
    q: V3<T>,
    /// flag to signaling if the Quaternion is normalized
    normalized: bool
}

impl<T> Quaternion<T> {

    /// Construct a new Quaternion from a number(real part) and a vector(imag part)
    #[inline(always)]
    pub const fn new(q0: T, q: V3<T>) -> Self {
        Self{q0, q, normalized: false}
    }

    /// Construct a new Quaternion from four numbers
    #[inline(always)]
    pub const fn new_from(q0: T, q1: T, q2: T, q3: T) -> Self {
        Self{q0, q: V3::new([q1, q2, q3]), normalized: false}
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
        Self{q0, q: V3::zeros(), normalized: false}
    }

    /// construct a pure "imaginary" Quaternion
    pub fn new_imag(q: &V3<T>) -> Self {
        Self{q0: T::zero(), q: *q, normalized: false}
    }

    /// calculate the abs2 of the Quaternion
    pub fn abs2(&self) -> T {
        self.q0 * self.q0 + self.q[0] * self.q[0] + self.q[1] * self.q[1] + self.q[2] * self.q[2]
    }

    /// Construct a new Quternion from a V4
    pub fn new_from_vec(v: &V4<T>) -> Self {
        Self{q0: v[0], q: V3::new_from(v[1], v[2], v[3]), normalized: false}
    }

    /// convert the Quaternion to a rotation matrix
    pub fn to_rotation(&self) -> M33<T> {
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
}

// create the zero Quaternion
impl<T: Num + Copy> Zero for Quaternion<T> {
    fn zero() -> Self {
        Self::new(T::zero(), V3::zeros())
    }

    fn is_zero(&self) -> bool {
        *self == Quaternion::zero()
    }
}

// create the unit Quaternion
impl<T: Num + Copy> One for Quaternion<T> {
    /// Create an identity Quaternion
    fn one() -> Self {
        let one = T::one();
        Self{q0: one, q: V3::zeros(), normalized: true}
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
        self * rhs.inverse().expect("the input has to be a non zero vector")
    }
}

// q * q
impl<T: Num + Copy> Mul for Quaternion<T> {
    type Output = Self;

    #[inline(always)]
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
    #[inline(always)]
    fn mul(self, rhs: V3<T>) -> Self::Output {
        let one = T::one();
        let two = one + one;
        let t = (self.q * two).cross(rhs);
        rhs + t * self.q0 + self.q.cross(t)
    }
}

// q * const
impl<T: Num + Copy> Mul<T> for Quaternion<T> {
    type Output = Quaternion<T>;
    fn mul(self, rhs: T) -> Self::Output {
        Self {q0: self.q0 * rhs, q: self.q * rhs, normalized: false}
    }
}

// -q
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

impl<T: Float + core::iter::Sum> Quaternion<T> {

    pub fn from_rotation(rot: &M33<T>) -> Self {
        let m = rot.transpose();
        let zero = T::zero();
        let one  = T::one();
        let half = one / (one + one);
        if m[(2, 2)] < zero {
            if m[(0, 0)] > m[(1, 1)] {
                let t = one + m[(0, 0)] - m[(1, 1)] - m[(2, 2)];
                let q = Self::new_from(m[(1, 2)]-m[(2, 1)],  t,  m[(0, 1)]+m[(1, 0)],  m[(2, 0)]+m[(0, 2)]);
                return q * half / t.sqrt();
            } else {
                let t = one - m[(0, 0)] + m[(1, 1)] - m[(2, 2)];
                let q = Self::new_from(m[(2, 0)]-m[(0, 2)],  m[(0, 1)]+m[(1, 0)],  t,  m[(1, 2)]+m[(2, 1)]);
                return q * half / t.sqrt();
            }
        } else {
            if m[(0, 0)] < -m[(1, 1)] {
                let t = one - m[(0, 0)] - m[(1, 1)] + m[(2, 2)];
                let q = Self::new_from(m[(0, 1)]-m[(1, 0)],  m[(2, 0)]+m[(0, 2)],  m[(1, 2)]+m[(2, 1)],  t);
                return q * half / t.sqrt();
            } else {
                let t = one + m[(0, 0)] + m[(1, 1)] + m[(2, 2)];
                let q = Self::new_from(t,  m[(1, 2)]-m[(2, 1)],  m[(2, 0)]-m[(0, 2)],  m[(0, 1)]-m[(1, 0)]);
                return q * half / t.sqrt();
            }
        }
    }
}

// NOTE(elsuizo:2021-04-23): we only need the Float for the sqrt function
impl<T: Float> Quaternion<T> {

    /// the euclidean norm of the Quaternion
    pub fn norm2(&self) -> T {
        self.dot(*self).sqrt()
    }

    /// normalize the Quaternion only if necessary
    pub fn normalize(&self) -> Option<Self> {
        if self.normalized {
            Some(*self)
        } else {
            let norm_sqr = self.norm2();
            if !nearly_zero(norm_sqr) {
                let mut result = *self / norm_sqr;
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

    /// generate a Quaternion that represents a rotation of a angle `theta`
    /// around the axis(normalized) `v`
    pub fn rotation(theta: T, v: &V3<T>) -> Self {
        let one = T::one();
        let two = one + one;
        let n   = v.normalize().expect("the input has to be a non zero vector");
        let (s, c) = (theta / two).sin_cos();
        let q0  = c;
        let q   = n * s;
        Self{q0, q, normalized: true}
    }

    /// generate a Quaternion that represents a rotation of a angle `theta`
    /// around the axis(normalized) `v`, the angle `theta` is encoded in the
    /// norm of the vector `v`
    pub fn rotation_norm_encoded(v: &V3<T>) -> Self {
        let one = T::one();
        let two = T::from(2.0).unwrap();
        let theta = v.norm2();
        if !nearly_zero(theta) {
            let (s, c) = (theta / two).sin_cos();
            Self{q0: c, q: *v * (s / theta), normalized: true}
        } else {
            Self::new(one, V3::zeros())
        }
    }

    /// create a quaternion that represents the rotation from a Euler angles
    /// with the roll-pitch-yay convention
    pub fn from_euler_angles(yay: T, pitch: T, roll: T) -> Self {
        let one = T::one();
        let two = one + one;
        let (roll_sin, roll_cos)   = (roll / two).sin_cos();
        let (pitch_sin, pitch_cos) = (pitch / two).sin_cos();
        let (yay_sin, yay_cos)     = (yay / two).sin_cos();
        let q0 = roll_cos * pitch_cos * yay_cos + roll_sin * pitch_sin * yay_sin;
        let q1 = roll_sin * pitch_cos * yay_cos - roll_cos * pitch_sin * yay_sin;
        let q2 = roll_cos * pitch_sin * yay_cos + roll_sin * pitch_cos * yay_sin;
        let q3 = roll_cos * pitch_cos * yay_sin - roll_sin * pitch_sin * yay_cos;

        Self{q0, q: V3::new_from(q1, q2, q3), normalized: true}
    }

    /// get the angle of representation from this Quaternion
    pub fn get_angle(&self) -> T {
        let one = T::one();
        let two = one + one;
        let n = self.q.norm2();

        two * T::atan2(n, self.q0)
    }

    /// get the axis of rotation from which this Quaternion represent
    pub fn get_axis(&self) -> Option<V3<T>> {
        let qn = self.normalize()?;
        let s = T::sin(qn.get_angle() / T::from(2.0)?);
        (s.abs() > T::epsilon()).then(|| qn.q / s)
    }

    /// combine the two previous methods: `get_axis` and `get_angle`
    pub fn axis_angle(&self) -> (Option<V3<T>>, T) {
        (self.get_axis(), self.get_angle())
    }

    // TODO(elsuizo:2021-05-20): this epsilon comparison could be wrong maybe we need a
    // nearly_equal here
    /// normalize the Quaternion
    pub fn normalize_q(&self) -> Self {
        let a = self.dot(*self);
        if a > T::epsilon() {
            let mut result = *self / a.sqrt();
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

    /// get the argument of the Quaternion
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

    // TODO(elsuizo:2021-04-24): maybe here its better a error for the corner cases

    /// Brief.
    ///
    /// Spherical Linear Interpolation between two Quaternions
    /// this implementation follow this implementations:
    /// https://www.mrpt.org/tutorials/programming/maths-and-geometry/slerp-interpolation/
    ///
    /// Function arguments:
    /// `a`: Quaternion(normalized)
    ///
    /// `b`: Quaternion(normalized)
    ///
    /// `t`: Float in the closed interval [0.0, 1.0]
    ///
    pub fn slerp(a: Self, b: Self, t: T) -> Self {
        let one = T::one();
        // calculate the angle betwen two unit Quaternions via dot product
        let mut cos_half_theta = a.dot(b);
        // if a = b or a = -b then theta(the angle between) = 0 then we can return a
        if cos_half_theta.abs() >= one {
            return a;
        }
        let mut reverse_a = false;
        // allways follow the shortest path
        if cos_half_theta < T::zero() {
            reverse_a = true;
            cos_half_theta = -cos_half_theta;
        }
        let half_theta     = T::acos(cos_half_theta);
        let sin_half_theta = T::sqrt(one - cos_half_theta * cos_half_theta);
        // TODO(elsuizo:2021-04-24): maybe here the comparison could be with epsilon
        if sin_half_theta.abs() < T::from(0.001).unwrap() {
            if !reverse_a {
                let mut result = Quaternion::zero();
                result.q0   = (one - t) * a.q0 + t * b.q0;
                result.q[0] = (one - t) * a.q[0] + t * b.q[0];
                result.q[1] = (one - t) * a.q[1] + t * b.q[2];
                result.q[2] = (one - t) * a.q[2] + t * b.q[1];
                return result
            }
        }
        let aux1 = T::sin((one - t) * half_theta) / sin_half_theta;
        let aux2 = T::sin(t * half_theta) / sin_half_theta;

        // this part handle the correct orientation
        if !reverse_a {
            let mut result = Quaternion::zero();
            result.q0   = aux1 * a.q0   + aux2 * b.q0;
            result.q[0] = aux1 * a.q[0] + aux2 * b.q[0];
            result.q[1] = aux1 * a.q[1] + aux2 * b.q[2];
            result.q[2] = aux1 * a.q[2] + aux2 * b.q[1];
            return result
        } else {
            let mut result = Quaternion::zero();
            result.q0   = aux1 * a.q0   - aux2 * b.q0;
            result.q[0] = aux1 * a.q[0] - aux2 * b.q[0];
            result.q[1] = aux1 * a.q[1] - aux2 * b.q[2];
            result.q[2] = aux1 * a.q[2] - aux2 * b.q[1];
            return result
        }
    }

    /// Brief.
    ///
    /// calculate the instantaneous Quaternion derivative representing a Quaternion rotating at
    /// rate given by a vector rate
    ///
    /// Function arguments:
    /// `rate`: V3<Float>
    ///
    pub fn derivative(&self, rate: &V3<T>) -> Self {
        let one = T::one();
        let two = one + one;
        Self::new_imag(rate) * (one / two) * (*self)
    }
}

impl<T: Float + Signed> Quaternion<T> {
    /// Calculate the inverse of the Quaternion
    pub fn inverse(&self) -> Option<Self> {
        if !self.normalized {
            let norm_sqr = self.abs2();
            if !nearly_zero(norm_sqr) {
                Some(self.conj() / norm_sqr)
            } else {
                None
            }
        } else {
            Some(self.conj())
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
}

impl<T: Float + FloatConst> Quaternion<T> {
    /// get the euler angles from the Quaternion
    pub fn to_euler_angles(&self) -> (T, T, T) {
        rotation_to_euler(&self.to_rotation())
    }
}

// convert from array to Quaternion
impl<T: Copy> From<[T; 4]> for Quaternion<T> {
    fn from(data: [T; 4]) -> Quaternion<T> {
        Quaternion::new_from(data[0], data[1], data[2], data[3])
    }
}

//-------------------------------------------------------------------------
//                        Display for Quaternion
//-------------------------------------------------------------------------
impl<T: Num + fmt::Display> fmt::Display for Quaternion<T> {
    fn fmt(&self, dest: &mut fmt::Formatter) -> fmt::Result {
        write!(dest, "q0: {0:^3.2}, q:{1:^3.2}", self.q0, self.q)
    }
}

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod test_quaternion {
    use crate::vector3::V3;
    use crate::quaternion::Quaternion;
    use crate::utils::{nearly_equal};

    // NOTE(elsuizo:2021-04-23): this could be more small but the rotation accumulates error in
    // sucesives runs
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
        let q1 = Quaternion::rotation(90.0f32.to_radians(), &V3::new_from(0.0, 0.0, 1.0));
        let x = V3::new_from(1.0, 0.0, 0.0);
        // rotate x around z 90 degrees
        let result = q1 * x;
        let expected = V3::new_from(0.0, 1.0, 0.0);
        assert!(nearly_equal(result[0], expected[0], EPS));
        assert!(nearly_equal(result[1], expected[1], EPS));
        assert!(nearly_equal(result[2], expected[2], EPS));
    }

    #[test]
    fn rotate_vec_composition_360() {
        let q1 = Quaternion::rotation(90.0f32.to_radians(), &V3::new_from(0.0, 0.0, 1.0));
        let x = V3::new_from(1.0, 0.0, 0.0);
        // rotate x around z (90 * 4 = 360) degrees
        let result = q1 * q1 * q1 * q1 * x;
        assert!(nearly_equal(result[0], x[0], EPS));
        assert!(nearly_equal(result[1], x[1], EPS));
        assert!(nearly_equal(result[2], x[2], EPS));
    }

    #[test]
    fn rotate_vec_angle_encode() {
        let q = Quaternion::rotation_norm_encoded(&V3::new_from(0.0, 0.0, 90.0f32.to_radians()));
        let x = V3::x_axis();
        let result = q * x;
        let expected = V3::new_from(0.0, 1.0, 0.0);
        assert!(nearly_equal(result[0], expected[0], EPS));
        assert!(nearly_equal(result[1], expected[1], EPS));
        assert!(nearly_equal(result[2], expected[2], EPS));
    }

    #[test]
    fn convert_rotation_test() {
        let q = Quaternion::rotation_norm_encoded(&V3::new_from(0.0, 0.0, 90.0f32.to_radians()));
        let x = V3::x_axis();
        // rotate the x around z axis 360 degrees
        let expected = q * q * q * q * x;
        // convert the quaternion to a rotation matrix
        let m = q.to_rotation();
        // rotate the x around z axis 360 degrees with the rotation matrix
        let result = m * m * m * m * x;

        assert!(nearly_equal(result[0], expected[0], EPS));
        assert!(nearly_equal(result[1], expected[1], EPS));
        assert!(nearly_equal(result[2], expected[2], EPS));
    }

    #[test]
    fn inverse_test() {
        let q = Quaternion::new_from(1.0, 1.0, 1.0, 10.0);
        if let Some(inv) = q.inverse() {
            let result = q * inv;
            let expected = Quaternion::one();
            assert!(nearly_equal(result.q0, expected.q0, EPS));
            assert!(nearly_equal(result.q[0], expected.q[0], EPS));
            assert!(nearly_equal(result.q[1], expected.q[1], EPS));
            assert!(nearly_equal(result.q[2], expected.q[2], EPS));
        }
    }

    #[test]
    fn euler_and_quaternions() {
        let expected = (0.1, 0.2, 0.3);
        let q = Quaternion::from_euler_angles(expected.0, expected.1, expected.2);
        let result = q.to_euler_angles();
        assert!(nearly_equal(result.0, expected.0, EPS));
        assert!(nearly_equal(result.1, expected.1, EPS));
        assert!(nearly_equal(result.2, expected.2, EPS));
    }

    #[test]
    fn slerp_test() {
        let a = Quaternion::rotation(1.78, &V3::new_from(1.0, 2.0, 3.0));
        let b = Quaternion::rotation(1.78, &V3::x_axis());
        let result = Quaternion::slerp(a, b, 0.3);
        // NOTE(elsuizo:2021-04-24): this result is from julia language
        let expected = Quaternion::new_from(0.6995922116669001, 0.42947374679735195, 0.31677365769795535, 0.475160486546933);
        assert!(nearly_equal(result.q0, expected.q0, EPS));
        assert!(nearly_equal(result.q[0], expected.q[0], EPS));
        assert!(nearly_equal(result.q[1], expected.q[1], EPS));
        assert!(nearly_equal(result.q[2], expected.q[2], EPS));
    }
}
