//-------------------------------------------------------------------------
// @file dual_quaternion.rs
//
// @date 05/17/21 10:01:13
// @author Martin Noblia
// @email mnoblia@disroot.org
//
// @brief
//
// @detail
//
// Licence MIT:
// Copyright <2021> <Martin Noblia>
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
//--------------------------------------------------------------------------
use std::fmt;
use crate::quaternion::Quaternion;
use crate::vector4::V4;
use crate::matrix4x4::M44;
use crate::matrix3x3::M33;
use crate::vector3::V3;
use crate::utils::{nearly_zero, nearly_equal};
use crate::transformations::{get_parts_raw, homogeneous_from_rotation};
use std::ops::{Mul, Add, Sub, Neg, Div};
use num::{Num, Float, Signed, Zero, One};
use num::traits::FloatConst;

// TODO(elsuizo:2021-05-17): maybe we need a flag to the norm
#[derive(Copy, Debug, Clone, PartialEq)]
pub struct DualQuaternion<T> {
    /// the real part
    q_real: Quaternion<T>,
    /// the dual part
    q_dual: Quaternion<T>,
    /// normalized flag
    normalized: bool
}

impl<T> DualQuaternion<T> {
    pub const fn new(q_real: Quaternion<T>, q_dual: Quaternion<T>) -> Self {
        Self{q_real, q_dual, normalized: false}
    }

    pub const fn new_from(a: T, b: T, c: T, d: T, e: T, f: T, g: T, h: T) -> Self {
        Self::new(Quaternion::new_from(a, b, c, d), Quaternion::new_from(e, f, g, h))
    }
}

impl<T: Num + Copy> DualQuaternion<T> {
    /// Get the real Quaternion part
    pub fn real(&self) -> Quaternion<T> {
        self.q_real
    }

    /// Get the dual Quaternion part
    pub fn dual(&self) -> Quaternion<T> {
        self.q_dual
    }

    /// Create a pure rotation transformation from a given Quaternion
    pub fn new_from_rotation(r: &Quaternion<T>) -> Self {
        Self{q_real: *r, q_dual: Quaternion::zero(), normalized: true}
    }

    /// Create a Dual Quaternion that represent a point when the real part is a unit and the dual
    /// part is a pure Quaternion
    pub fn new_from_point(v: &V3<T>) -> Self {
        Self{q_real: Quaternion::one(), q_dual: Quaternion::new_imag(v), normalized: true}
    }

}

// dq + dq
impl<T: Num + Copy> Add for DualQuaternion<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.q_real + rhs.q_real, self.q_dual + rhs.q_dual)
    }
}

// dq - dq
impl<T: Num + Copy> Sub for DualQuaternion<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.q_real - rhs.q_real, self.q_dual - rhs.q_dual)
    }
}

// -dq
impl<T: Num + Copy + Signed> Neg for DualQuaternion<T> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.q_real, -self.q_dual)
    }
}

// dq * dq
impl<T: Num + Copy> Mul for DualQuaternion<T> {
    type Output = Self;

    #[inline(always)]
    fn mul(self, rhs: Self) -> Self::Output {
        let q_real = self.q_real * rhs.q_real;
        let q_dual = self.q_real * rhs.q_dual + self.q_dual * rhs.q_real;
        Self::new(q_real, q_dual)
    }
}

// dq * constant
impl<T: Num + Copy> Mul<T> for DualQuaternion<T> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self::Output {
        Self::new(self.q_real * rhs, self.q_dual * rhs)
    }
}

// dq / dq
impl<T: Float + Signed> Div for DualQuaternion<T> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let rhs_real_sqr = rhs.q_real * rhs.q_real;
        let prod_real    = self.q_real * rhs.q_real / rhs_real_sqr;
        let prod_dual    = (rhs.q_real * self.q_dual - self.q_real * rhs.q_dual) / rhs_real_sqr;
        Self::new(prod_real, prod_dual)
    }
}

impl<T: Float + Signed> DualQuaternion<T> {
    /// Convert the `DualQuaternion` to a homogeneus transformation matrix `M44<Float>` in SE(3)
    pub fn to_homogeneous(&self) -> M44<T> {
        let (r, p) = self.to_rotation_translation();
        homogeneous_from_rotation(&r, &p)
    }

    /// Calculate the inverse of the `DualQuaternion`
    pub fn inverse(&self) -> Self {
        if self.normalized {
            self.conj()
        } else {
            let q_real_inv = self.q_real.inverse().expect("the zero Quaternion cannot be inverted");
            Self::new(q_real_inv, -q_real_inv * self.q_dual * q_real_inv)
        }
    }

    /// Get the screw parameters from this `DualQuaternion` the implementation follow the following
    /// paper: "Dual Quaternions" by Yan-Bin Jia
    ///
    /// Output:
    /// l: V3<Float> a unit 3d vector that represent one of the plucker coordinates
    /// m: V3<Float> a vector in 3d that represent the moment of the line l and his norm represent
    /// the distance from the origin to the line
    /// theta: Float the amount of rotation around the screw axis l
    /// d: Float the amount of translation along the screw axis l
    ///
    pub fn get_screw_parameters(&self) -> (V3<T>, V3<T>, T, T) {
        let one = T::one();
        let two = one + one;
        let half = one / two;
        let theta = self.real().get_angle();
        let t = self.get_translation();

        if !nearly_zero(theta) {
            let l = self.real().imag() / T::sin(theta / two);
            let d = t * l;
            let cot = one / T::tan(theta / two);
            let m = (t.cross(l) + (t - l * d) * cot) * half;
            return (l, m, theta, d)
        } else {
            let d = t.norm2();
            let mut l = V3::zero();
            if !nearly_zero(d) {
                l = t / d;
            }
            let inf = T::infinity();
            let m = V3::new_from(inf, inf, inf);
            return (l, m, theta, d)
        }
    }

    // TODO(elsuizo:2021-05-23): maybe here we need to validate the norm of l
    /// Create a `DualQuaternion` from the screw parameters
    ///
    /// Function arguments:
    /// `l`: a unit vector that represent the screw axis
    /// `m`: screw axis moment that its perpendicular to l (m * l = 0) and the norm represent the
    /// actual moment
    /// `theta`: screw angle, represent the amount of rotation around the screw axis
    /// `d`: screw displacement, represent the amount of displacement along the screw axis
    ///
    pub fn new_from_screw_parameters(l: &V3<T>, m: &V3<T>, theta: T, d: T) -> Self {
        let one = T::one();
        let two = one + one;
        let (s, c) = (theta / two).sin_cos();
        let q_real = Quaternion::new(c, *l * s);
        let q_dual = Quaternion::new(s * (-d / two), *m * s + *l * (d / two) * c);
        Self::new(q_real, q_dual)
    }

    // NOTE(elsuizo:2021-05-23): the multiplication looks funky because we dont do products with
    // constants to the lef size of the multiplications

    /// Compute the power of the `DualQuaternion`
    ///
    /// Function arguments:
    /// `exponent`: Float the exponent to raise the power
    ///
    /// Output: Self^(exponent)
    ///
    pub fn pow(&self, exponent: T) -> Self {
        let one = T::one();
        let two = one + one;
        let theta = two * T::acos(self.real().real());
        let (s, c) = (theta / two).sin_cos();
        if nearly_zero(theta) {
            Self::new_from_point(&(self.get_translation() * exponent))
        } else {
            let s0 = self.real().imag() / s;
            let d  = self.dual().real() * -two / s;
            let se = (self.dual().imag() - s0 * (d / two) * c) / s;

            let (s_exp, c_exp) = (exponent * theta / two).sin_cos();
            let q_real = Quaternion::new(c_exp, s0 * s_exp);
            let q_dual = Quaternion::new(s_exp * -exponent * d / two, s0 * c_exp * exponent * d / two + se * s_exp);
            Self::new(q_real, q_dual)
        }
    }
    // TODO(elsuizo:2021-05-23): maybe this clone is not necesary...

    /// Screw Linear Interpolation: is an extension of the spherical linear interpolation (SLERP)
    /// over Quaternions. It performs a screw motion with constant rotation and translation speeds
    /// from `begin` to `end`
    ///
    /// Function arguments:
    /// `begin`: `DualQuaternion<Float>` the first "point" for the interpolation
    /// `end`: `DualQuaternion<Float>` the second "point" for the interpolation
    /// `tau`: Float in [0, 1] representing how far along and around the screw axis to interpolate
    ///
    pub fn screw_lerp(begin: &Self, end: &Self, tau: T) -> Self {
        let one = T::one();
        let mut start = begin.clone();
        // TODO(elsuizo:2021-05-23): this is from the python implementation that refers to a paper
        // that "ensure we always find closest solution, See Kavan and Zara 2005"
        if (start.real() * end.real()).real() < T::zero() {
            start.q_real = begin.real() * -one;
        }
        start * (start.inverse() * *end).pow(tau)
    }
}

impl<T: Float> DualQuaternion<T> {

    /// Create a `DualQuaternion` from an array that encodes a Quaternion and a 3d vecto (V3<Float>)
    pub fn new_from_array(array: [T; 7]) -> Self {
        let one = T::one();
        let half = one / (one + one);
        let v_q = V4::new_from(array[0], array[1], array[2], array[3]);
        let q_real = Quaternion::new_from_vec(&v_q);
        let v_d = V3::new_from(array[4], array[5], array[6]);
        let q_dual = Quaternion::new_imag(&v_d) * half * q_real;
        Self::new(q_real, q_dual)
    }

    /// Create a `DualQuaternion` that represent a pure translation transformation
    pub fn new_from_translation(t: &V3<T>) -> Self {
        let one = T::one();
        let two = one + one;
        Self{q_real: Quaternion::one(), q_dual: Quaternion::new_imag(&V3::new_from(t[0]/two, t[1]/two, t[2]/two)), normalized: true}
    }

    /// Create a new DualQuaternion from a rotation(Quaternion) and translation(V3)
    pub fn new_from_rot_trans(rot: &Quaternion<T>, translation: &V3<T>) -> Self {
        let one = T::one();
        let half = one / (one + one);
        let q_real = rot.normalize().expect("the quaternion it can't be zero!!!");
        let q_dual = (Quaternion::new_imag(translation) * half) * q_real;
        Self::new(q_real, q_dual)
    }

    pub fn is_normalized(&self) -> bool {
        if self.normalized {
            return true;
        } else {
            if nearly_zero(self.real().norm2()) {
                return true;
            }
            let check1 = nearly_equal(self.real().norm2(), T::one(), T::epsilon());
            let dual_norm = self.dual() / self.real().norm2();
            let check2 = nearly_equal(dual_norm.real(), self.dual().real(), T::epsilon()) &&
                         nearly_equal(dual_norm.imag()[0], self.dual().imag()[0], T::epsilon()) &&
                         nearly_equal(dual_norm.imag()[1], self.dual().imag()[1], T::epsilon()) &&
                         nearly_equal(dual_norm.imag()[2], self.dual().imag()[2], T::epsilon());
            return check1 && check2;
        }
    }

}

impl<T: Float + FloatConst + std::iter::Sum> DualQuaternion<T> {

    /// Normalize the `DualQuaternion` only if necessary
    pub fn normalize(&self) -> Option<Self> {
        if self.normalized {
            Some(*self)
        } else {
            let norm_q_real = self.q_real.norm2();
            if !nearly_zero(norm_q_real) {
                let mut result = Self::one();
                result.q_real = self.q_real / norm_q_real;
                result.q_dual = self.q_dual / norm_q_real;
                result.normalized = true;
                Some(result)
            } else {
                None
            }
        }
    }

    // TODO(elsuizo:2021-05-18): maybe here is better with the function `get_parts` that get the
    // Quaternion from euler angles
    /// Create a new `DualQuaternion` from a homogeneous transformation matrix in SE(3)
    pub fn new_from_homogeneous(t: &M44<T>) -> Self {
        let (rot, p) = get_parts_raw(t);
        let q = Quaternion::from_rotation(&rot).normalize().expect("the quaternion it can't be zero!!!");
        Self::new_from_rot_trans(&q, &p)
    }

    /// Create a `DualQuaternion` that represent a rotation pure transformation
    pub fn new_from_rotation_matrix(m: &M33<T>) -> Self {
        Self{q_real: Quaternion::from_rotation(m), q_dual: Quaternion::zero(), normalized: true}
    }
}

impl<T: Num + Copy + Signed> DualQuaternion<T> {
    /// Calculate the conjugate of the `DualQuaternion`
    pub fn conj(&self) -> Self {
        Self::new(self.q_real.conj(), self.q_dual.conj())
    }

    /// Calculate the dual number conjugation of the `DualQuaternion`
    pub fn conj_as_dual(&self) -> Self {
        Self::new(self.q_real, -self.q_dual)
    }

    /// Calculate the combined(as a dual number and Quaternion) conjugation of the `DualQuaternion`
    pub fn conj_combined(&self) -> Self {
        Self::new(self.q_real.conj(), -self.q_dual.conj())
    }

    /// Calculate the norm of the DualQuaternion
    pub fn norm(&self) -> Self {
        *self * self.conj()
    }

    /// Get the underlying rotation and translation from the `DualQuaternion`
    pub fn to_rotation_translation(&self) -> (M33<T>, V3<T>) {
        let r = self.real().to_rotation();
        let t = self.get_translation();
        (r, t)
    }

    /// Get the underlying translation from the `DualQuaternion`
    pub fn get_translation(&self) -> V3<T> {
        let one = T::one();
        let two = one + one;
        ((self.dual() * two) * self.real().conj()).imag()
    }

    /// Get the underlying Quaternion and translation from the `DualQuaternion`
    pub fn to_quaternion_translation(&self) -> (Quaternion<T>, V3<T>) {
        let one = T::one();
        let two = one + one;
        let q = self.real();
        let t = (self.dual() * two) * self.real().conj();
        (q, t.imag())
    }

    /// Transform the given point in 3D with the transformation that represent this
    /// `DualQuaternion` like a homogeneous transformation in SE(3)
    pub fn transform_point(&self, point: &V3<T>) -> V3<T> {
        let dq_point = Self::new_from_point(point);
        (*self * dq_point * self.conj()).dual().imag()
    }
}

impl<T: Num + Copy> DualQuaternion<T> {

    /// Construct a new DualQuaternion from two V4
    pub fn new_from_vecs(q_real: &V4<T>, q_dual: &V4<T>) -> Self {
        Self::new(Quaternion::new_from_vec(q_real), Quaternion::new_from_vec(q_dual))
    }

    /// construct a zero DualQuaternion
    pub fn zero() -> Self {
        <DualQuaternion<T> as Zero>::zero()
    }

    /// construct a unit DualQuaternion
    pub fn one() -> DualQuaternion<T> {
        <DualQuaternion<T> as One>::one()
    }
}

// create the zero DualQuaternion
impl<T: Num + Copy> Zero for DualQuaternion<T> {
    fn zero() -> Self {
        Self::new(Quaternion::zero(), Quaternion::zero())
    }

    fn is_zero(&self) -> bool {
        *self == DualQuaternion::zero()
    }
}

// create the unit DualQuaternion
impl<T: Num + Copy> One for DualQuaternion<T> {
    /// Create an identity DualQuaternion
    fn one() -> Self {
        Self{q_real: Quaternion::one(), q_dual: Quaternion::zero(), normalized: true}
    }
}

// TODO(elsuizo:2021-05-17): this could be better...
//-------------------------------------------------------------------------
//                      Display for DualQuaternion
//-------------------------------------------------------------------------
impl<T: Num + fmt::Display> fmt::Display for DualQuaternion<T> {
    fn fmt(&self, dest: &mut fmt::Formatter) -> fmt::Result {
        println!();
        write!(dest, "real: {}\ndual: {}", self.q_real, self.q_dual)
    }
}

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod test_dual_quaternion {
    use crate::dual_quaternion::DualQuaternion;
    use crate::quaternion::Quaternion;
    use crate::vector3::V3;
    use crate::vector4::V4;
    use crate::matrix3x3::M33;
    use crate::matrix4x4::M44;
    use crate::m44_new;
    use crate::utils::{nearly_equal, nearly_zero, compare_vecs, compare_dual_quaternions};
    use crate::transformations::{homogeneous_from_quaternion, euler_to_rotation};
    const EPS: f32 = 1e-4;

    #[test]
    fn basic_dual_quaternion_tests() {
        let dq1 = DualQuaternion::new_from(0.7071067811, 0.7071067811, 0.0, 0.0, -3.535533905, 3.535533905, 1.767766952, -1.767766952);
        let dq2 = DualQuaternion::new_from_vecs(&V4::new_from(0.7071067811, 0.7071067811, 0.0, 0.0), &V4::new_from(-3.535533905, 3.535533905, 1.767766952, -1.767766952));

        assert!(compare_dual_quaternions(dq1, dq2, EPS))
    }

    #[test]
    fn unity_product_dual_quaternion_test() {
        let unit = DualQuaternion::one();
        let expected: DualQuaternion<f32> = DualQuaternion::new_from_translation(&V3::z_axis());
        let result = unit * expected;
        assert!(compare_dual_quaternions(result, expected, EPS))
    }

    #[test]
    fn dual_quaternion_creation_tests() {
        // create a new DualQuaternion that represents a translation pure transformation
        let dq: DualQuaternion<f32> = DualQuaternion::new_from_translation(&V3::z_axis());
        // get the vector of translation and the rotation matrix
        let result = dq.to_rotation_translation();
        let result2 = dq.to_quaternion_translation();
        let expected1 = V3::z_axis();

        assert_eq!(
            &result.1[..],
            &expected1[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result.1[..],
            &expected1[..]
        );

        assert_eq!(
            &result2.1[..],
            &expected1[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result2.1[..],
            &expected1[..]
        );

        let expected2 = M33::identity();
        assert!(compare_vecs(&result.0.as_vec(), &expected2.as_vec(), EPS));

        let expected3 = Quaternion::one();
        assert!(nearly_equal(result2.0.real(), expected3.real(), EPS));
        assert!(nearly_equal(result2.0.imag()[0], expected3.imag()[0], EPS));
        assert!(nearly_equal(result2.0.imag()[1], expected3.imag()[1], EPS));
        assert!(nearly_equal(result2.0.imag()[2], expected3.imag()[2], EPS));
    }

    #[test]
    fn dual_quaternion_test_norm() {
        // let q = Quaternion::rotation(1.78, &V3::x_axis());
        // let dq = DualQuaternion::new_from_rot_trans(&q, &V3::y_axis());
        let dq = DualQuaternion::new_from_vecs(&V4::new_from(1.0, 2.0, 3.0, 4.0), &V4::new_from(5.0, 6.0, 7.0, 8.0));
        let one:DualQuaternion<f32> = DualQuaternion::one();
        let normalized = dq.normalize().unwrap();
        assert!(normalized.is_normalized());
        assert!(one.is_normalized());
        assert!(compare_dual_quaternions(one, one.normalize().unwrap(), EPS))
    }

    // NOTE(elsuizo:2021-05-18): homogeneous ---> DualQuaternion ---> homogeneous
    #[test]
    fn homogeneous_conversions() {
        let q  = Quaternion::from_euler_angles(10f32.to_radians(), 10f32.to_radians(), 10f32.to_radians());
        let expected = homogeneous_from_quaternion(&q, &V3::new_from(1.0, 2.0, 3.0));
        let result = DualQuaternion::new_from_homogeneous(&expected).to_homogeneous();

        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn conjugate_product_test() {
        let q1  = Quaternion::from_euler_angles(10f32.to_radians(), 10f32.to_radians(), 10f32.to_radians());
        let q2  = Quaternion::from_euler_angles(30f32.to_radians(), 10f32.to_radians(), 30f32.to_radians());
        let dq1 = DualQuaternion::new_from_rot_trans(&q1, &V3::x_axis());
        let dq2 = DualQuaternion::new_from_rot_trans(&q2, &V3::z_axis());

        let result = (dq1 * dq2).conj();
        let expected = dq2.conj() * dq1.conj();

        assert!(compare_dual_quaternions(result, expected, EPS));
    }

    // NOTE(elsuizo:2021-05-19): rot_pure * traslation_pure == normal_dual
    #[test]
    fn combined_transformations_tests() {
        let rot = euler_to_rotation(10f32.to_radians(), 10f32.to_radians(), 10f32.to_radians(), None);
        let q  = Quaternion::from_euler_angles(10f32.to_radians(), 10f32.to_radians(), 10f32.to_radians());

        let t_pure = DualQuaternion::new_from_translation(&V3::x_axis());
        let r_pure = DualQuaternion::new_from_rotation(&q);

        let expected = DualQuaternion::new_from_rot_trans(&q, &V3::x_axis());
        let result1 = t_pure * r_pure;

        let r_pure2 = DualQuaternion::new_from_rotation_matrix(&rot);
        let result2 = t_pure * r_pure2;

        assert!(compare_dual_quaternions(result1, expected, EPS));
        assert!(compare_dual_quaternions(result2, expected, EPS));
    }

    // NOTE(elsuizo:2021-05-20): if the DualQuaternion is not normalized the accumulate rounding
    // error increase and we need a EPS = 1e-6
    #[test]
    fn dual_quaternion_inverse_test() {
        // let dq = DualQuaternion::new_from_vecs(&V4::new_from(1.0, 2.0, 3.0, 4.0), &V4::new_from(5.0, 6.0, 7.0, 8.0));
        let q  = Quaternion::from_euler_angles(10f32.to_radians(), 10f32.to_radians(), 10f32.to_radians());
        let dq = DualQuaternion::new_from_rot_trans(&q, &V3::x_axis());
        let result = dq.inverse() * dq;
        let expected = DualQuaternion::one();
        assert!(compare_dual_quaternions(result, expected, EPS));
    }

    // NOTE(elsuizo:2021-05-20): the inverse of the transformation represented by a DualQuaternion
    // is a simple conjugation!!!
    #[test]
    fn dual_quaternion_transformation_inverse() {
        let t_1_2 = m44_new!(0.0, 1.0, 0.0, 2.0;
                            -1.0, 0.0, 0.0, 4.0;
                             0.0, 0.0, 1.0, 6.0;
                             0.0, 0.0, 0.0, 1.0);

        let t_2_1 = m44_new!(0.0, -1.0,  0.0,  4.0;
                             1.0,  0.0,  0.0, -2.0;
                             0.0,  0.0,  1.0, -6.0;
                             0.0,  0.0,  0.0,  1.0);

        let dq_1_2 = DualQuaternion::new_from_homogeneous(&t_1_2);
        let dq_2_1 = DualQuaternion::new_from_homogeneous(&t_2_1);

        assert!(compare_vecs(&dq_2_1.to_homogeneous().as_vec(), &dq_1_2.inverse().to_homogeneous().as_vec(), EPS));
        assert!(compare_vecs(&dq_1_2.to_homogeneous().as_vec(), &dq_2_1.inverse().to_homogeneous().as_vec(), EPS));
    }

    #[test]
    fn dual_quaternion_transform_point_test() {
        let q  = Quaternion::from_euler_angles(10f32.to_radians(), 10f32.to_radians(), 10f32.to_radians());
        let t = homogeneous_from_quaternion(&q, &V3::new_from(1.0, 2.0, 3.0));

        let p = V4::new_from(1.0, 2.0, 3.0, 0.0);
        let expected = t * p;

        let p = V3::new_from(1.0, 2.0, 3.0);
        let dq = DualQuaternion::new_from_homogeneous(&t);
        let result = dq.transform_point(&p);

        // NOTE(elsuizo:2021-05-22): the transform point with the DualQuaternion should be the same
        // that the with the homogeneus transformation(the first three elements)
        assert!(nearly_equal(result[0], expected[0], EPS));
        assert!(nearly_equal(result[1], expected[1], EPS));
        assert!(nearly_equal(result[2], expected[2], EPS));
        assert!(nearly_zero(expected[3]));
    }

    #[test]
    fn dual_quaternion_get_screw_params_translation_test() {
        use num::Float;
        let trans = &V3::new_from(10.0, 3.7, 7.3);
        let t_pure = DualQuaternion::new_from_translation(&trans);
        let (l_result, m_result, theta_result, d_result) = t_pure.get_screw_parameters();
        let inf = f32::infinity();
        // NOTE(elsuizo:2021-05-22): l should be unit vector
        assert!(nearly_equal(l_result.norm2(), 1.0, EPS));
        // NOTE(elsuizo:2021-05-22): we use inf to signaling that m is tooo far
        assert!(compare_vecs(&*m_result, &*V3::new_from(inf, inf, inf), EPS));
        // NOTE(elsuizo:2021-05-22): theta should be zero because no rotation in the transformation
        assert!(nearly_zero(theta_result));
        // NOTE(elsuizo:2021-05-22): d should be finite
        assert!(d_result.is_finite());
        // NOTE(elsuizo:2021-05-22): the amount of distance along the screw axis should be the norm
        // of the translation vector because this is a translation pure transformation
        assert!(nearly_equal(d_result, trans.norm2(), EPS))
    }

    #[test]
    fn dual_quaternion_screw_params_rotation_test() {
        let q  = Quaternion::from_euler_angles(73f32.to_radians(), 10f32.to_radians(), 10f32.to_radians());
        let r_pure = DualQuaternion::new_from_rotation(&q);
        let (l_result, m_result, theta_result, d_result) = r_pure.get_screw_parameters();
        // NOTE(elsuizo:2021-05-22): l should be a unit vector
        assert!(nearly_equal(l_result.norm2(), 1.0, EPS));
        // NOTE(elsuizo:2021-05-22): m should be zero because the moment should be zero
        assert!(compare_vecs(&*m_result, &*V3::zeros(), EPS));
        // NOTE(elsuizo:2021-05-22): theta should be equal to the angle of rotation that represent
        // the Quaternion transformation
        assert!(nearly_equal(theta_result, q.get_angle(), EPS));
        // NOTE(elsuizo:2021-05-22): d should be zero because this is a pure rotation
        // transformation
        assert!(nearly_zero(d_result));
    }

    // NOTE(elsuizo:2021-05-22): rotate around the z axis 45 degrees and translate along the z axis
    // 7 units, this example the l vector should be z_hat, d should be 7, theta should be 45
    // degrees and m should be zero because the moment is zero
    // like a screw in the floor that rotate 45 degrees and go 7 units up :)
    #[test]
    fn dual_quternion_screw_paras_full() {
        let v = V3::z_axis() * 45f32.to_radians();
        let q = Quaternion::rotation_norm_encoded(&v);
        let trans = V3::new_from(0.0, 0.0, 7.0);
        let dq_full = DualQuaternion::new_from_rot_trans(&q, &trans);
        let (l, m, theta, d) = dq_full.get_screw_parameters();
        assert!(compare_vecs(&*l, &*V3::z_axis(), EPS));
        assert!(compare_vecs(&*m, &*V3::zeros(), EPS));
        assert!(nearly_equal(theta, 45f32.to_radians(), EPS));
        assert!(nearly_equal(d, 7.0, EPS));
    }

    #[test]
    fn dual_quaternion_screw_params_test() {
        let v = V3::z_axis() * 45f32.to_radians();
        let q = Quaternion::rotation_norm_encoded(&v);
        let trans = V3::new_from(0.0, 0.0, 7.0);
        let dq_full = DualQuaternion::new_from_rot_trans(&q, &trans);
        let (l, m, theta, d) = dq_full.get_screw_parameters();
        let result = DualQuaternion::new_from_screw_parameters(&l, &m, theta, d);

        assert!(compare_dual_quaternions(dq_full, result, EPS));
    }

    #[test]
    fn dual_quaternion_pow_test() {
        let q  = Quaternion::from_euler_angles(73f32.to_radians(), 10f32.to_radians(), 10f32.to_radians());
        let trans = V3::new_from(0.0, 0.0, 7.0);
        let dq_full = DualQuaternion::new_from_rot_trans(&q, &trans);
        let expected = dq_full * dq_full * dq_full;
        let result = dq_full.pow(3.0);

        assert!(compare_dual_quaternions(expected, result, EPS));
    }

    #[test]
    fn dual_quaternion_interpolation_test() {
        let taus = [0.0, 0.37, 0.73, 1.0];
        let q  = Quaternion::from_euler_angles(73f32.to_radians(), 10f32.to_radians(), 10f32.to_radians());
        let trans = V3::new_from(0.0, 0.0, 7.0);
        let dq = DualQuaternion::new_from_rot_trans(&q, &trans).normalize().unwrap();
        let (l, m, theta, d) = dq.get_screw_parameters();

        for tau in &taus {
            let result = DualQuaternion::screw_lerp(&DualQuaternion::one(), &dq, *tau);
            let expected = DualQuaternion::new_from_screw_parameters(&l, &m, *tau * theta, *tau * d);
            assert!(compare_dual_quaternions(expected, result, EPS));
        }

    }
}
