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
use crate::utils::nearly_zero;
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
    pub fn new(q_real: Quaternion<T>, q_dual: Quaternion<T>) -> Self {
        Self{q_real, q_dual, normalized: false}
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

impl<T:Float + Signed> DualQuaternion<T> {

    /// Convert the `DualQuaternion` to a homogeneus transformation matrix `M44<Float>` in SE(3)
    pub fn to_homogeneous(&self) -> M44<T> {
        let (r, p) = self.to_rotation_translation();
        homogeneous_from_rotation(&r, &p)
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
                let mut result = Self::zero();
                result.q_real = self.q_real / norm_q_real;
                result.q_dual = self.q_dual / norm_q_real;
                result.normalized = true;
                Some(result)
            } else {
                None
            }
        }
    }

    /// Create a new DualQuaternion from a rotation(Quaternion) and translation(V3)
    pub fn new_from_rot_trans(rot: &Quaternion<T>, translation: &V3<T>) -> Self {
        let one = T::one();
        let half = one / (one + one);
        let q_real = rot.normalize().expect("the quaternion it can't be zero!!!");
        let q_dual = (Quaternion::new_imag(translation) * half) * q_real;
        Self::new(q_real, q_dual)
    }

    // TODO(elsuizo:2021-05-18): maybe here is better with the function `get_parts` that get the
    // Quaternion from euler angles
    pub fn new_from_homogeneous(t: &M44<T>) -> Self {
        let (rot, p) = get_parts_raw(t);
        let q = Quaternion::from_rotation(&rot).normalize().expect("the quaternion it can't be zero!!!");
        Self::new_from_rot_trans(&q, &p)
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

    pub fn to_rotation_translation(&self) -> (M33<T>, V3<T>) {
        let one = T::one();
        let two = one + one;
        let r = self.q_real.to_rotation();
        let t = (self.q_dual * two) * self.q_real.conj();
        (r, t.imag())
    }

}

impl<T: Num + Copy> DualQuaternion<T> {

    /// Construct a new DualQuaternion from two V4
    pub fn new_from(q_real: &V4<T>, q_dual: &V4<T>) -> Self {
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
        Self::new(Quaternion::one(), Quaternion::zero())
    }
}

// TODO(elsuizo:2021-05-17): this could be better...
//-------------------------------------------------------------------------
//                      Display for DualQuaternion
//-------------------------------------------------------------------------
impl<T: Num + fmt::Display> fmt::Display for DualQuaternion<T> {
    fn fmt(&self, dest: &mut fmt::Formatter) -> fmt::Result {
        println!("");
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
    use crate::utils::{nearly_equal, compare_vecs};
    use crate::transformations::{homogeneous_from_quaternion};
    const EPS: f32 = 1e-7;

    #[test]
    fn test_norm() {
        let q = Quaternion::rotation(1.78, V3::x_axis());
        let dq = DualQuaternion::new_from_rot_trans(&q, &V3::y_axis());
        let normalized = dq.normalize().unwrap();
        let result = normalized.norm();
        let expected: DualQuaternion<f32> = DualQuaternion::one();
        assert!(nearly_equal(result.q_real.real(), expected.q_real.real(), EPS));
        assert!(nearly_equal(result.q_real.imag()[0], expected.q_real.imag()[0], EPS));
        assert!(nearly_equal(result.q_real.imag()[1], expected.q_real.imag()[1], EPS));
        assert!(nearly_equal(result.q_real.imag()[2], expected.q_real.imag()[2], EPS));

        assert!(nearly_equal(result.q_dual.real(),    expected.q_dual.real(), EPS));
        assert!(nearly_equal(result.q_dual.imag()[0], expected.q_dual.imag()[0], EPS));
        assert!(nearly_equal(result.q_dual.imag()[1], expected.q_dual.imag()[1], EPS));
        assert!(nearly_equal(result.q_dual.imag()[2], expected.q_dual.imag()[2], EPS));

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

        assert!(nearly_equal(result.q_real.real(), expected.q_real.real(), EPS));
        assert!(nearly_equal(result.q_real.imag()[0], expected.q_real.imag()[0], EPS));
        assert!(nearly_equal(result.q_real.imag()[1], expected.q_real.imag()[1], EPS));
        assert!(nearly_equal(result.q_real.imag()[2], expected.q_real.imag()[2], EPS));

        assert!(nearly_equal(result.q_dual.real(),    expected.q_dual.real(), EPS));
        assert!(nearly_equal(result.q_dual.imag()[0], expected.q_dual.imag()[0], EPS));
        assert!(nearly_equal(result.q_dual.imag()[1], expected.q_dual.imag()[1], EPS));
        assert!(nearly_equal(result.q_dual.imag()[2], expected.q_dual.imag()[2], EPS));
    }
}
