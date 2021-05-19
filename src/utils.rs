//-------------------------------------------------------------------------
// @file utils.rs
//
// @date 06/02/20 11:06:12
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
use num::{Float};
use crate::matrix3x3::M33;
use crate::traits::LinearAlgebra;
use crate::matrix4x4::M44;
use crate::dual_quaternion::DualQuaternion;
// NOTE(elsuizo:2020-06-02): the following function
// is a translation of the implementation that is here:
// https://floating-point-gui.de/errors/comparison/
//
/// a comparison function for floating point values
///
pub fn nearly_equal<T: Float>(a: T, b: T, epsilon: T) -> bool {
    let abs_a = a.abs();
    let abs_b = b.abs();
    let abs_diff = (a - b).abs();
    let zero = T::zero();
    // short-cut, handles infinity
    if a == b { true }

    else if a == zero || b == zero || (abs_a + abs_b < T::min_value()) {
        // a or b is zero or both are extremely close to it
        // relative error is less meaningful here
        abs_diff < epsilon
    } else {
      abs_diff / T::min(abs_a + abs_b, T::max_value()) < epsilon
    }
}

pub fn nearly_zero<T: Float>(a: T) -> bool {
    nearly_equal(a, T::zero(), T::epsilon())
}

// TODO(elsuizo:2020-06-16): this is the only place when i use a Vec
/// utility function to compare vectors of Floats
pub fn compare_vecs<T: Float>(v1: &[T], v2: &[T], epsilon: T) -> bool {
    let v_result: Vec<bool> = v1
        .iter()
        .zip(v2)
        .map(|(a, b)| nearly_equal(*a, *b, epsilon))
        .collect();
    v_result.iter().all(|&x| x)
}

pub fn compare_floats<T: Float>(num1: T, num2: T, tol: T) -> bool {
    Float::abs(num1 - num2) < tol
}

/// utility function to verify if a Matrix is a propper rotation matrix
pub fn is_rotation<T: Float + std::iter::Sum>(r: M33<T>) -> bool {
    let r2 = r * r;
    if nearly_equal(r.det(), T::one(), T::epsilon()) && nearly_equal(r2.det(), T::one(), T::epsilon()) {
        true
    } else {
        false
    }
}

pub fn is_rotation_h<T: Float + std::iter::Sum>(r: M44<T>) -> bool {
    let r2 = r * r;
    let eps = T::from(1e-6).unwrap();
    if nearly_equal(r.det(), T::one(), eps) && nearly_equal(r2.det(), T::one(), eps) {
        true
    } else {
        false
    }
}

pub fn compare_dual_quaternions<T: Float>(a: DualQuaternion<T>, b: DualQuaternion<T>, epsilon: T) -> bool {
    nearly_equal(a.real().real(), b.real().real(), epsilon) &&
    nearly_equal(a.real().imag()[0], b.real().imag()[0], epsilon) &&
    nearly_equal(a.real().imag()[1], b.real().imag()[1], epsilon) &&
    nearly_equal(a.real().imag()[2], b.real().imag()[2], epsilon) &&

    nearly_equal(a.dual().real(),    b.dual().real(), epsilon) &&
    nearly_equal(a.dual().imag()[0], b.dual().imag()[0], epsilon) &&
    nearly_equal(a.dual().imag()[1], b.dual().imag()[1], epsilon) &&
    nearly_equal(a.dual().imag()[2], b.dual().imag()[2], epsilon)
}

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod test_utils {
    use super::*;
    #[test]
    fn test_nearly_equal() {
        let a = 0.15 + 0.15;
        let b = 0.1 + 0.2;
        assert_eq!(nearly_equal(a, b, 1e-10), true);
    }

    #[test]
    fn test_zero_signs() {
        let a = 0.0;
        let b = -0.0;
        assert_eq!(nearly_equal(a, b, 1e-10), true);
    }

    #[test]
    fn test_one_zero() {
        let b = -0.0000000005561918225744;
        let a = 0.0;
        assert_eq!(nearly_equal(a, b, 1e-8), true);
    }
}
