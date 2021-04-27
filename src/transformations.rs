//-------------------------------------------------------------------------
// @file transformations.rs
//
// @date 09/14/20 09:11:48
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
use crate::matrix2x2::M22;
use crate::matrix3x3::M33;
use crate::matrix4x4::M44;
use crate::m44_new;
use crate::vector3::V3;
use crate::quaternion::Quaternion;

use crate::utils::nearly_equal;
use num::{Float};
use num::traits::FloatConst;

///
/// Euler sequences conventions of rotations
///
pub enum EulerSeq {
    XYX,
    XYZ,
    XZX,
    XZY,
    YXY,
    YXZ,
    YZX,
    YZY,
    ZXY,
    ZXZ
}
//-------------------------------------------------------------------------
//                        transformations
//-------------------------------------------------------------------------
/// Compute rotation matrix from a angle in radians
pub fn rot2<T: Float>(angle: T) -> M22<T> {
    let (s, c) = angle.sin_cos();
    m22_new!(c, -s;
             s,  c)
}

/// brief.
///
/// compute the rotation around the `x` axis(in cartesian coordinates)
///
/// description
///
/// * `angle` - angle of rotation in radians
///
pub fn rotx<T: Float>(angle: T) -> M33<T> {
    let one = T::one();
    let zero = T::zero();
    let (s, c) = angle.sin_cos();
    m33_new!( one, zero, zero;
             zero,    c,   -s;
             zero,    s,    c)
}

/// Brief.
///
/// Compute the rotation around the `y` axis(in cartesian coordinates)
///
/// Description
///
/// * `angle` - Angle of rotation in radians
///
pub fn roty<T: Float>(angle: T) -> M33<T> {
    let one = T::one();
    let zero = T::zero();
    let (s, c) = angle.sin_cos();
    m33_new!(   c, zero,    s;
             zero,  one, zero;
               -s, zero,    c)
}

/// Brief.
///
/// Compute the rotation around the `z` axis(in cartesian coordinates)
///
/// Description
///
/// * `angle` - Angle of rotation in radians
///
pub fn rotz<T: Float>(angle: T) -> M33<T> {
    let one = T::one();
    let zero = T::zero();
    let (s, c) = angle.sin_cos();
    m33_new!(   c,   -s, zero;
                s,    c, zero;
             zero, zero,  one)
}

/// Brief.
///
/// Compute the rotation matrix from euler angles with the following conventions:
/// XYX, XYZ, XZX, XZY, YXY, YXZ, YZX, YZY, ZXY, ZXZ
///
/// Function arguments:
///
/// yay: first euler angle in radians (Float number)
///
/// pitch: second euler angle in radians (Float number)
///
/// roll: third euler angle in radians (Float number)
///
/// s: `Option<EulerSeq>:` Optional Euler sequence if is None compute ZYX
///
/// Output:
/// r: Rotation matrix(`M33<Float>`)
///
pub fn euler_to_rotation<T: Float>(yay: T, pitch: T, roll: T, s: Option<EulerSeq>) -> M33<T> {
    match s {
        Some(EulerSeq::XYX) => rotx(yay) * roty(pitch) * rotx(roll),
        Some(EulerSeq::XYZ) => rotx(yay) * roty(pitch) * rotz(roll),
        Some(EulerSeq::XZX) => rotx(yay) * rotz(pitch) * rotx(roll),
        Some(EulerSeq::XZY) => rotx(yay) * rotz(pitch) * roty(roll),
        Some(EulerSeq::YXY) => roty(yay) * rotx(pitch) * roty(roll),
        Some(EulerSeq::YXZ) => roty(yay) * rotx(pitch) * rotz(roll),
        Some(EulerSeq::YZX) => rotx(yay) * roty(pitch) * rotz(roll),
        Some(EulerSeq::YZY) => roty(yay) * rotz(pitch) * roty(roll),
        Some(EulerSeq::ZXY) => rotz(yay) * rotx(pitch) * roty(roll),
        Some(EulerSeq::ZXZ) => rotz(yay) * rotx(pitch) * rotz(roll),
        None                => rotz(yay) * roty(pitch) * rotx(roll)
    }
}

// TODO(elsuizo:2021-04-23): handle all the rotation sequences
/// Brief.
///
/// get the euler angles from a rotation matrix comming from the convention ZYX
///
/// Function arguments:
///
/// `r`: Rotation matrix(M33<Float>)
///
/// Output:
///
/// Euler angles: (yay, pitch, roll)
///
pub fn rotation_to_euler<T: Float + FloatConst>(r: M33<T>) -> (T, T, T) {
    let one   = T::one();
    let pitch = T::atan2(-r[(2, 0)], (r[(0, 0)] * r[(0, 0)] + r[(1, 0)] * r[(1, 0)]).sqrt());

    // singularity
    if nearly_equal(pitch, -one * T::FRAC_PI_2(), T::epsilon()) {
        let yay  = T::atan2(-r[(1, 2)], -r[(0, 2)]);
        let roll = T::zero();
        (yay, pitch, roll)
    // singularity
    } else if nearly_equal(pitch, T::FRAC_PI_2(), T::epsilon()) {
        let yay  = T::atan2(r[(1, 2)], r[(0, 2)]);
        let roll = T::zero();
        (yay, pitch, roll)
    // normal case
    } else {
        let yay   = T::atan2(r[(1, 0)], r[(0, 0)]);
        let roll  = T::atan2(r[(2, 1)], r[(2, 2)]);
        (yay, pitch, roll)
    }
}

// TODO(elsuizo:2021-04-27): i think that the names are too long...
pub fn homogeneous_from_rotation<T: Float>(r: M33<T>, p: V3<T>) -> M44<T> {
    let zero = T::zero();
    let one  = T::one();
    m44_new!(r[(0, 0)], r[(0, 1)], r[(0, 2)], p[0];
             r[(1, 0)], r[(1, 1)], r[(1, 2)], p[1];
             r[(2, 0)], r[(2, 1)], r[(2, 2)], p[2];
                zero  ,   zero   ,    zero  , one)
}

pub fn homogeneous_from_quaternion<T: Float>(q: Quaternion<T>, p: V3<T>) -> M44<T> {
    homogeneous_from_rotation(q.to_rotation(), p)
}

// TODO(elsuizo:2020-09-20): do the tests
//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod test_transformations {

    use super::{rotation_to_euler, euler_to_rotation, rotx, roty, rotz};
    use crate::utils::{nearly_equal, is_rotation};
    const EPS: f32 = 1e-6;

    #[test]
    fn rotation_and_euler_test() {
        let expected = (0.1, 0.2, 0.3);
        let r = euler_to_rotation(expected.0, expected.1, expected.2, None);
        let result = rotation_to_euler(r);
        assert!(nearly_equal(result.0, expected.0, EPS));
        assert!(nearly_equal(result.1, expected.1, EPS));
        assert!(nearly_equal(result.2, expected.2, EPS));
    }

    #[test]
    fn rotation_x() {
        let r = rotx(20f32.to_radians());
        assert!(is_rotation(r));
    }

    #[test]
    fn rotation_y() {
        let r = roty(20f32.to_radians());
        assert!(is_rotation(r));
    }

    #[test]
    fn rotation_z() {
        let r = rotz(20f32.to_radians());
        assert!(is_rotation(r));
    }

}
