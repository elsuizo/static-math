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
use crate::vector3::V3;

use num::{Float, Zero};

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
/// Compute rotation matrix from a angle in degrees
pub fn rot2<T: Float>(angle: T) -> M22<T> {
    let c = angle.to_radians().cos();
    let s = angle.to_radians().sin();
    m22_new!(c, -s;
             s,  c)
}

/// brief.
///
/// compute the rotation around the `x` axis(in cartesian coordinates)
///
/// description
///
/// * `angle` - angle of rotation in degrees
///
pub fn rotx<T: Float>(angle: T) -> M33<T> {
    let one = T::one();
    let zero = T::zero();
    let c = angle.to_radians().cos();
    let s = angle.to_radians().sin();
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
/// * `angle` - Angle of rotation in degrees
///
pub fn roty<T: Float>(angle: T) -> M33<T> {
    let one = T::one();
    let zero = T::zero();
    let c = angle.to_radians().cos();
    let s = angle.to_radians().sin();
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
/// * `angle` - Angle of rotation in degrees
///
pub fn rotz<T: Float>(angle: T) -> M33<T> {
    let one = T::one();
    let zero = T::zero();
    let c = angle.to_radians().cos();
    let s = angle.to_radians().sin();
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
/// phi: first euler angle in degrees (Float number)
/// theta: second euler angle in degrees (Float number)
/// psi: third euler angle in degrees (Float number)
/// s: Option<EulerSeq>: Optional Euler sequence if is None compute ZYZ
///
/// Output:
/// R: Rotation matrix(M33<Float>)
///
pub fn euler2rot<T: Float>(phi: T, theta: T, psi: T, s: Option<EulerSeq>) -> M33<T> {
    match s {
        Some(EulerSeq::XYX) => rotx(phi) * roty(theta) * rotx(psi),
        Some(EulerSeq::XYZ) => rotx(phi) * roty(theta) * rotz(psi),
        Some(EulerSeq::XZX) => rotx(phi) * rotz(theta) * rotx(psi),
        Some(EulerSeq::XZY) => rotx(phi) * rotz(theta) * roty(psi),
        Some(EulerSeq::YXY) => roty(phi) * rotx(theta) * roty(psi),
        Some(EulerSeq::YXZ) => roty(phi) * rotx(theta) * rotz(psi),
        Some(EulerSeq::YZX) => roty(phi) * rotz(theta) * rotx(psi),
        Some(EulerSeq::YZY) => roty(phi) * rotz(theta) * roty(psi),
        Some(EulerSeq::ZXY) => rotz(phi) * rotx(theta) * roty(psi),
        Some(EulerSeq::ZXZ) => rotz(phi) * rotx(theta) * rotz(psi),
        None                => rotz(phi) * roty(theta) * rotz(psi)
    }
}

// TODO(elsuizo:2020-09-20): do the tests
//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod test_transformations {

}
