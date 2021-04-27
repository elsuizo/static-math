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
use crate::vector4::V4;
use crate::quaternion::Quaternion;
use crate::traits::LinearAlgebra;
use crate::utils::{nearly_equal};
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

// TODO(elsuizo:2021-04-27): handle only valid rotations
// TODO(elsuizo:2021-04-23): handle all the rotation sequences
/// Brief.
///
/// get the euler angles from a rotation matrix comming from the convention ZYX
///
/// Function arguments:
///
/// `r`: a reference to a Rotation matrix(M33<Float>)
///
/// Output:
///
/// Euler angles: (yay, pitch, roll)
///
pub fn rotation_to_euler<T: Float + FloatConst>(r: &M33<T>) -> (T, T, T) {
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
/// Brief.
///
/// generate a homogeneous matrix from a rotation represented by a Matrix and a translation
/// represented by a vector
///
/// Function arguments:
///
/// q: a reference to a M33<Float> (Rotation part)
///
/// p: a reference to a  V3<Float> (Translation part)
///
///
pub fn homogeneous_from_rotation<T: Float>(r: &M33<T>, p: &V3<T>) -> M44<T> {
    let zero = T::zero();
    let one  = T::one();
    m44_new!(r[(0, 0)], r[(0, 1)], r[(0, 2)], p[0];
             r[(1, 0)], r[(1, 1)], r[(1, 2)], p[1];
             r[(2, 0)], r[(2, 1)], r[(2, 2)], p[2];
                zero  ,   zero   ,    zero  , one)
}

/// Brief.
///
/// generate a homogeneous matrix from a rotation represented by a quaternion and a translation
/// represented by a vector
///
/// Function arguments:
///
/// q: a reference to a Quaternion<Float> (Rotation part)
///
/// p: a reference to a V3<Float> (Translation part)
///
///
pub fn homogeneous_from_quaternion<T: Float>(q: &Quaternion<T>, p: &V3<T>) -> M44<T> {
    homogeneous_from_rotation(&q.to_rotation(), p)
}

/// Brief.
///
/// get the parts of a homogeneous transformation, the rotation(expresed by a Quaternion) and the
/// translation (expresed by a vector)
///
/// Function arguments:
/// r: Homogeneus transformation reference (&M44<Float>)
///
pub fn get_parts<T: Float + FloatConst>(r: &M44<T>) -> (Quaternion<T>, V3<T>) {
    let rot = m33_new!(r[(0, 0)], r[(0, 1)], r[(0, 2)];
                       r[(1, 0)], r[(1, 1)], r[(1, 2)];
                       r[(2, 0)], r[(2, 1)], r[(2, 2)]);

    let p = V3::new_from(r[(0, 3)], r[(1, 3)], r[(2, 3)]);

    let (yay, pitch, roll) = rotation_to_euler(&rot);
    let q = Quaternion::from_euler_angles(yay, pitch, roll);

    (q, p)
}

/// Brief.
///
/// get the parts of a homogeneous transformation, the rotation(expresed by a Matrix) and the
/// translation (expresed by a vector)
///
/// Function arguments:
/// r: Homogeneus transformation reference(&M44<Float>)
///
pub fn get_parts_raw<T: Float>(r: &M44<T>) -> (M33<T>, V3<T>) {
    let rot = m33_new!(r[(0, 0)], r[(0, 1)], r[(0, 2)];
                       r[(1, 0)], r[(1, 1)], r[(1, 2)];
                       r[(2, 0)], r[(2, 1)], r[(2, 2)]);

    let p = V3::new_from(r[(0, 3)], r[(1, 3)], r[(2, 3)]);

    (rot, p)
}

/// Brief.
///
/// get the inverse of the homogeneous transformation
///
/// Function arguments:
/// r: Homogeneus transformation reference (&M44<Float>)
///
pub fn homogeneous_inverse<T: Float + std::iter::Sum>(r: &M44<T>) -> M44<T> {
    let one = T::one();
    let (rot, p) = get_parts_raw(r);
    let rot_new = rot.transpose();
    // TODO(elsuizo:2021-04-27): this is because M33<T> and all the matrix dont impl Neg
    let p_new   = rot_new * -one * p;
    homogeneous_from_rotation(&rot_new, &p_new)
}

/// Brief.
///
/// transform a vector with the inverse of a given homogeneous transformation
///
/// Function arguments:
/// r: Homogeneus transformation reference (&M44<Float>)
///
/// p: reference to vector to trasnform(&V3<T>)
///
pub fn homogeneous_inverse_transform<T>(r: &M44<T>, p: &V3<T>) -> V3<T>
where
    T: Float + std::iter::Sum + FloatConst
{
    let (_, translation) = get_parts(r);
    let p_trans = *p - translation;
    let p_new = V4::new_from(p_trans[0], p_trans[1], p_trans[2], T::zero());
    let result = homogeneous_inverse(r) * p_new;
    V3::new_from(result[0], result[1], result[2])
}

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod test_transformations {

    use super::{rotation_to_euler, euler_to_rotation, rotx, roty, rotz, homogeneous_from_quaternion, homogeneous_inverse, homogeneous_inverse_transform};
    use crate::utils::{nearly_equal, is_rotation, is_rotation_h, compare_vecs};
    use crate::quaternion::Quaternion;
    use crate::vector3::V3;
    use crate::matrix4x4::M44;
    const EPS: f32 = 1e-6;

    #[test]
    fn rotation_and_euler_test() {
        let expected = (0.1, 0.2, 0.3);
        let r = euler_to_rotation(expected.0, expected.1, expected.2, None);
        let result = rotation_to_euler(&r);
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

    #[test]
    fn homogeneous_transform_test() {
        let q = Quaternion::rotation_norm_encoded(V3::y_axis() * std::f32::consts::FRAC_PI_2);
        let iso_s = homogeneous_from_quaternion(&q, &V3::new_from(0.0, 0.0, 3.0));
        assert!(is_rotation_h(iso_s));
    }

    #[test]
    fn homogeneous_inverse_test() {
        let q = Quaternion::rotation_norm_encoded(V3::y_axis() * std::f32::consts::FRAC_PI_2);
        let iso_s = homogeneous_from_quaternion(&q, &V3::new_from(0.0, 0.0, 3.0));
        let result = iso_s * homogeneous_inverse(&iso_s);
        let expected = M44::identity();
        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    // NOTE(elsuizo:2021-04-27): this is the same example from nalgebra:
    // https://docs.rs/nalgebra/0.26.2/nalgebra/geometry/struct.Isometry.html#method.inverse_transform_point
    #[test]
    fn inverse_point_transform_test() {
        let q = Quaternion::rotation_norm_encoded(V3::y_axis() * std::f32::consts::FRAC_PI_2);
        let iso_s = homogeneous_from_quaternion(&q, &V3::new_from(0.0, 0.0, 3.0));
        let result = homogeneous_inverse_transform(&iso_s, &V3::new_from(1.0, 2.0, 3.0));
        let expected = V3::new_from(0.0, 2.0, 1.0);
        assert!(compare_vecs(&*result, &*expected, EPS));
    }
}
