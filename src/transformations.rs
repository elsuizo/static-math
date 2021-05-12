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
use crate::matrix6x6::M66;
use crate::m44_new;
use crate::vector3::V3;
use crate::vector4::V4;
use crate::vector6::V6;
use crate::quaternion::Quaternion;
use crate::traits::LinearAlgebra;
use crate::utils::{nearly_zero, nearly_equal};
use num::{Float, Signed};
use num::traits::FloatConst;
use crate::slices_methods::norm2;

/// Euler sequences conventions of rotations
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

/// Compute the rotation around the `x` axis(in cartesian coordinates)
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
/// Get the euler angles from a rotation matrix comming from the convention ZYX
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
/// Generate a homogeneous matrix from a rotation represented by a Matrix and a translation
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

pub fn se3_from_parts<T: Float>(r: &M33<T>, p: &V3<T>) -> M44<T> {
    let zero = T::zero();
    m44_new!(r[(0, 0)], r[(0, 1)], r[(0, 2)], p[0];
             r[(1, 0)], r[(1, 1)], r[(1, 2)], p[1];
             r[(2, 0)], r[(2, 1)], r[(2, 2)], p[2];
                zero  ,   zero   ,    zero  , zero)
}

pub fn se3_get_parts<T: Float>(se3: &M44<T>) -> (M33<T>, V3<T>) {
    let skew_omega = m33_new!(se3[(0, 0)], se3[(0, 1)], se3[(0, 2)];
                              se3[(1, 0)], se3[(1, 1)], se3[(1, 2)];
                              se3[(2, 0)], se3[(2, 1)], se3[(2, 2)]);
    let v = V3::new_from(se3[(0, 3)], se3[(1, 3)], se3[(2, 3)]);
    (skew_omega, v)
}

/// Generate a homogeneous matrix from a rotation represented by a quaternion and a translation
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

/// Get the parts of a homogeneous transformation, the rotation(expresed by a Quaternion) and the
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

/// Get the parts of a homogeneous transformation, the rotation(expresed by a Matrix) and the
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

/// Get the inverse of the homogeneous transformation
///
/// Function arguments:
/// r: Homogeneus transformation reference (&M44<Float>)
///
pub fn homogeneous_inverse<T: Float + std::iter::Sum + Signed>(r: &M44<T>) -> M44<T> {
    let (rot, p) = get_parts_raw(r);
    let rot_new = rot.transpose();
    let p_new   = -rot_new * p;
    homogeneous_from_rotation(&rot_new, &p_new)
}

/// Transform a vector with the inverse of a given homogeneous transformation
///
/// Function arguments:
/// r: Homogeneus transformation reference (&M44<Float>)
///
/// p: reference to vector to trasnform(&V3<T>)
///
pub fn homogeneous_inverse_transform<T>(r: &M44<T>, p: &V3<T>) -> V3<T>
where
    T: Float + std::iter::Sum + FloatConst + Signed
{
    let (_, translation) = get_parts(r);
    let p_trans = *p - translation;
    let p_new = V4::new_from(p_trans[0], p_trans[1], p_trans[2], T::zero());
    let result = homogeneous_inverse(r) * p_new;
    V3::new_from(result[0], result[1], result[2])
}

/// Convert a 3d Vector to a so(3) representation
///
/// Function arguments:
/// `v`: V3<Float>
///
/// Output: A M33<Float> that is the skew symetric representation of `v`
///
pub fn skew_from_vec<T: Float>(v: &V3<T>) -> M33<T> {
    let zero = T::zero();
    m33_new!( zero, -v[2],  v[1];
              v[2],  zero, -v[0];
             -v[1],  v[0],  zero)
}

/// Convert an so(3) representation to a 3d vector
///
/// Function arguments:
/// `r`: M33<Float>
///
/// Output: A V3<Float>
///
pub fn skew_to_vec<T: Float>(r: &M33<T>) -> V3<T> {
    V3::new_from(r[(2, 1)], r[(0, 2)], r[(1, 0)])
}

pub fn skew_scalar<T: Float>(number: T) -> M22<T> {
    let zero = T::zero();
    m22_new!(  zero, -number;
             number,    zero)
}

pub fn vex_m22<T: Float>(m: &M22<T>) -> T {
    m[(1, 0)]
}

/// Create a pose in 2D from a angle(in radians) and a cartesian position (x, y) values
pub fn ksi<T: Float>(angle: T, x: T, y: T) -> M33<T> {
    let zero = T::zero();
    let one = T::one();
    let (s, c) = angle.sin_cos();
    m33_new!(   c,   -s,   x;
                s,    c,   y;
             zero, zero, one)
}

/// Create augmented skew-symmetric matrix
pub fn skew_from_vec_aug<T: Float>(v: V3<T>) -> M33<T> {
    let zero = T::zero();
    m33_new!(zero, -v[2], v[0];
             v[2],  zero, v[1];
             zero,  zero, zero)
}

/// Create augmented skew-symmetric matrix
pub fn skew_v6<T: Float>(v: V6<T>) -> M44<T> {
    let zero = T::zero();
    m44_new!( zero, -v[5],  v[4], v[0];
              v[5],  zero, -v[3], v[1];
             -v[4],  v[3],  zero, v[2];
              zero,  zero,  zero, zero)
}

/// Convert a 3d vector of exponential coordinates for rotation into axis-angle
/// form
///
/// Function arguments:
/// `exp`: V3<Float> 3d vector of exponential coordinates
///
/// Output:
/// (omega_hat, theta)
///
pub fn axis_angle<T: Float>(exp: &V3<T>) -> (V3<T>, T) {
    (exp.normalize().expect("empty vector"), exp.norm2())
}

/// Create a pure translation homogeneous transformation
pub fn translation_2d<T: Float>(x: T, y: T) -> M33<T> {
    let zero = T::zero();
    let one = T::one();
    m33_new!( one, zero,   x;
             zero,  one,   y;
             zero, zero, one)
}

/// Create a pure translation homogeneous transformation in 3d
pub fn translation_3d<T: Float>(x: T, y: T, z: T) -> M44<T> {
    let zero = T::zero();
    let one = T::one();
    m44_new!( one, zero, zero,  x;
             zero,  one, zero,  y;
             zero, zero,  one,  z;
             zero, zero, zero, one)
}

// NOTE(elsuizo:2021-05-04): wide code is better code
/// Compute the matrix exponential for omega theta(exponential coordinates): so(3) ---> SO(3)
/// with the Rodriguez formula
pub fn rotation_from_axis_angle<T: Float>(omega: &V3<T>, theta: T) -> M33<T> {
    let skew       = skew_from_vec(omega);
    let (sin, cos) = theta.sin_cos();
    M33::identity() + skew * sin + skew * skew * (T::one() - cos)
}

// NOTE(elsuizo:2021-05-09): remember we have problems to Mul to left side with a constant...
/// Compute the matrix exponential of a matrix in so(3)
///
/// Function arguments:
/// `so3`: M33<T> a screw symetric matrix
///
/// Output:
/// M33<T>: representing the matrix exponential of the input
///
pub fn matrix_exponential<T: Float>(so3: &M33<T>) -> M33<T> {
    let omega_theta = skew_to_vec(so3);
    let one = T::one();
    if nearly_zero(omega_theta.norm2()) {
        return M33::identity();
    } else {
        let theta = axis_angle(&omega_theta).1;
        let omega = *so3 / theta;
        let (s, c) = theta.sin_cos();
        return M33::identity() + omega * s + omega * omega * (one - c);
    }
}

/// Brief.
///
/// Compute the matrix logarithm of a rotation matrix
///
/// Function arguments:
/// `r`: &M33<Float> rotation matrix
///
/// Output:
/// M33<Float>: representing the matrix logarithm of the input
///
pub fn matrix_log<T: Float + std::iter::Sum + FloatConst>(r: &M33<T>) -> M33<T> {
    let one = T::one();
    let two = one + one;
    let angle = (r.trace() - one) / two;
    if angle >= one {
        return M33::zeros();
    } else if angle <= -one {
        if !nearly_zero(one + r[(2, 2)]) {
            let omega = V3::new_from(r[(0, 2)], r[(1, 2)], one + r[(2, 2)]) * (one / T::sqrt(two * (one + r[(2, 2)])));
            return skew_from_vec(&(omega * T::PI()))
        } else if !nearly_zero(one + r[(1, 1)]) {
            let omega = V3::new_from(r[(0, 1)], one + r[(1, 1)], r[(2, 1)]) * (one / T::sqrt(two * (one + r[(1, 1)])));
            return skew_from_vec(&(omega * T::PI()))
        } else {
            let omega = V3::new_from(one + r[(0, 0)], one + r[(1, 0)], r[(2, 0)]) * (one / T::sqrt(two * (one + r[(0, 0)])));
            return skew_from_vec(&(omega * T::PI()))
        }
    } else {
        let theta = T::acos(angle);
        // TODO(elsuizo:2021-05-10): esto esta bien???
        return (*r - r.transpose()) * (theta / two / T::sin(theta));
    }
}

/// Inverse of a Rotation matrix, where R: SO(3)
pub fn rotation_inverse<T: Float + std::iter::Sum>(r: &M33<T>) -> M33<T> {
    r.transpose()
}

/// Convert a spatial velocity vector to a M44 matrix in se3 space
///
/// Function arguments:
/// `twist`: V6<Float> representing the "twist"
///
/// Output:
/// M44<Float>: a matrix in the se3 space
///
pub fn twist_to_se3<T: Float>(twist: &V6<T>) -> M44<T> {
    let skew = skew_from_vec(&V3::new_from(twist[0], twist[1], twist[2]));
    se3_from_parts(&skew, &V3::new_from(twist[3], twist[4], twist[5]))
}

/// Convert a se3 matrix representation to a spatial velocity vector known as "twist"
///
/// Function arguments:
/// `se3`: M44<Float> a matrix in the se3 space
///
/// Output:
/// V3<Float>: representing spatial velocity vector("twist")
///
pub fn se3_to_twist<T: Float>(se3: &M44<T>) -> V6<T> {
    V6::new_from(se3[(2, 1)], se3[(0, 2)], se3[(1, 0)], se3[(0, 3)], se3[(1, 3)], se3[(2, 3)])
}

/// Compute the adjoint representation of a homogeneous transformation matrix
///
/// Function arguments:
/// `transform`: M44<Float> a homogeneous transformation matrix
///
/// Output:
/// M66<Float>: the 6x6 matrix representation of the Adjoint of T
///
pub fn adjoint<T: Float>(transform: &M44<T>) -> M66<T> {
    let mut result = M66::zeros();
    let (r, p) = get_parts_raw(transform);
    let skew = skew_from_vec(&p);
    let transform_p = skew * r;
    result.copy_elements_from(&r, 1);
    result.copy_elements_from(&r, 4);
    result.copy_elements_from(&transform_p, 3);
    result
}

pub fn screw_to_axis<T: Float>(q: &V3<T>, s: &V3<T>, h: T) -> V6<T> {
    let v = q.cross(*s) + (*s) * h;
    V6::new_from(s[0], s[1], s[2], v[0], v[1], v[2])
}

/// Convert a 6D vector of exponential coordinates into screw axis angle
///
/// Function arguments:
/// `exp`: V6<Float> A 6D vector of exponential coordinates for rigid-body motion S * theta
///
/// Output:
/// (V6<Float, Float) a tuple with the first element the normalized "screw" axis, and the second
/// tuple element called "theta" representing the distance traveled along/about S
///
pub fn axis_angle6<T: Float>(exp: &V6<T>) -> (V6<T>, T) {
    // NOTE(elsuizo:2021-05-11): here we treat the vector like a slice(to obtain the first three
    // elements)
    let theta = norm2(&exp[0..3]);
    if nearly_zero(theta) {
        let theta = norm2(&exp[3..6]);
        return (*exp / theta, theta)
    }
    return (*exp / theta, theta)
}

pub fn matrix_exponential6<T: Float>(se3: &M44<T>) -> M44<T> {
    let zero = T::zero();
    let one  = T::one();
    let (skew_omega, v) = se3_get_parts(se3);
    let omega_theta = skew_to_vec(&skew_omega);
    // NOTE(elsuizo:2021-05-11): if w = 0 ---> "pitch" == infinity and |v| = 1 ---> theta represents only a linear distance
    if nearly_zero(omega_theta.norm2()) {
        m44_new!( one,  zero,  one, v[0];
                 zero,   one, zero, v[1];
                 zero,  zero,  one, v[2];
                 zero,  zero, zero, one)
    } else {
        // NOTE(elsuizo:2021-05-11): |w| = 1 ---> theta represents the distance along side the
        // "screw" axis S
        let theta = axis_angle(&omega_theta).1;
        let omega = skew_omega / theta;
        let mat_exp = matrix_exponential(&skew_omega);
        let (s, c) = theta.sin_cos();
        let v_t = (M33::identity() * theta + omega * (one - c) + omega * omega * (theta - s)) * v / theta;

        m44_new!(mat_exp[(0, 0)], mat_exp[(0, 1)], mat_exp[(0, 2)], v_t[0];
                 mat_exp[(1, 0)], mat_exp[(1, 1)], mat_exp[(1, 2)], v_t[1];
                 mat_exp[(2, 0)], mat_exp[(2, 1)], mat_exp[(2, 2)], v_t[2];
                        zero    ,     zero       ,        zero    ,  one)
    }
}
//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod test_transformations {

    use super::{rotation_to_euler, euler_to_rotation, rotx, roty, rotz,
                homogeneous_from_quaternion, homogeneous_inverse,
                homogeneous_inverse_transform, matrix_log, matrix_exponential, skew_from_vec,
                twist_to_se3, se3_to_twist, adjoint, screw_to_axis, axis_angle6, matrix_exponential6};
    use crate::utils::{nearly_equal, is_rotation, is_rotation_h, compare_vecs};
    use crate::quaternion::Quaternion;
    use crate::vector3::V3;
    use crate::matrix4x4::M44;
    use crate::m44_new;
    use crate::matrix6x6::M66;
    use crate::m66_new;
    use crate::vector6::V6;
    use crate::matrix3x3::M33;
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

    #[test]
    fn matrix_log_test() {
        let r = m33_new!(0.0, 0.0, 1.0;
                         1.0, 0.0, 0.0;
                         0.0, 1.0, 0.0);
        let result = matrix_log(&r);
        let expected = m33_new!(        0.0, -1.20919958,  1.20919958;
                                 1.20919958,         0.0, -1.20919958;
                                -1.20919958,  1.20919958,         0.0);
        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn matrix_exponential_test() {
        let so3 = skew_from_vec(&V3::new_from(1.0, 2.0, 3.0));
        let result = matrix_exponential(&so3);
        let expected  = m33_new!(-0.69492056,  0.71352099,  0.08929286;
                                 -0.19200697, -0.30378504,  0.93319235;
                                  0.69297817,  0.6313497 ,  0.34810748);
        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn twist_to_se3_test() {
        let result = twist_to_se3(&V6::new_from(1.0, 2.0, 3.0, 4.0, 5.0, 6.0));
        let expected = m44_new!( 0.0, -3.0,  2.0, 4.0;
                                 3.0,  0.0, -1.0, 5.0;
                                -2.0,  1.0,  0.0, 6.0;
                                 0.0,  0.0,  0.0, 0.0);
        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn se3_to_twist_test() {
        let se3 = m44_new!( 0.0, -3.0,  2.0, 4.0;
                            3.0,  0.0, -1.0, 5.0;
                           -2.0,  1.0,  0.0, 6.0;
                            0.0,  0.0,  0.0, 0.0);
        let result = se3_to_twist(&se3);
        let expected = V6::new_from(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn adjoint_test() {
        let t = m44_new!(1.0, 0.0,  0.0, 0.0;
                         0.0, 0.0, -1.0, 0.0;
                         0.0, 1.0,  0.0, 3.0;
                         0.0, 0.0,  0.0, 1.0);

        let result = adjoint(&t);
        let expected = m66_new!(1.0, 0.0,  0.0, 0.0, 0.0,  0.0;
                                0.0, 0.0, -1.0, 0.0, 0.0,  0.0;
                                0.0, 1.0,  0.0, 0.0, 0.0,  0.0;
                                0.0, 0.0,  3.0, 1.0, 0.0,  0.0;
                                3.0, 0.0,  0.0, 0.0, 0.0, -1.0;
                                0.0, 0.0,  0.0, 0.0, 1.0,  0.0);

        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }

    #[test]
    fn screw_to_axis_test() {
        let q = V3::new_from(3.0, 0.0, 0.0);
        let s = V3::new_from(0.0, 0.0, 1.0);
        let h = 2.0;
        let result = screw_to_axis(&q, &s, h);
        let expected = V6::new_from(0.0, 0.0, 1.0, 0.0, -3.0, 2.0);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn axis_angle6_test() {
        let v = V6::new_from(1.0, 0.0, 0.0, 1.0, 2.0, 3.0);
        let result = axis_angle6(&v);
        let expected = (V6::new_from(1.0, 0.0, 0.0, 1.0, 2.0, 3.0), 1.0);
        assert_eq!(
            &result.0[..],
            &expected.0[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result.0[..],
            &expected.0[..]
        );
        assert!(nearly_equal(result.1, expected.1, EPS));
    }

    #[test]
    fn matrix_exponential6_test() {
        let se3mat = m44_new!(0.0,          0.0,           0.0,          0.0;
                              0.0,          0.0,   -1.57079632,   2.35619449;
                              0.0,   1.57079632,           0.0,   2.35619449;
                              0.0,          0.0,           0.0,          0.0);

        let result = matrix_exponential6(&se3mat);
        let expected = m44_new!(1.0, 0.0,  0.0, 0.0;
                                0.0, 0.0, -1.0, 0.0;
                                0.0, 1.0,  0.0, 3.0;
                                0.0, 0.0,  0.0, 1.0);
        assert!(compare_vecs(&result.as_vec(), &expected.as_vec(), EPS));
    }
}
