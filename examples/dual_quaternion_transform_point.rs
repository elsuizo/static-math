//-------------------------------------------------------------------------
// @file dual_quaternion_transform_point.rs
//
// @date 05/20/21 19:39:58
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
//-------------------------------------------------------------------------
use static_math::{DualQuaternion, Quaternion, V3, V4};
use static_math::transformations::{homogeneous_from_quaternion};

// In this example we transform a point with a homogeneous matrix in SE(3) and with
// a DualQuaternion to obtain the same result
// This representation is more compact since there are only 8 elements instead of the 16 of the
// homogeneous matrix
fn main() {

    let q  = Quaternion::from_euler_angles(10f32.to_radians(), 30f32.to_radians(), 45f32.to_radians());
    let t = homogeneous_from_quaternion(&q, &V3::new_from(1.0, 2.0, 3.0));

    let p = V4::new_from(1.0, 2.0, 3.0, 0.0);
    let expected = t * p;

    let p = V3::new_from(1.0, 2.0, 3.0);
    let dq = DualQuaternion::new_from_homogeneous(&t);
    let result = dq.transform_point(&p);

    println!("expected: {}", expected);
    println!("result: {}", result);
}

