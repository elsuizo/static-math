//-------------------------------------------------------------------------
// @file dual_quaternions.rs
//
// @date 05/19/21 15:14:54
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
use static_math::{DualQuaternion, Quaternion, V3};
use static_math::transformations::{euler_to_rotation, homogeneous_from_quaternion};

// this example show some of the most common methods and function to work with `DualQuaternion`s
fn main() {
    let rot = euler_to_rotation(10f32.to_radians(), 10f32.to_radians(), 10f32.to_radians(), None);
    let q  = Quaternion::from_euler_angles(10f32.to_radians(), 10f32.to_radians(), 10f32.to_radians());
    // create a DualQuaternion and back again
    let t = homogeneous_from_quaternion(&q, &V3::new_from(1.0, 2.0, 3.0));
    let double = DualQuaternion::new_from_homogeneous(&t).to_homogeneous();

    let t_pure = DualQuaternion::new_from_translation(&V3::x_axis());
    let r_pure = DualQuaternion::new_from_rotation(&q);

    let normal = DualQuaternion::new_from_rot_trans(&q, &V3::x_axis());
    let combined = t_pure * r_pure;

    let r_pure2 = DualQuaternion::new_from_rotation_matrix(&rot);
    let combined2 = t_pure * r_pure2;

    println!("t: {}", t);
    println!("double: {}", double);

    println!("normal: {}", normal);
    println!("combined: {}", combined);
    println!("combined2: {}", combined2);
}


