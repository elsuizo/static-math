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
//  Licence:
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
//--------------------------------------------------------------------------
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


