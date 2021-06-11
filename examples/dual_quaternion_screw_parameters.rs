//-------------------------------------------------------------------------
// @file dual_quaternion_screw_parameters.rs
//
// @date 05/22/21 20:49:24
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
use static_math::{Quaternion, DualQuaternion, V3};

// In this example we rotate around the z axis 45 degrees and translate along the z axis
// 7 units, so the l vector should be z_hat, d should be 7, theta should be 45
// degrees and m should be zero because the moment is zero
fn main() {
    let v = V3::z_axis() * 45f32.to_radians();
    let q = Quaternion::rotation_norm_encoded(&v);
    let trans = V3::new_from(0.0, 0.0, 7.0);
    let dq_full = DualQuaternion::new_from_rot_trans(&q, &trans);
    let (l, m, theta, d) = dq_full.get_screw_parameters();
    println!("l: {}", l);
    println!("m: {}", m);
    println!("theta: {}", theta.to_degrees());
    println!("d: {}", d);
}


