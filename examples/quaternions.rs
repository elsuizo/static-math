//-------------------------------------------------------------------------
// @file quaternions.rs
//
// @date 09/11/20 11:48:54
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
extern crate static_math;

use static_math::quaternion::Quaternion;
use static_math::vector3::V3;

// In this example we rotate the x axis around the z axis 360 degrees
// to obtain the x axis again, but the rotation is via a composition of
// rotations of 90 degrees
fn main() {

    // vector to rotate: x axis: [1, 0, 0]
    let x = V3::x_axis();
    // quaternion represent the rotation around the z axis 90 degrees, the angle
    // is encoded in the vector norm: [0, 0, 90]
    let q = Quaternion::rotation_norm_encoded(V3::z_axis() * 90.0);
    let r = q * q * q * q * x;
    println!("r: {:}", r);
}
