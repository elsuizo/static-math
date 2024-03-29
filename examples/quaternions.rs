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
use static_math::{V3, Quaternion};
use static_math::vector3::X_AXIS;
// In this example we rotate the x axis around the z axis 360 degrees
// to obtain the x axis again, but the rotation is via a composition of
// rotations of 90 degrees
fn main() {

    // vector to rotate: x axis: [1, 0, 0]
    let x = X_AXIS;
    // quaternion represent the rotation around the z axis 90 degrees, the angle
    // is encoded in the vector norm: [0, 0, 90]
    let v = V3::z_axis() * 90f32.to_radians();
    let q = Quaternion::rotation_norm_encoded(&v);
    let r = q * q * q * q * x;
    println!("r: {:}", r);
    //-------------------------------------------------------------------------
    //                        Quaternions and euler angles
    //-------------------------------------------------------------------------
    let q = Quaternion::from_euler_angles(0.1, 0.2, 0.3);
    println!("q: {}", q);
    let euler_angles = q.to_euler_angles();
    // this would have to give the same value :)
    println!("euler_angles: {:?}", euler_angles);
}
