//-------------------------------------------------------------------------
// @file eigenvalues.rs
//
// @date 08/11/20 12:00:53
// @author Martin Noblia
// @email mnoblia@disroot.org
//
// @brief
//
// @detail
//
// Licence MIT:
// Copyright <year> <Martin Noblia>
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

/// Example for calculate the real eigenvalues of a Matrix
/// This is a implementation of the amazing lecture of the profesor
/// Gilbert Strang: https://youtu.be/d32WV1rKoVk

#[macro_use]
extern crate static_math;

use static_math::matrix3x3::M33;
use static_math::traits::LinearAlgebra;
use static_math::slices_methods::check_elements;

fn convert_to_similar(m: &mut M33<f32>) {
    if let Some((q, r)) = m.qr() {
        *m = r * q;
    }
}

fn main() {

    let mut m = m33_new!(5.0, 2.0, 0.0;
                         2.0, 5.0, 0.0;
                         -3.0, 4.0, 6.0);

    let mut result = false;
    let mut counter = 0;
    let max_iterations = 200;

    while !result && counter < max_iterations {
        convert_to_similar(&mut m);
        let lower = m.get_lower_triangular();
        result = check_elements(&lower, 1e-8);
        counter += 1;
    }

    let eigenvalues = m.get_diagonal();
    println!("number of iterations: {:}", counter);
    println!("eigenvalues: {:}", eigenvalues);

}

