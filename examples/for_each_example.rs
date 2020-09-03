//-------------------------------------------------------------------------
// @file for_each_example.rs
//
// @date 09/03/20 15:08:39
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
#[macro_use]
extern crate static_math;

use static_math::matrix6x6::M66;

// Apply the closure f to all the elements of the matrix
fn main() {

    let m66 = m66_new!( 1.0,  1.0, 3.0,  4.0,  9.0, 3.0;
                       10.0, 10.0, 1.0,  2.0,  2.0, 5.0;
                        2.0,  9.0, 6.0, 10.0, 10.0, 9.0;
                       10.0,  9.0, 9.0,  7.0,  3.0, 6.0;
                        7.0,  6.0, 6.0,  2.0,  9.0, 5.0;
                        3.0,  8.0, 1.0,  4.0,  1.0, 5.0);

    let result = m66.for_each(|element: f32| element.cos() + element.sin());
    println!("result: {:}", result);
}
