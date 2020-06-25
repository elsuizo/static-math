//-------------------------------------------------------------------------
// @file inverse.rs
//
// @date 06/25/20 19:32:54
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
#[macro_use]
extern crate static_math;

use static_math::matrix6x6::M66;
use static_math::traits::LinearAlgebra;

fn main() {

    let m66 = m66_new!( 1.0,  1.0, 3.0,  4.0,  9.0, 3.0;
                       10.0, 10.0, 1.0,  2.0,  2.0, 5.0;
                        2.0,  9.0, 6.0, 10.0, 10.0, 9.0;
                       10.0,  9.0, 9.0,  7.0,  3.0, 6.0;
                        7.0,  6.0, 6.0,  2.0,  9.0, 5.0;
                        3.0,  8.0, 1.0,  4.0,  1.0, 5.0);

    if let Some(inv) = m66.inverse() {
        println!("inverse: {}", inv);
        println!("check: {}", inv * m66);
    }
}


