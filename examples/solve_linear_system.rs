//-------------------------------------------------------------------------
// @file solve_linear_system.rs
//
// @date 06/27/20 22:43:12
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
use static_math::{M66, V6, m66_new};
use static_math::traits::LinearAlgebra;

fn main() {

    let a = m66_new!( 1.0,  1.0, 3.0,  4.0,  9.0, 3.0;
                     10.0, 10.0, 1.0,  2.0,  2.0, 5.0;
                      2.0,  9.0, 6.0, 10.0, 10.0, 9.0;
                     10.0,  9.0, 9.0,  7.0,  3.0, 6.0;
                      7.0,  6.0, 6.0,  2.0,  9.0, 5.0;
                      3.0,  8.0, 1.0,  4.0,  1.0, 5.0);

    let b = V6::new_from(0.0, 1.0, 3.0, 0.0, 1.0, 2.0);

    if let Some(inv) = a.inverse() {
        let solution = inv * b;
        println!("the solution is: {}", solution);
        println!("verification: a * solution = b?: {}", a * solution);
        println!("b: {}", b);
    }
}
