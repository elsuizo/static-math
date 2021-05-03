//-------------------------------------------------------------------------
// @file qr_example.rs
//
// @date 06/25/20 19:30:44
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
use static_math::{M22, M33, M44, M55, M66, m22_new, m33_new, m44_new, m55_new, m66_new};
use static_math::traits::LinearAlgebra;

fn main() {

    let m22 = m22_new!(1.0, 2.0;
                       3.0, 4.0);

    let m33 = m33_new!(10.0, -1.0, 2.0;
                        3.0,  3.0, 1.0;
                        1.0,  2.0, 5.0);

    let m44 = m44_new!( 0.0,  1.0,  2.0,  3.0;
                        4.0,  0.0,  1.0,  7.0;
                       10.0,  9.0, 10.0, 11.0;
                       12.0, 13.0, 14.0, 15.0);

    let m55 = m55_new!(10.0, 1.0, 7.0,  1.0,  5.0;
                        2.0, 4.0, 8.0,  3.0,  2.0;
                        5.0, 1.0, 2.0,  9.0, 10.0;
                        6.0, 9.0, 1.0,  7.0,  3.0;
                        1.0, 8.0, 8.0, 10.0,  5.0);

    let m66 = m66_new!( 1.0,  1.0, 3.0,  4.0,  9.0, 3.0;
                       10.0, 10.0, 1.0,  2.0,  2.0, 5.0;
                        2.0,  9.0, 6.0, 10.0, 10.0, 9.0;
                       10.0,  9.0, 9.0,  7.0,  3.0, 6.0;
                        7.0,  6.0, 6.0,  2.0,  9.0, 5.0;
                        3.0,  8.0, 1.0,  4.0,  1.0, 5.0);

    if let Some((q, r)) = m22.qr() {
        println!("q: {}", q);
        println!("r: {}", r);
        println!("m22: {}", q * r);
        println!("q.det(): {}", q.det());
    }

    if let Some((q, r)) = m33.qr() {
        println!("q: {}", q);
        println!("r: {}", r);
        println!("m33: {}", q * r);
        println!("q.det(): {}", q.det());
    }

    if let Some((q, r)) = m44.qr() {
        println!("q: {}", q);
        println!("r: {}", r);
        println!("m44: {}", q * r);
        println!("q.det(): {}", q.det());
    }

    if let Some((q, r)) = m55.qr() {
        println!("q: {}", q);
        println!("r: {}", r);
        println!("m55: {}", q * r);
        println!("q.det(): {}", q.det());
    }

    if let Some((q, r)) = m66.qr() {
        println!("q: {}", q);
        println!("r: {}", r);
        println!("m66: {}", q * r);
        println!("q.det(): {}", q.det());
    }

}
