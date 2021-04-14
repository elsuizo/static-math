//-------------------------------------------------------------------------
// @file lib.rs
//
// @date 06/18/20 11:26:55
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
//--------------------------------------------------------------------------
extern crate num;

pub mod traits;
pub mod matrix2x2;
pub mod matrix3x3;
pub mod matrix4x4;
pub mod matrix5x5;
pub mod matrix6x6;
pub mod vector2;
pub mod vector3;
pub mod vector4;
pub mod vector5;
pub mod vector6;
pub mod utils;
pub mod errors;
pub mod slices_methods;
pub mod quaternion;
pub mod transformations;
pub mod constants;
