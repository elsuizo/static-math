//-------------------------------------------------------------------------
// @file slices_methods.rs
//
// @date 06/24/20 21:48:21
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
use num::{Num, Float};

/// calculate the euclidean norm of the slice
pub fn norm2<T: Num + Copy + Float>(slice: &[T]) -> T {
    slice.iter().fold(T::zero(), |n, &i| (i * i) + n).sqrt()
}

/// calculate the dot product of two slices
pub fn dot<T: Num + Copy + std::iter::Sum>(slice1: &[T], slice2: &[T]) -> T {
    slice1.iter().zip(slice2).map(|(&a, &b)| a * b).sum()
}

/// normalize the slice
pub fn normalize<T: Float>(slice: &mut [T]) {
    let n = norm2(slice);
    for element in slice.iter_mut() {
        *element = *element / n;
    }
}

/// project x in the direction of y
pub fn project_x_over_y<T: Float + std::iter::Sum>(x: &[T], y: &[T]) -> T {
    dot(x, y) / dot(y, y)
}
