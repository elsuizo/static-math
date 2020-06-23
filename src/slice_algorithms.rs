//-------------------------------------------------------------------------
// @file slice_algorithms.rs
//
// @date 06/22/20 20:50:15
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
use num::{Float, Num};

/// Max and Min values in a slice
#[derive(Copy, Clone, Debug)]
pub struct MaxMin<T> {
    max: T,
    min: T,
}

/// find the max and min in slice
pub fn find_max_min<T: std::cmp::PartialOrd + Copy>(slice: &[T]) -> MaxMin<T> {
    let mut max = &slice[0];
    let mut min = &slice[0];

    for index in 1..slice.len() {
        if slice[index] < *min { min = &slice[index];}
        if slice[index] > *max { max = &slice[index];}
    }

    MaxMin{max: *max, min: *min}
}

/// dot product between two slices
fn dot<T: Num + Copy + std::iter::Sum>(slice1: &[T], slice2: &[T]) -> T {
    slice1.iter().zip(slice2).map(|(&a, &b)| a * b).sum()
}

pub fn norm2<T: Num + Copy + Float>(slice: &[T]) -> T {
    slice.iter().fold(T::one(), |n, &i| (i * i) + n).sqrt()
}

// TODO(elsuizo:2020-06-22): faltaria esta solamente
/// normalize a slice
pub fn normalize<T: Num + Copy + Float>(slice: &mut [T]) -> Result<(), VectorErrors> {
    let n = norm2(slice);
    if n != T::zero() {
        //slice.iter_mut().map(|element| element / n);
        for i in 0..slice.len() {
            slice[i] = slice[i] / n;
        }
        Ok(())
    } else {
        Err(VectorErrors::Norm2IsZero)
    }
}

/// orthogonal projection of slice `x` over the slide `y`
fn projection_x_over_y<T: Num + Copy + std::iter::Sum>(x: &[T], y: &[T]) -> T {
    dot(x, y) / dot(y, y)
}
