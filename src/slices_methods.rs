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

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MaxMin<T> {
    max: T,
    min: T,
}

/// generic function to fin min and max values in a slice
pub fn find_max_min<T: std::cmp::PartialOrd + Copy>(slice: &[T]) -> MaxMin<T> {
    let mut max = &slice[0];
    let mut min = &slice[0];

    for index in 1..slice.len() {
        if slice[index] < *min { min = &slice[index];}
        if slice[index] > *max { max = &slice[index];}
    }

    MaxMin{max: *max, min: *min}
}

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

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod test_slides_methods {

    use crate::vector3::V3;
    use crate::slices_methods::*;

    #[test]
    fn find_max_min_test() {
        let v = V3::new([1, 10, 37]);

        let result = find_max_min(&*v);

        let expected = MaxMin{max: 37, min: 1};

        assert_eq!(result, expected);

    }
}
