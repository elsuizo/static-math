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
use num::{Num, Float};
use crate::utils::compare_floats;

// TODO(elsuizo:2020-08-31): implement the Display trait for this type
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MaxMin<T> {
    max: (T, usize),
    min: (T, usize),
}

/// generic function to fin min and max values and the position in a slice
pub fn find_max_min<T: std::cmp::PartialOrd + Copy>(slice: &[T]) -> MaxMin<T> {
    let mut max = &slice[0];
    let mut min = &slice[0];

    let mut max_pos: usize = 0;
    let mut min_pos: usize = 0;

    for index in 1..slice.len() {
        if slice[index] < *min { min = &slice[index]; min_pos = index;}
        if slice[index] > *max { max = &slice[index]; max_pos = index;}
    }

    MaxMin{max: (*max, max_pos), min: (*min, min_pos)}
}

/// calculate the euclidean norm of the slice
pub fn norm2<T: Float>(slice: &[T]) -> T {
    slice.iter().fold(T::zero(), |n, &i| (i * i) + n).sqrt()
}

/// calculate the dot product of two slices
pub fn dot<T: Num + Copy + std::iter::Sum>(slice1: &[T], slice2: &[T]) -> T {
    slice1.iter().zip(slice2).map(|(&a, &b)| a * b).sum()
}

/// normalize the slice
pub fn normalize<T: Float>(slice: &mut [T]) {
    let n = norm2(slice);
    slice.iter_mut().for_each(|element| {
        *element = *element / n;
    })
}

/// project x in the direction of y
pub fn project_x_over_y<T: Float + std::iter::Sum>(x: &[T], y: &[T]) -> T {
    dot(x, y) / dot(y, y)
}

pub fn check_elements<T: Float>(v: &[T], tol: T) -> bool {
    let mut result = false;
    for num in v.iter() {
        result |= compare_floats(*num, T::zero(), tol);
    }
    result
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

        let expected = MaxMin{max: (37, 2), min: (1, 0)};

        assert_eq!(result, expected);

    }

    #[test]
    fn dot_tests() {
       let v1 = V3::new([1, 1, 1]);
       let v2 = V3::new([1, 1, 3]);

       let result = dot(&*v1, &*v2);
       let expected = 5;

       assert_eq!(result, expected);
    }

    #[test]
    fn normalize_test() {
        let mut v1 = V3::new([1.0, 1.0, 1.0]);
        normalize(&mut *v1);

        let expected = V3::new([0.5773502691896258, 0.5773502691896258, 0.5773502691896258]);

        assert_eq!(
            &v1[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &v1[..],
            &expected[..]
        );
    }
}
