//-------------------------------------------------------------------------
// @file utils.rs
//
// @date 06/02/20 11:06:12
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
use num::{Float};

// NOTE(elsuizo:2020-06-02): the following function
// is a translation of the implementation that is here:
// https://floating-point-gui.de/errors/comparison/

/// a comparison function for floating point values
///
pub fn nearly_equal<T: Float>(a: T, b: T, epsilon: T) -> bool {
    let abs_a = a.abs();
    let abs_b = b.abs();
    let abs_diff = (a - b).abs();
    let zero = T::zero();
    // short-cut, handles infinity
    if a == b { true }

    else if a == zero || b == zero || (abs_a + abs_b < T::min_positive_value()) {
        // a or b is zero or both are extremely close to it
        // relative error is less meaningful here
        abs_diff < (epsilon * T::min_positive_value())
    } else {
      abs_diff / T::min(abs_a + abs_b, T::max_value()) < epsilon
    }
}

// TODO(elsuizo:2020-06-16): this is the only place when i use a Vec
/// utility function to compare vectors of Floats
pub fn compare_vecs<T: Float>(v1: &[T], v2: &[T], epsilon: T) -> bool {
    let v_result: Vec<bool> = v1
        .iter()
        .zip(v2)
        .map(|(a, b)| nearly_equal(*a, *b, epsilon))
        .collect();
    v_result.iter().all(|&x| x == true)
}


//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod test_utils {
    use super::*;
    #[test]
    fn test_nearly_equal() {
        let a = 0.15 + 0.15;
        let b = 0.1 + 0.2;
        assert_eq!(nearly_equal(a, b, 1e-9), true);
    }

    #[test]
    fn test_zero_signs() {
        let a = 0.0;
        let b = -0.0;
        assert_eq!(nearly_equal(a, b, 1e-9), true);
    }
}
