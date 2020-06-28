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
