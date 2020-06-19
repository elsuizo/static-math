//-------------------------------------------------------------------------
// @file slices_algorithms.rs
//
// @date 06/19/20 10:27:59
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
use crate::errors::VectorErrors;
use crate::vector2::V2;
use crate::vector3::V3;
use crate::vector4::V4;
use crate::vector5::V5;
use crate::vector6::V6;

enum VectorTypes {
    Vector2(V2),
    Vector3(V3),
    Vector4(V4),
    Vector5(V5),
    Vector6(V6)
}

impl VectorTypes

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
pub fn dot<T: Num + Copy>(slice1: &[T], slice2: &[T]) -> T {
    let mut result = T::zero();
    for i in 0..slice1.len() {
        result = result + slice1[i] * slice2[i];
    }
    result
}

pub fn norm2<T: Num + Copy + Float>(slice: &[T]) -> T {
    let mut result = T::zero();
    for element in 0..slice.len() {
        result = result + slice[element] * slice[element];
    }
    T::sqrt(result)
}

/// normalize a slice
pub fn normalize<T: Num + Copy + Float>(slice: &[T]) -> Result<(), VectorErrors> {
    let n = norm2(slice);
    let result = match slice.norm() {
        2 => V2::zeros(),
        3 => V3::zeros(),
        4 => V4::zeros(),
        5 => V5::zeros(),
        6 => V6::zeros()
    }

    if n != T::zero() {
        for i in 0..3 {
            slice[i] = slice[i] / n;
        }
        Ok(())
    } else {
        Err(VectorErrors::Norm2IsZero)
    }
}
