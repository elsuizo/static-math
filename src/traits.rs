//-------------------------------------------------------------------------
// @file traits.rs
//
// @date 06/01/20 22:19:00
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
/// Generic Trait for Matrix operations and Linear Algebra methods
pub trait LinearAlgebra<T> {
    /// get the rows of the matrix
    fn rows(&self) -> usize;

    /// get the columns of the matrix
    fn cols(&self) -> usize;

    /// get the overal shape of the matrix
    fn shape(&self) -> (usize, usize) {
        (self.rows(), self.cols())
    }

    /// transpose dimentions of the matrix
    fn transpose(&self) -> Self;

    /// get the trace of the matrix
    fn trace(&self) -> T;

    // FIXME(elsuizo:2020-06-01): me parece que no se puede hacer generica porque
    // su implementacion lleva sqrt
    /// compute the euclidean norm of the matrix
    fn norm2(&self) -> T;

    /// compute the determinant of the matrix
    fn det(&self) -> T;

    /// compute the inverse of the matrix
    fn inverse(&self) -> Option<Self>
    where
        Self: Sized;
}
