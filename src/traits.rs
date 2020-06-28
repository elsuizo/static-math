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

    /// compute the euclidean norm of the matrix
    fn norm2(&self) -> T;

    /// compute the determinant of the matrix
    fn det(&self) -> T;

    /// compute the inverse of the matrix
    fn inverse(&self) -> Option<Self>
    where
        Self: Sized;

    /// compute the QR factorization of the matrix(if has inverse)
    fn qr(&self) -> Option<(Self, Self)>
    where
        Self: Sized;
}
