//-------------------------------------------------------------------------
// @file bench_inverse.rs
//
// @date 06/07/20 20:33:40
// @author Martin Noblia
// @email mnoblia@disroot.org
//
// @brief
// Some benchs for speed reference
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
use criterion::{black_box, criterion_group, criterion_main, Criterion};

use static_math::traits::LinearAlgebra;
use static_math::matrix6x6::M66;

fn inverse_test() {

    let m = M66::new([
        [1.0, 1.0, 3.0, 4.0, 9.0, 3.0],
        [10.0, 10.0, 1.0, 2.0, 2.0, 5.0],
        [2.0, 9.0, 6.0, 10.0, 10.0, 9.0],
        [10.0, 9.0, 9.0, 7.0, 3.0, 6.0],
        [7.0, 6.0, 6.0, 2.0, 9.0, 5.0],
        [3.0, 8.0, 1.0, 4.0, 1.0, 5.0],
    ]);

    if let Some(_result) = m.inverse() {
    }
}

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("inverse 6x6", |b| b.iter(|| inverse_test()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
