use criterion::{black_box, criterion_group, criterion_main, Criterion};
use sketchlib::distances::jaccard::{jaccard_index, jaccard_index2, jaccard_index3, jaccard_index4};
use sketchlib::sketch::multisketch::MultiSketch;
use sketchlib::sketch::sketch_datafile::SketchArrayReader;
use simd_aligned::{VecSimd, arch::u16x16};

pub fn criterion_benchmark(c: &mut Criterion) {
    let prefix = "/Users/jlees/Documents/EBI/sketchlib_rust/tests/test_files_in/legacy_db";
    let mut sketches = MultiSketch::load_metadata(prefix).unwrap();
    sketches.read_sketch_data(prefix);
    let sketch1 = sketches.get_sketch_slice(0, 0);
    let sketch2 = sketches.get_sketch_slice(1, 0);
    let sketchsize64 = sketches.sketchsize64;
    c.bench_function("jaccard1", |b| b.iter(|| black_box(jaccard_index(sketch1, sketch2, sketchsize64))));
    c.bench_function("jaccard4", |b| b.iter(|| black_box(jaccard_index4(sketch1, sketch2, sketchsize64))));

    let skq_filename = "/Users/jlees/Documents/EBI/sketchlib_rust/tests/test_files_in/inverted.skq";
    let sketchsize = 1000;
    let mut skq_reader = SketchArrayReader::open(
        skq_filename,
        false,
        1,
        1,
        sketchsize,
    );
    let skq_bins =
        skq_reader.read_all_from_skq(None);

    let sketch1= &skq_bins[0..sketchsize];
    let sketch2 = &skq_bins[sketchsize..(2 * sketchsize)];
    c.bench_function("jaccard2", |b| b.iter(|| black_box(jaccard_index2(sketch1, sketch2, sketchsize))));


    let mut s1 = VecSimd::<u16x16>::with(0_u16, sketchsize);
    s1.flat_mut().clone_from_slice(&skq_bins[0..sketchsize]);
    let mut s2 = VecSimd::<u16x16>::with(0_u16, sketchsize);
    s2.flat_mut().clone_from_slice(&skq_bins[sketchsize..(2 * sketchsize)]);

    c.bench_function("jaccard3", |b| b.iter(|| black_box(jaccard_index3(&s1, &s2, sketchsize))));

}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);