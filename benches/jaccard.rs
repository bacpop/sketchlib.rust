use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use sketchlib::hashing::nthash_iterator::NtHashIterator;
use sketchlib::sketch::{num_bins, Sketch, BBITS};

// ─── Alignment helper ────────────────────────────────────────────────────────

/// Return a Vec<u64> with at least `align` bytes alignment on its data pointer.
/// The returned Vec may be longer than `src`; use `[offset..offset+src.len()]`.
fn to_aligned_u64(src: &[u64], align: usize) -> (Vec<u64>, usize) {
    let elem = std::mem::size_of::<u64>();
    let overhead = align / elem;
    let mut buf = vec![0u64; src.len() + overhead];
    let addr = buf.as_ptr() as usize;
    let aligned = (addr + align - 1) / align * align;
    let offset = (aligned - addr) / elem;
    buf[offset..offset + src.len()].copy_from_slice(src);
    (buf, offset)
}

fn to_aligned_u16(src: &[u16], align: usize) -> (Vec<u16>, usize) {
    let elem = std::mem::size_of::<u16>();
    let overhead = align / elem;
    let mut buf = vec![0u16; src.len() + overhead];
    let addr = buf.as_ptr() as usize;
    let aligned = (addr + align - 1) / align * align;
    let offset = (aligned - addr) / elem;
    buf[offset..offset + src.len()].copy_from_slice(src);
    (buf, offset)
}

// ─── Scalar / auto-vec baselines ─────────────────────────────────────────────

fn jaccard_u16(a: &[u16], b: &[u16]) -> u32 {
    std::iter::zip(a, b)
        .map(|(a, b)| (a != b) as u32)
        .sum()
}

fn jaccard_packed_u64(a: &[u64], b: &[u64]) -> u32 {
    a.chunks_exact(BBITS as usize)
        .zip(b.chunks_exact(BBITS as usize))
        .map(|(chunk1, chunk2)| {
            let mut bits: u64 = !0;
            chunk1.iter().zip(chunk2.iter()).for_each(|(&s1, &s2)| {
                bits &= !(s1 ^ s2);
            });
            bits.count_ones()
        })
        .sum()
}

// ─── ARM NEON (aarch64) ───────────────────────────────────────────────────────

/// NEON: NEON NOT (mvn.16b) + AND-fold — no scalar round-trip for NOT.
#[cfg(target_arch = "aarch64")]
fn jaccard_neon(a: &[u64], b: &[u64]) -> u32 {
    unsafe { jaccard_neon_inner(a, b) }
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn jaccard_neon_inner(a: &[u64], b: &[u64]) -> u32 {
    use std::arch::aarch64::*;
    let n_chunks = a.len() / (BBITS as usize);
    let mut total = 0u32;
    for i in 0..n_chunks {
        let base = i * BBITS as usize;
        let ap = a.as_ptr().add(base) as *const u8;
        let bp = b.as_ptr().add(base) as *const u8;
        let a0 = vld1q_u8(ap);         let a1 = vld1q_u8(ap.add(16));
        let a2 = vld1q_u8(ap.add(32)); let a3 = vld1q_u8(ap.add(48));
        let a4 = vld1q_u8(ap.add(64)); let a5 = vld1q_u8(ap.add(80));
        let a6 = vld1q_u8(ap.add(96)); let a7 = vld1q_u8(ap.add(112));
        let b0 = vld1q_u8(bp);         let b1 = vld1q_u8(bp.add(16));
        let b2 = vld1q_u8(bp.add(32)); let b3 = vld1q_u8(bp.add(48));
        let b4 = vld1q_u8(bp.add(64)); let b5 = vld1q_u8(bp.add(80));
        let b6 = vld1q_u8(bp.add(96)); let b7 = vld1q_u8(bp.add(112));
        let x0 = veorq_u8(a0, b0); let x1 = veorq_u8(a1, b1);
        let x2 = veorq_u8(a2, b2); let x3 = veorq_u8(a3, b3);
        let x4 = veorq_u8(a4, b4); let x5 = veorq_u8(a5, b5);
        let x6 = veorq_u8(a6, b6); let x7 = veorq_u8(a7, b7);
        let or01 = vorrq_u8(x0, x1); let or23 = vorrq_u8(x2, x3);
        let or45 = vorrq_u8(x4, x5); let or67 = vorrq_u8(x6, x7);
        let or0123 = vorrq_u8(or01, or23);
        let or4567 = vorrq_u8(or45, or67);
        let or_all = vorrq_u8(or0123, or4567);
        let not_or = vmvnq_u8(or_all);
        let folded = vand_u8(vget_low_u8(not_or), vget_high_u8(not_or));
        total += vaddv_u8(vcnt_u8(folded)) as u32;
    }
    total
}

/// NEON unroll-2: process two independent chunks per loop for ILP.
#[cfg(target_arch = "aarch64")]
fn jaccard_neon_unroll2(a: &[u64], b: &[u64]) -> u32 {
    unsafe { jaccard_neon_unroll2_inner(a, b) }
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn jaccard_neon_unroll2_inner(a: &[u64], b: &[u64]) -> u32 {
    use std::arch::aarch64::*;
    let chunk = BBITS as usize;
    let n_chunks = a.len() / chunk;
    let n_pairs = n_chunks / 2;
    let mut total = 0u32;

    #[inline(always)]
    unsafe fn one_chunk(ap: *const u8, bp: *const u8) -> u8 {
        let a0 = vld1q_u8(ap);         let a1 = vld1q_u8(ap.add(16));
        let a2 = vld1q_u8(ap.add(32)); let a3 = vld1q_u8(ap.add(48));
        let a4 = vld1q_u8(ap.add(64)); let a5 = vld1q_u8(ap.add(80));
        let a6 = vld1q_u8(ap.add(96)); let a7 = vld1q_u8(ap.add(112));
        let b0 = vld1q_u8(bp);         let b1 = vld1q_u8(bp.add(16));
        let b2 = vld1q_u8(bp.add(32)); let b3 = vld1q_u8(bp.add(48));
        let b4 = vld1q_u8(bp.add(64)); let b5 = vld1q_u8(bp.add(80));
        let b6 = vld1q_u8(bp.add(96)); let b7 = vld1q_u8(bp.add(112));
        let x0=veorq_u8(a0,b0); let x1=veorq_u8(a1,b1);
        let x2=veorq_u8(a2,b2); let x3=veorq_u8(a3,b3);
        let x4=veorq_u8(a4,b4); let x5=veorq_u8(a5,b5);
        let x6=veorq_u8(a6,b6); let x7=veorq_u8(a7,b7);
        let or01=vorrq_u8(x0,x1); let or23=vorrq_u8(x2,x3);
        let or45=vorrq_u8(x4,x5); let or67=vorrq_u8(x6,x7);
        let or_all=vorrq_u8(vorrq_u8(or01,or23),vorrq_u8(or45,or67));
        let not_or=vmvnq_u8(or_all);
        vaddv_u8(vcnt_u8(vand_u8(vget_low_u8(not_or),vget_high_u8(not_or))))
    }

    for i in 0..n_pairs {
        let b0 = i * 2 * chunk;
        total += one_chunk(a.as_ptr().add(b0) as _, b.as_ptr().add(b0) as _) as u32;
        total += one_chunk(a.as_ptr().add(b0+chunk) as _, b.as_ptr().add(b0+chunk) as _) as u32;
    }
    if n_chunks % 2 != 0 {
        let base = (n_chunks - 1) * chunk;
        total += one_chunk(a.as_ptr().add(base) as _, b.as_ptr().add(base) as _) as u32;
    }
    total
}

/// NEON u16: compare 16 pairs per iteration, accumulate in u16×8 register.
/// Counts EQUAL pairs (same semantics as jaccard_index samebits).
#[cfg(target_arch = "aarch64")]
fn jaccard_u16_neon_wide(a: &[u16], b: &[u16]) -> u32 {
    unsafe { jaccard_u16_neon_wide_inner(a, b) }
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn jaccard_u16_neon_wide_inner(a: &[u16], b: &[u16]) -> u32 {
    use std::arch::aarch64::*;
    let n = a.len();
    // Production data: n = sketchsize64 × 64, always a multiple of 64.
    debug_assert_eq!(n % 16, 0, "u16 length must be a multiple of 16");
    let n16 = n / 16;
    let mut vacc = vdupq_n_u16(0);
    for i in 0..n16 {
        let ap = a.as_ptr().add(i * 16);
        let bp = b.as_ptr().add(i * 16);
        let ones0 = vshrq_n_u16(vceqq_u16(vld1q_u16(ap),     vld1q_u16(bp)),     15);
        let ones1 = vshrq_n_u16(vceqq_u16(vld1q_u16(ap.add(8)), vld1q_u16(bp.add(8))), 15);
        vacc = vaddq_u16(vaddq_u16(vacc, ones0), ones1);
    }
    vaddvq_u16(vacc) as u32
}

// ─── Intel AVX2 (x86_64) ─────────────────────────────────────────────────────

/// AVX2 packed u64: 256-bit registers → 4 u64 per load, depth-2 OR tree.
/// Dispatches via runtime feature detection.
#[cfg(target_arch = "x86_64")]
fn jaccard_packed_avx2(a: &[u64], b: &[u64]) -> u32 {
    if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("popcnt") {
        unsafe { jaccard_packed_avx2_inner(a, b) }
    } else {
        jaccard_packed_u64(a, b)
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
unsafe fn jaccard_packed_avx2_inner(a: &[u64], b: &[u64]) -> u32 {
    use std::arch::x86_64::*;

    let n_chunks = a.len() / (BBITS as usize);
    let mut total = 0u32;
    // Use XOR with all-ones for NOT (AVX2 has no vpnot)
    let all_ones = _mm256_set1_epi8(-1i8);

    for i in 0..n_chunks {
        let base = i * BBITS as usize;
        // Each __m256i holds 4 u64 values (32 bytes).
        // 4 registers × 4 u64 = 16 u64 = one full BBITS chunk.
        let ap = a.as_ptr().add(base) as *const __m256i;
        let bp = b.as_ptr().add(base) as *const __m256i;

        let a0 = _mm256_loadu_si256(ap);
        let a1 = _mm256_loadu_si256(ap.add(1));
        let a2 = _mm256_loadu_si256(ap.add(2));
        let a3 = _mm256_loadu_si256(ap.add(3));

        let b0 = _mm256_loadu_si256(bp);
        let b1 = _mm256_loadu_si256(bp.add(1));
        let b2 = _mm256_loadu_si256(bp.add(2));
        let b3 = _mm256_loadu_si256(bp.add(3));

        // XOR all 4 pairs
        let x0 = _mm256_xor_si256(a0, b0);
        let x1 = _mm256_xor_si256(a1, b1);
        let x2 = _mm256_xor_si256(a2, b2);
        let x3 = _mm256_xor_si256(a3, b3);

        // OR-reduce: depth-2 balanced tree (4 → 2 → 1)
        let or01 = _mm256_or_si256(x0, x1);
        let or23 = _mm256_or_si256(x2, x3);
        let or_all = _mm256_or_si256(or01, or23);

        // NOT (XOR with all-ones since vpnot doesn't exist)
        let not_or = _mm256_xor_si256(or_all, all_ones);

        // Fold 256 → 128: AND upper and lower 128-bit halves.
        // Each half holds contributions from different bit-planes; AND combines them.
        let hi128 = _mm256_extracti128_si256(not_or, 1); // upper 128 bits
        let lo128 = _mm256_castsi256_si128(not_or);       // lower 128 bits (free op)
        let r128 = _mm_and_si128(lo128, hi128);

        // Fold 128 → 64: AND upper and lower 64-bit halves
        let hi64 = _mm_unpackhi_epi64(r128, r128); // copy upper 64 to lower lane
        let r64 = _mm_and_si128(r128, hi64);

        // Extract to u64 and popcount
        let val = _mm_cvtsi128_si64(r64) as u64;
        total += val.count_ones(); // compiles to popcnt with target_feature = "popcnt"
    }
    total
}

/// AVX2 packed u64 unroll-2: process two independent chunks per iteration.
#[cfg(target_arch = "x86_64")]
fn jaccard_packed_avx2_unroll2(a: &[u64], b: &[u64]) -> u32 {
    if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("popcnt") {
        unsafe { jaccard_packed_avx2_unroll2_inner(a, b) }
    } else {
        jaccard_packed_u64(a, b)
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
unsafe fn jaccard_packed_avx2_unroll2_inner(a: &[u64], b: &[u64]) -> u32 {
    use std::arch::x86_64::*;

    let chunk = BBITS as usize;
    let n_chunks = a.len() / chunk;
    let n_pairs = n_chunks / 2;
    let mut total = 0u32;
    let all_ones = _mm256_set1_epi8(-1i8);

    #[inline(always)]
    unsafe fn one_chunk(ap: *const __m256i, bp: *const __m256i, all_ones: __m256i) -> u32 {
        use std::arch::x86_64::*;
        let a0 = _mm256_loadu_si256(ap);
        let a1 = _mm256_loadu_si256(ap.add(1));
        let a2 = _mm256_loadu_si256(ap.add(2));
        let a3 = _mm256_loadu_si256(ap.add(3));
        let b0 = _mm256_loadu_si256(bp);
        let b1 = _mm256_loadu_si256(bp.add(1));
        let b2 = _mm256_loadu_si256(bp.add(2));
        let b3 = _mm256_loadu_si256(bp.add(3));
        let or_all = _mm256_or_si256(
            _mm256_or_si256(_mm256_xor_si256(a0,b0), _mm256_xor_si256(a1,b1)),
            _mm256_or_si256(_mm256_xor_si256(a2,b2), _mm256_xor_si256(a3,b3)),
        );
        let not_or = _mm256_xor_si256(or_all, all_ones);
        let hi128 = _mm256_extracti128_si256(not_or, 1);
        let lo128 = _mm256_castsi256_si128(not_or);
        let r128 = _mm_and_si128(lo128, hi128);
        let r64 = _mm_and_si128(r128, _mm_unpackhi_epi64(r128, r128));
        (_mm_cvtsi128_si64(r64) as u64).count_ones()
    }

    for i in 0..n_pairs {
        let base0 = i * 2 * chunk;
        let base1 = base0 + chunk;
        let ap0 = a.as_ptr().add(base0) as *const __m256i;
        let bp0 = b.as_ptr().add(base0) as *const __m256i;
        let ap1 = a.as_ptr().add(base1) as *const __m256i;
        let bp1 = b.as_ptr().add(base1) as *const __m256i;
        total += one_chunk(ap0, bp0, all_ones);
        total += one_chunk(ap1, bp1, all_ones);
    }
    if n_chunks % 2 != 0 {
        let base = (n_chunks - 1) * chunk;
        let ap = a.as_ptr().add(base) as *const __m256i;
        let bp = b.as_ptr().add(base) as *const __m256i;
        total += one_chunk(ap, bp, all_ones);
    }
    total
}

/// AVX2 u16: compare 32 u16 pairs per iteration with u16×16 accumulator.
/// Counts EQUAL pairs. Requires SSSE3 for the final hadd-based reduction.
#[cfg(target_arch = "x86_64")]
fn jaccard_u16_avx2(a: &[u16], b: &[u16]) -> u32 {
    if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("ssse3") {
        unsafe { jaccard_u16_avx2_inner(a, b) }
    } else {
        // Scalar fallback: count equal pairs
        a.iter().zip(b.iter()).map(|(x, y)| (x == y) as u32).sum()
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,ssse3")]
unsafe fn jaccard_u16_avx2_inner(a: &[u16], b: &[u16]) -> u32 {
    use std::arch::x86_64::*;

    let n = a.len();
    // Production data: n = sketchsize64 × 64, always a multiple of 64 ≥ 32.
    debug_assert_eq!(n % 32, 0, "u16 length must be a multiple of 32 for AVX2");
    let n32 = n / 32;

    // u16×16 accumulator.
    // Max per lane per iteration = 2 (two add_epi16 calls, each adding 0 or 1).
    // For sketch_size=10000: n32 ≤ 314 → max 628 per lane ≪ 65535. ✓
    let mut vacc = _mm256_setzero_si256();

    for i in 0..n32 {
        let ap = a.as_ptr().add(i * 32) as *const __m256i;
        let bp = b.as_ptr().add(i * 32) as *const __m256i;

        let va0 = _mm256_loadu_si256(ap);
        let va1 = _mm256_loadu_si256(ap.add(1));
        let vb0 = _mm256_loadu_si256(bp);
        let vb1 = _mm256_loadu_si256(bp.add(1));

        // Compare: 0xFFFF if equal, 0 otherwise
        // Shift right by 15 to get 0 or 1 per u16 lane
        let ones0 = _mm256_srli_epi16(_mm256_cmpeq_epi16(va0, vb0), 15);
        let ones1 = _mm256_srli_epi16(_mm256_cmpeq_epi16(va1, vb1), 15);

        vacc = _mm256_add_epi16(_mm256_add_epi16(vacc, ones0), ones1);
    }

    // Horizontal reduce: sum all 16 u16 lanes.
    // Step 1: fold 256 → 128 (add upper and lower lanes)
    let hi128 = _mm256_extracti128_si256(vacc, 1);
    let lo128 = _mm256_castsi256_si128(vacc);
    let sum128 = _mm_add_epi16(lo128, hi128); // 8 u16 values

    // Step 2: three hadd_epi16 (SSSE3 phaddw) reduce 8 → 4 → 2 → 1 u16
    let s1 = _mm_hadd_epi16(sum128, sum128);
    let s2 = _mm_hadd_epi16(s1, s1);
    let s3 = _mm_hadd_epi16(s2, s2);

    (_mm_cvtsi128_si32(s3) & 0xFFFF) as u32
}

// ─── Benchmark ────────────────────────────────────────────────────────────────

fn bench_jaccard(c: &mut Criterion) {
    let mut group = c.benchmark_group("jaccard");

    let files1 = vec!["tests/test_files_in/14412_3#82.contigs_velvet.fa.gz".to_string()];
    let files2 = vec!["tests/test_files_in/14412_3#84.contigs_velvet.fa.gz".to_string()];
    let mut it1 = NtHashIterator::new(&files1, true, 0);
    let mut it2 = NtHashIterator::new(&files2, true, 0);

    for sketch_size in [100u64, 1000, 10000] {
        let (_sketchsize64, signs_size, usigs_size) = num_bins(sketch_size);

        let (signs1, _) = Sketch::get_signs(&mut it1[0], 21, &mut None, signs_size);
        let (signs2, _) = Sketch::get_signs(&mut it2[0], 21, &mut None, signs_size);

        let mut usigs1 = vec![0u64; usigs_size as usize];
        let mut usigs2 = vec![0u64; usigs_size as usize];
        Sketch::fill_usigs(&mut usigs1, &signs1);
        Sketch::fill_usigs(&mut usigs2, &signs2);

        let signs1_u16: Vec<u16> = signs1.iter().map(|&s| s as u16).collect();
        let signs2_u16: Vec<u16> = signs2.iter().map(|&s| s as u16).collect();

        // ── Alignment report (printed once per sketch_size, not per iteration) ──
        let u64_align = usigs1.as_ptr() as usize % 32;
        let u16_align = signs1_u16.as_ptr() as usize % 32;
        println!(
            "sketch_size={sketch_size}: usigs ptr%32={u64_align}  signs_u16 ptr%32={u16_align}  \
             (0 = already 32-byte aligned)"
        );

        // ── 32-byte-aligned copies for alignment benchmarks ──────────────────
        let (abuf1, aoff1) = to_aligned_u64(&usigs1, 32);
        let (abuf2, aoff2) = to_aligned_u64(&usigs2, 32);
        let aligned_usigs1 = &abuf1[aoff1..aoff1 + usigs1.len()];
        let aligned_usigs2 = &abuf2[aoff2..aoff2 + usigs2.len()];
        debug_assert_eq!(aligned_usigs1.as_ptr() as usize % 32, 0);
        debug_assert_eq!(aligned_usigs2.as_ptr() as usize % 32, 0);

        let (sbuf1, soff1) = to_aligned_u16(&signs1_u16, 32);
        let (sbuf2, soff2) = to_aligned_u16(&signs2_u16, 32);
        let aligned_signs1 = &sbuf1[soff1..soff1 + signs1_u16.len()];
        let aligned_signs2 = &sbuf2[soff2..soff2 + signs2_u16.len()];
        debug_assert_eq!(aligned_signs1.as_ptr() as usize % 32, 0);
        debug_assert_eq!(aligned_signs2.as_ptr() as usize % 32, 0);

        // ── Correctness check (runs once per sketch_size) ────────────────────
        let ref_count = jaccard_packed_u64(&usigs1, &usigs2);
        #[cfg(target_arch = "aarch64")]
        assert_eq!(jaccard_neon(&usigs1, &usigs2), ref_count, "neon mismatch");
        #[cfg(target_arch = "aarch64")]
        assert_eq!(jaccard_neon_unroll2(&usigs1, &usigs2), ref_count, "neon_unroll2 mismatch");
        // u16 compares the lower 16 bits only; verify self-consistency across u16 variants
        let _u16_ref = jaccard_u16_neon_wide_or_scalar(&signs1_u16, &signs2_u16);
        #[cfg(target_arch = "x86_64")]
        {
            assert_eq!(jaccard_packed_avx2(&usigs1, &usigs2), ref_count, "avx2 mismatch");
            assert_eq!(jaccard_packed_avx2_unroll2(&usigs1, &usigs2), ref_count, "avx2_unroll2 mismatch");
            assert_eq!(jaccard_u16_avx2(&signs1_u16, &signs2_u16), _u16_ref, "u16_avx2 mismatch");
        }

        // ── packed_u64 variants ───────────────────────────────────────────────
        group.bench_with_input(
            BenchmarkId::new("packed_u64", sketch_size),
            &sketch_size,
            |b, _| b.iter(|| jaccard_packed_u64(black_box(&usigs1), black_box(&usigs2))),
        );
        group.bench_with_input(
            BenchmarkId::new("packed_u64_aligned32", sketch_size),
            &sketch_size,
            |b, _| b.iter(|| jaccard_packed_u64(black_box(aligned_usigs1), black_box(aligned_usigs2))),
        );

        #[cfg(target_arch = "aarch64")]
        group.bench_with_input(
            BenchmarkId::new("neon_intrinsics", sketch_size),
            &sketch_size,
            |b, _| b.iter(|| jaccard_neon(black_box(&usigs1), black_box(&usigs2))),
        );
        #[cfg(target_arch = "aarch64")]
        group.bench_with_input(
            BenchmarkId::new("neon_aligned32", sketch_size),
            &sketch_size,
            |b, _| b.iter(|| jaccard_neon(black_box(aligned_usigs1), black_box(aligned_usigs2))),
        );
        #[cfg(target_arch = "aarch64")]
        group.bench_with_input(
            BenchmarkId::new("neon_unroll2", sketch_size),
            &sketch_size,
            |b, _| b.iter(|| jaccard_neon_unroll2(black_box(&usigs1), black_box(&usigs2))),
        );

        // ── AVX2 packed_u64 variants ──────────────────────────────────────────
        #[cfg(target_arch = "x86_64")]
        group.bench_with_input(
            BenchmarkId::new("avx2_packed", sketch_size),
            &sketch_size,
            |b, _| b.iter(|| jaccard_packed_avx2(black_box(&usigs1), black_box(&usigs2))),
        );
        #[cfg(target_arch = "x86_64")]
        group.bench_with_input(
            BenchmarkId::new("avx2_packed_aligned32", sketch_size),
            &sketch_size,
            |b, _| b.iter(|| jaccard_packed_avx2(black_box(aligned_usigs1), black_box(aligned_usigs2))),
        );
        #[cfg(target_arch = "x86_64")]
        group.bench_with_input(
            BenchmarkId::new("avx2_packed_unroll2", sketch_size),
            &sketch_size,
            |b, _| b.iter(|| jaccard_packed_avx2_unroll2(black_box(&usigs1), black_box(&usigs2))),
        );

        // ── u16 variants ──────────────────────────────────────────────────────
        group.bench_with_input(
            BenchmarkId::new("u16_signs", sketch_size),
            &sketch_size,
            |b, _| b.iter(|| jaccard_u16(black_box(&signs1_u16), black_box(&signs2_u16))),
        );
        #[cfg(target_arch = "aarch64")]
        group.bench_with_input(
            BenchmarkId::new("u16_neon_wide", sketch_size),
            &sketch_size,
            |b, _| b.iter(|| jaccard_u16_neon_wide(black_box(&signs1_u16), black_box(&signs2_u16))),
        );
        #[cfg(target_arch = "aarch64")]
        group.bench_with_input(
            BenchmarkId::new("u16_neon_wide_aligned32", sketch_size),
            &sketch_size,
            |b, _| b.iter(|| jaccard_u16_neon_wide(black_box(aligned_signs1), black_box(aligned_signs2))),
        );
        #[cfg(target_arch = "x86_64")]
        group.bench_with_input(
            BenchmarkId::new("u16_avx2", sketch_size),
            &sketch_size,
            |b, _| b.iter(|| jaccard_u16_avx2(black_box(&signs1_u16), black_box(&signs2_u16))),
        );
        #[cfg(target_arch = "x86_64")]
        group.bench_with_input(
            BenchmarkId::new("u16_avx2_aligned32", sketch_size),
            &sketch_size,
            |b, _| b.iter(|| jaccard_u16_avx2(black_box(aligned_signs1), black_box(aligned_signs2))),
        );
    }

    group.finish();
}

/// Helper for correctness check: call u16_neon_wide on aarch64, scalar on x86
fn jaccard_u16_neon_wide_or_scalar(a: &[u16], b: &[u16]) -> u32 {
    #[cfg(target_arch = "aarch64")]
    return jaccard_u16_neon_wide(a, b);
    #[cfg(not(target_arch = "aarch64"))]
    a.iter().zip(b.iter()).map(|(x, y)| (x == y) as u32).sum()
}

criterion_group!(benches, bench_jaccard);
criterion_main!(benches);
