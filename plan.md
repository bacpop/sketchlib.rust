# CUDA self-distance planning

## Current round: shared dense/kNN path and completeness

Status: initial implementation, remote release build, smoke test, and agreed
integration validation complete. The follow-up removal of redundant dense text
batching is implemented locally, with remote verification pending.
The user authorised this round with "Yes, implement it".
Remote commands were approved individually and run sequentially as `jlees`.
The completed tiling phase and its validation remain recorded below.

Requested scope:

- Use a common synchronous GPU tile pipeline for dense self-distances and kNN,
  selecting the consumer with `knn: Option<usize>`. Dense output goes to a writer;
  kNN continues to populate heaps.
- Upload completeness values once, using f32 on the device, and apply completeness
  correction consistently across dense and kNN modes.
- Retain the current user-controlled block-size default of 10000 and the GPU
  minimum dataset size of 256. No automatic GPU-memory sizing is requested.
- Keep execution synchronous. Concurrent tiles, Tokio scheduling, pinned sketch
  buffers, parallel consumers, prefiltering, shared-memory query caching, and
  reference/query cross-comparisons remain deferred. Core/accessory regression
  remains deferred from the previous phase.
- Preserve canonical transfer/allocation helpers unless a necessary change is
  separately approved; type-only variants are permitted. Avoid unnecessary
  validation and extra wrapper abstractions.
- Run compilation and runtime tests only on the remote machine, as `jlees`.
  Ask before every remote command and run one command at a time. A reviewable
  Python runner may group sequential smoke and CPU/GPU integration tests into
  one approved remote invocation; do not add tests to the Cargo harness.

Resolved decisions for this round:

- Dense GPU output adopts the CPU streaming contract: the same tab-separated
  line format, each unordered sample pair exactly once, and no guaranteed global
  record ordering. The user explicitly accepted tile-by-tile streaming to keep
  memory bounded; do not assemble or reorder a complete dense matrix.
- Replace the two GPU entry points with one writer-based function,
  `self_dists_stream(writer, ..., knn: Option<usize>, ...) -> Result<()>`.
  `None` streams dense output; `Some(k)` retains the existing heap selection and
  sparse output format. Construct the sparse matrix and perform its final write
  inside this function. The user accepted the slightly less tidy caller in
  exchange for avoiding additional optional arguments. CPU entry points remain
  unchanged.
- Preserve the CPU completeness rule on the GPU: after random-match correction,
  apply correction when `c_i * c_j >= completeness_cutoff`, cap the corrected
  Jaccard index at 1, and then derive distance or ANI. Use f32 completeness
  values, cutoff, and correction arithmetic; retain the existing arithmetic for
  uncorrected Jaccard and ANI. Upload completeness once for all tiles and use
  each tile's global sample offsets to select values from it.
- Document that f32 rounding can choose a different correction branch from CPU
  f64 for inputs extremely close to the cutoff. The user explicitly excluded
  testing that precision boundary; do not add a boundary test for it.
- Defer progress reporting entirely until the asynchronous tile-launch phase,
  as explicitly requested by the user. Do not add a polling thread, managed
  counter, device atomic updates, capability checks, or progress/quiet tests in
  this round. The proposed polling-thread design was not accepted; future work
  should reconsider progress alongside asynchronous scheduling.

Findings used for this round:

- Before this round, GPU dense calculation allocated all distances, while the
  kNN route already had the bounded synchronous tile pipeline now shared by
  both output modes.
- CPU completeness correction runs after random-match correction. Its gate is
  `c_i * c_j >= completeness_cutoff` (default 0.64), rather than an individual
  threshold on each sample. The corrected Jaccard index is capped at 1 before
  conversion to Jaccard distance or ANI. Missing file entries are already filled
  with completeness 1 by the existing reader.
- Dense ANI needs the existing kernel's raw-ANI mode. kNN must retain its current
  inverse-ANI ranking and sparse finalisation, avoiding an additional float
  conversion through inverse ANI for dense output.
- Dense formatting reuses the tile's pair decoder and global offsets and writes
  each record directly to the supplied writer. `set_ostream()` already wraps
  file and stdout destinations in `std::io::BufWriter`; do not add a second
  text batching layer in `consume()`. Writer errors propagate through the
  unified entry point.

Agreed remote validation for this round:

- Compile using `CUDA_HOME=/usr/local/cuda-13.4 cargo oxide build` and run the
  supplied 1000-sample GPU kNN smoke command against the resulting executable.
- Use a 256-sample subset for the integration matrix: dense and kNN output,
  Jaccard and ANI, with no completeness file and with a synthetic mixed
  completeness file. Choose values comfortably above/below the default product
  cutoff, exercising both corrected and uncorrected pairs without testing the
  excluded precision boundary.
- Compare CPU results with both single-tile and forced-multiple-tile GPU output.
  Check dense unordered-pair coverage, exact pair count, no duplicates or
  self-pairs, finite values, and absolute distance tolerance 1e-6. Keep the
  previously agreed sparse structural checks and CPU-verified cutoff-tie rule;
  dense CPU output also supplies distances for verifying alternative neighbours.
- Include a singleton remainder in a dense tiled run to exercise the newly
  shared consumer and skipped empty diagonal tile. Verify that mixed
  completeness actually changes results for some eligible pairs.
- Keep the runner outside the Cargo test harness and execute its subprocesses
  sequentially in one approved remote invocation. Request separate approval for
  code copies and compilation as needed; never run a remote command silently.

Completed implementation sequence:

1. Replace the standalone GPU dense implementation and GPU kNN entry point with
   `self_dists_stream()`. Keep the existing `DistTile` lifecycle:
   construct, build pipeline, `.sync()`, consume, then start the next tile.
   Initialise sparse storage/heaps only for `Some(k)`; dense mode must not
   construct a `DistanceMatrix` or retain previous tiles' results.
2. Adapt `DistTile::consume()` to the two output modes while retaining the shared
   kernel/host pair decoder. Dense mode writes each record to the supplied
   buffered writer before moving to the next tile. kNN still adds each unordered pair to
   both heaps and uses the existing shared finalisation before writing.
   Select raw ANI in dense mode and inverse ANI for sparse ranking.
3. Convert the existing optional host completeness values to f32 once, upload
   them once using a permitted type-only `h2d_f32()` variant, and share the
   device buffer across tile operations. Retain all buffer owners through
   `.sync()`, as in the validated pipeline. Without a completeness file, retain
   the existing uncorrected arithmetic. Add only the kernel arguments and
   correction logic needed for the agreed gate/formula, using global sample
   offsets. Preserve the canonical `h2d()`, `d2h()`, and `zeros()` functions.
4. Update CLI self-distance GPU dispatch to call the unified writer-based entry
   point and pass the block limit in both modes. Preserve CPU routes and
   existing CLI kNN clamping. Update block-limit help text to cover dense and
   sparse GPU self-distances. Keep the default 10000, dataset minimum 256,
   current-format restriction, and deferred core/accessory regression.
5. Format and inspect the local diff without local compilation or tests.
   Prepare the one-off sequential remote runner, request approval for every
   remote command, then compile and run the agreed validation. Resolve failures
   in scope and record actual outcomes and any remaining limits in this plan.

No new wrapper hierarchy or generic consumer framework was introduced. Changes
are confined to the GPU implementation, required dispatch/help updates,
planning/glossary documentation, and the one-off remote validation script. No
ADR was needed for these readily reversible interface and execution choices.

Implementation and validation record for this round:

- `gpu::self_dists_stream()` is the single writer-based entry point. Dense and
  sparse calculations use the same `DistTile` pipeline and shared pair decoder.
  The sparse matrix is constructed and written inside this function. CPU
  calculation functions are unchanged; the caller retains its kNN clamping.
- Dense mode does not allocate a `DistanceMatrix` or heaps. It writes records
  directly with `writeln!(writer, ...)`, relying on the existing `BufWriter`.
  The user requested removal of the initial redundant 64 KiB `String` batching
  layer, including its manual flush logic and unused formatting-trait import.
  Dense ANI is calculated directly, while sparse ANI retains inverse-ANI
  ranking and the existing shared heap finalisation.
- Completeness is converted to f32 and uploaded once with `h2d_f32()`, a
  type-only variant of the canonical transfer helper. Tiles share this buffer
  with Arc and use global sample offsets for lookup. When correction is absent,
  a single unused placeholder value supplies a valid typed launch buffer and
  the uncorrected arithmetic remains unchanged.
- The correction gate, formula, and cap match the agreed CPU semantics using
  f32 arithmetic. The public entry-point documentation and kernel comment
  describe possible cutoff-branch differences from f64. No precision-boundary
  test was added or run, as requested.
- The canonical `h2d()`, `d2h()`, and `zeros()` functions are unchanged. Device
  buffers remain retained through `.sync()`. Progress reporting and every other
  deferred feature remain unimplemented.
- Local Rust formatting and `git diff --check` passed. No local compilation or
  runtime tests ran, and no tests were added to the Cargo harness.

The approved remote build used:

```sh
CUDA_HOME=/usr/local/cuda-13.4 /home/jlees/.cargo/bin/cargo oxide build
```

It completed successfully in 19.15s with the release profile. The remote
toolchain remains `1.99.1`; no manifest, toolchain, or cuda-oxide revision pinning
changes were needed in this round. Runtime validation used that newly built
`target/release/sketchlib` executable on GPU 0.

The GPU kNN smoke run on `1k_subset` passed. The integration runner selected
256 names spread through the collection and reversed their order to exercise
logical sample indices independently of metadata order. Mixed completeness
cycled through 0.95, 0.85, 0.50, and 1.0, away from the default cutoff boundary.

| Completeness | Metric | Dense maximum absolute error | kNN maximum absolute error | Verified cutoff substitutions per kNN run |
| --- | --- | ---: | ---: | ---: |
| Absent | Jaccard | 0 | 0 | 0 |
| Absent | ANI | 0 | 0 | 0 |
| Mixed | Jaccard | 1.4e-7 | 1.3e-7 | 1 |
| Mixed | ANI | 1e-7 | 1e-7 | 1 |

Every table row passed for both the default single tile and two 128-sample
blocks. An additional mixed-completeness dense Jaccard run with block limit 255
passed with maximum error 1.4e-7, exercising the singleton remainder and skipped
empty diagonal tile. In total, all **17 CPU/GPU comparisons passed** at the
agreed absolute tolerance 1e-6.

Each dense run contained exactly 32640 distinct unordered pairs, with no
self-pairs and finite values. Each kNN run contained 256 rows with ten distinct
non-self neighbours in the correct order; alternative cutoff neighbours were
verified against dense CPU distances. The runner also confirmed that mixed
completeness changes eligible pairs and leaves below-cutoff pairs uncorrected.

The one-off runner is `/private/tmp/sketchlib-cuda-phase2.py` locally. The
successful invocation supplied it to remote `python3 -` over SSH standard input;
an earlier launch found no runner at the requested remote file path and ran no
tests. No Cargo test harness was used. All successful runtime outputs, per-run
logs, the subset/completeness fixtures, and `report.json` remain remotely in
`/tmp/sketchlib-cuda-phase2-sqre1woq/`.

Local logs are `/private/tmp/sketchlib-cuda-phase2-build.log` and
`/private/tmp/sketchlib-cuda-phase2-validation.log`.

## Completed round: CUDA self-distance kNN tiling

Status: implementation, remote build, and planned remote validation complete.
The user confirmed shared understanding and explicitly authorised implementation.

Implementation note: the cached cuda-oxide `DeviceBox` contract requires all
using streams to be synchronised before a buffer is dropped. The user explicitly
approved the additional `d2h()` lifetime fix: return its device buffer alongside
the host values, so it and all intermediate buffers survive pipeline completion.

## Confirmed brief

- Extend the GPU self-distance kNN path. The current GPU entry point is
  `src/gpu/mod.rs::self_dists_knn()`; the CPU implementation has the name
  `self_dists_knn_generic()` used in the brief.
- Generalise the distance kernel to rectangular reference/query comparisons,
  verify sketch interleaving, and consume tile results into per-sample heaps.
- Support single-k Jaccard distance and ANI; defer core/accessory regression.
- Introduce `DistTile`, with launch-range construction in `new()` and methods
  for running the tile pipeline. Keep additional abstractions minimal.
- Run tiles synchronously in this phase, with each tile's GPU work ordered on
  a CUDA stream. Connect the full H2D, interleaving, distance, and D2H pipeline;
  call `.sync()` and consume its results into the heaps before starting the
  next tile.
- Pass a `max_sketch_batch` CLI setting through to the GPU path, defaulting to
  10000. Use one batch when the sample count is at or below that limit. The
  limit is user-controlled; do not calculate block sizes from available GPU
  memory.
  Use the conventional CLI spelling `--max-sketch-batch` and require a positive
  value at the CLI boundary.
- Require at least 256 samples at the GPU entry point.
  This minimum applies to the entire dataset; individual blocks can be shorter.
- Store the existing aligned sketch vector in `MultiSketch` as `Arc<SketchVec>`.
  Each `DistTile` retains an Arc clone and its sample ranges, sharing sketch
  storage without copying host sketch payloads.
- Preserve the canonical `h2d()`, `d2h()`, and `zeros()` helpers. Type-only
  variants are allowed; other changes require the user's agreement.
- Avoid unnecessary input checks, particularly in device code and kernels.
- Leave `self_dists_all_stream()` outside the implementation scope while
  retaining compilation compatibility with shared kernel changes.

## Deferred

- Concurrent tile execution and Tokio dispatch through `run_gpu_batches()`.
  Connecting the complete pipeline for synchronous execution is in scope.
- `PinnedHostBuffer`, progress reporting, `quiet` handling, completeness
  vectors, and completeness correction.
- Core/accessory regression.
- Formal package tests for this work.

## Findings from the current checkout

- The unfinished loop fixes the number of sketch blocks at four. Its final
  block is calculated from the remainder, which is zero when the sample count
  divides evenly. Inclusive endpoint calculations are then used as exclusive
  Rust slice bounds, dropping sketches.
- The distance kernel currently accepts one sketch buffer, despite separate
  reference/query dimensions and strides. Separate tile inputs and matching
  host decoding need to be designed together.
- Interleaving maps host `[sample][k-mer][bin]` storage to device
  `[k-mer][bin][sample]` storage. Launch arguments must use the full number of
  stored u64 bins, including all bit planes.
- CPU kNN ranks ANI using `(1.0_f64 - ANI) as f32` and restores output ANI after
  sorting. GPU ranking and tie handling need to be checked against this path.
- The GPU module contains incomplete pipeline expressions and unresolved
  identifiers, including in the dense path. No compilation has been attempted.
- The rectangular decoding helper assumes `i < n_query` in a debug assertion,
  which does not hold for tall rectangular tiles. The inclusive-triangle row
  helper similarly asserts that the flattened tile index is below the number
  of blocks, rather than below the number of tiles. Neither assumption should
  be carried into the new tile mapping.
- CUDA dependencies currently reference unpinned `NVlabs/cuda-oxide` git
  sources. API compatibility must be checked during the authorised remote build.

## Resolved decisions

- Block sizing: the user controls `max_sketch_batch`, default 10000. There is
  no automatic adjustment for available GPU memory. At the default limit, a
  rectangular tile has at most 100 million f32 outputs (400 MB decimal); a
  strict-upper diagonal tile would have approximately 200 MB of outputs,
  before sketch buffers and the host result copy.
- Synchronous execution: build the complete per-tile H2D/interleave/distance/
  D2H pipeline, execute it with `.sync()`, then consume its results into heaps.
  Complete this lifecycle before starting the next tile. Tokio dispatch and
  overlapping tiles remain deferred.
- Pair coverage: partition sketches into contiguous blocks of at most the
  selected limit, with a shorter final block when necessary. Visit the upper
  triangle of block pairs, including diagonal tiles. Diagonal tiles contain
  compact strict-upper sample comparisons; off-diagonal tiles contain the full
  rectangle, with reference samples as rows and query samples as columns.
  Calculate each unordered sample pair once and insert it into both samples'
  heaps. A one-sample diagonal tile has no comparisons and is skipped.
- Consumer indexing: reuse the kernel's flattened-index-to-local-(i,j) logic
  when decoding each tile's output. Add the tile's reference/query sample
  offsets to select the global (i,j) and (j,i) heap entries. Keep host and device
  mapping consistent rather than introducing independent formulas.
- The 256-sample minimum applies to the complete GPU dataset. Final blocks may
  contain fewer than 256 samples, including a single sample.
- Host ownership: `MultiSketch` owns `Arc<SketchVec>`, preserving the existing
  aligned sketch allocation. Each tile holds a cheap Arc clone and ranges
  into that allocation. Continue exposing borrowed slices for CPU callers and
  add access to the shared owner for GPU tiles. Keep that owner alive through
  synchronous pipeline completion, including while H2D copies are in flight.
- Supported modes: single-k Jaccard distance and ANI. Core/accessory regression
  remains deferred. For ANI, match CPU heap ranking by computing
  `(1.0_f64 - ANI) as f32`, and restore output ANI after sorting. Share the
  existing heap insertion and final-row conversion logic where required by the
  current GPU TODO, rather than introducing a separate GPU ordering policy.

## Agreed DistTile interface

The user accepted this interface with compatibility for later concurrent tile
execution and a note about future dense streaming:

- One `DistTile` represents one pair of sketch blocks.
- `new()` takes a flattened tile index, sketch metadata/shared storage, and the
  user-selected block limit. It calculates block coordinates, half-open sample
  and bin ranges, diagonal/rectangular output length, and kernel launch ranges.
- `pipeline()` builds the complete `DeviceOperation` chain and returns its
  distance output. It does not synchronise internally or borrow the shared
  heaps. The driver executes that chain with `.sync()` in this phase.
- `consume()` decodes the returned distances using the same local pair mapping
  as the kernel, adds the tile offsets, and inserts into the two sample heaps.
- The driver owns the module and global heaps, visits one tile at a time,
  skips empty diagonal tiles, and finalises the sparse matrix after all tiles.
  Avoid an additional tile manager, batch wrapper, or eagerly retained tile
  collection.
- Future concurrency: independent tile operations can be scheduled through
  `StreamPoolRoundRobin`. Keep each tile and its shared host-storage owner alive
  through completion, and associate its result with the same tile metadata.
  A future task may own the tile, await its pipeline, and return the tile with
  its results. This allows the operation to borrow from its owning tile while
  the enclosing task owns the complete lifetime. Update the global heaps from
  a serial host consumer after completion, allowing GPU work to overlap without
  concurrent heap mutation. Bound the number of in-flight tiles independently
  of the stream-pool size; the pool only assigns streams. These are future
  compatibility requirements, not concurrent execution work for this phase.
- Add a comment in `consume()` that its tile-to-global pair decoding can later
  feed a dense-distance streaming writer instead of updating kNN heaps. Dense
  streaming and its output-order handling remain deferred; do not introduce
  a new generic consumer abstraction in this phase.

## Approved compilation compatibility edits

The user explicitly approved these narrow exceptions to the dense-path and
transfer-helper constraints:

- In `self_dists_all_stream()`, remove the stray `if knn >= n` error block:
  this function has no `knn` argument. Adapt its existing `module.dists()` call
  to pass `&sketch_dev` as both the reference and query inputs to the shared
  two-buffer kernel. Keep its algorithm, original-layout strides, output
  handling, and interleaving/streaming TODOs otherwise as they are.
- In the existing borrowed-u64 `h2d()` signature, explicitly tie the returned
  operation to the lifetime of its borrowed slice by adding `+ '_`:
  `fn h2d(host_data: &[u64]) -> impl DeviceOperation<Output = DeviceBox<[u64]>> + '_`.
  Keep its allocation and copy body as written. This is a Rust lifetime fix
  for the current edition-2021 signature, not an ownership or transfer-policy
  change. Preserve the `d2h()` and f32 `zeros()` bodies; a u64 zero-allocation
  variant is permitted by the original brief.
- Additional approved lifetime fix: `d2h()` returns
  `(Vec<f32>, DeviceBox<[f32]>)` and ends with `value((host, dev))`. Its copy
  remains as written. Carry the intermediate device buffers alongside the host
  result in the pipeline output, then drop those buffers after `.sync()`.
- Replace the unused `run_gpu_batches()` example body, which references
  undefined `prepare_batches()`, `process_batches()`, and `consume_results()`,
  with a clearly marked deferred `todo!()` placeholder. Do not implement the
  Tokio executor or call this placeholder from the synchronous path.
- Add the imports and shared-kernel argument adjustments required for the
  agreed implementation. Resolve any further nontrivial protected-helper or
  dense-path changes with the user if the remote compiler identifies them.

## Agreed comparison policy

The user explicitly accepted this CPU/GPU comparison rule:

- Compare GPU output against CPU output for Jaccard and ANI, and compare
  single-tile and forced-multiple-tile GPU runs.
- Require the expected number of neighbours for each sample, valid sample
  names, no self-pairs, no repeated neighbour within a row, finite values, and
  sorting in the correct direction.
- Compare distances for matching sample pairs with absolute tolerance 1e-6.
- Require matching neighbour identities except for ties at the kNN cutoff
  within that tolerance. Ignore ordering within tied groups. Verify any
  alternative cutoff neighbour against its actual CPU distance, using an
  additional CPU reference run if needed; do not infer a tie solely from the
  GPU-reported value.
- The CPU `SparseJaccard` comparator orders only by distance, with no sample-ID
  tiebreaker. Keep that production comparator as it is; the exception above is
  a validation rule, not a change to ranking behavior.

## Implementation sequence

The user authorised this sequence, now completed:

1. Mark implementation in progress in this plan. Introduce shared aligned
   sketch storage, adapting existing constructors/loaders and providing the
   shared-owner getter while preserving borrowed slice access and file formats.
2. Implement `DistTile::new()` and the shared pair-index decoding. Use
   `ceil(n / max_sketch_batch)` blocks with half-open sample ranges
   `[block * limit, min(n, (block + 1) * limit))`. Enumerate the upper triangle
   of block pairs, including its diagonal, with no empty remainder block.
3. Generalise the distance kernel to separate reference/query inputs and
   strides. Verify the interleaving index calculation against the actual host
   layout, including all bit planes and multiple stored k-mer lengths. Keep
   the necessary output-index guard; avoid redundant device validation.
4. Build `DistTile::pipeline()` using the canonical transfers and the permitted
   u64 zero-allocation variant. Carry the buffers needed by later stages and
   retain host ownership until completion. Implement `consume()` with shared
   pair decoding and two heap insertions, plus the requested dense-streaming
   comment. Share CPU heap/final-row conversion logic. Compute inverse ANI
   before the f32 cast for kNN ranking, while preserving raw ANI semantics for
   the existing dense caller.
5. Thread `--max-sketch-batch` through CLI dispatch into GPU kNN. Replace the
   unfinished tile loop with construct/pipeline/sync/consume processing of one
   tile at a time. Apply only the approved compatibility edits to the dense
   caller, `h2d()` lifetime annotation, and deferred async placeholder.
6. Review the local diff without local compilation or tests. Prepare the code
   transfer and concrete remote validation commands, requesting permission
   before each remote command. Run the agreed validation remotely, resolve
   failures within scope, and record actual outcomes and any remaining limits
   in this plan. Do not add package tests or start the deferred phases.

## Implementation authorisation

The final shared-understanding checkpoint is complete. The user instructed:
"We do, implement the plan". Remote commands still require individual approval.

## Remote validation plan

Ask the user before **each** remote command,
including copies and read-only inspections. Execute one command at a time as
`jlees`, never root, on
`ubu-m-14c7768bd2eb75f5.windows.ebi.ac.uk`, in
`/home/jlees/sketchlib.rust`.

- Build: `CUDA_HOME=/usr/local/cuda-13.4 cargo oxide build`.
- Baseline GPU smoke command:
  `target/release/sketchlib dist -o 1k_tmp -k 21 --knn 10 -v --gpu 0 1k_subset`.
- After agreeing the CLI spelling and tile semantics, repeat the smoke run
  with `--max-sketch-batch 256` to exercise diagonal tiles, rectangular tiles,
  and a short final block when the actual sample count is not divisible by
  256. Also cover an exactly divisible partition and a one-sample final block,
  choosing the limits/subsets after confirming the remote dataset size.
- Integration: run the same dataset and options on CPU with `--gpu 0` removed,
  writing to a separate output prefix. Compare neighbour identities, ordering,
  distances, self-pair exclusion, and the expected number of neighbours using
  the agreed comparison policy above. Cover both Jaccard and `--ani` output.
- Boundary coverage used a temporary 256-name subset of `1k_subset` and a
  one-sample final block obtained with `--max-sketch-batch 999`.
- Verify the build command produces the executable being tested before using
  `target/release/sketchlib` as build evidence.

## Implementation and validation record

- Implemented `DistTile::new()`, `pipeline()`, and `consume()`, with shared
  pair decoding, bounded tile ranges, stream-ordered operations, and both heap
  insertions per unordered pair.
- Wired positive `--max-sketch-batch`, default 10000, through GPU dispatch.
  Stored the aligned host allocation in `Arc<SketchVec>` and shared CPU/GPU
  heap insertion and Jaccard/ANI row finalisation.
- Generalised the kernel to two sketch inputs with independent strides.
  Interleaving includes all bit planes. kNN computes inverse ANI before the
  f32 cast; the existing dense caller retains raw ANI semantics.
- Applied the agreed compatibility edits and the separately approved D2H
  lifetime fix. All intermediate and output device buffers survive `.sync()`.
  Included the requested future dense-streaming note in `consume()`.
- Fixed E0277 in the inclusive-triangle decoder by converting the boolean
  square-root correction with `i64::from(...)` before adding it to `sqrt_floor`.
- Local formatting and `git diff --check` passed. No local compilation or tests
  ran, and no tests were added to the package. The one-off comparator lives
  locally at `/private/tmp/sketchlib-cuda-compare.py` and remotely at
  `/tmp/sketchlib-cuda-compare.py`.

### Remote build environment

All remote commands were requested individually and run sequentially as `jlees`.
Source copies retained backups. The remote manifest originally lacked the async
settings already present locally; the local manifest was copied with a backup.

cuda-oxide remains unpinned, as requested. Validation resolved revision
`9df9e7cef561e03a066230b3bb433b3d8285da81`. The compatible remote Rust channel is
`1.99.1`; its previous file was backed up as
`rust-toolchain.toml.before-cuda-1.99.1`. The user updated `cargo-oxide`, resolving
the compiler-flag environment failure. No explicit bootstrap override was used
for the successful build.

Successful build, from `/home/jlees/sketchlib.rust`:

```sh
CUDA_HOME=/usr/local/cuda-13.4 /home/jlees/.cargo/bin/cargo oxide build
```

It completed in 18.91s with the `release` profile. Tests used that freshly built
`target/release/sketchlib` executable on GPU 0.

### Runtime validation

The source dataset has 1000 new-format DNA sketches, sketch size 1024, and stored
k-mer list `[21]`. Every comparison used `-k 21 --knn 10`. ANI boundary runs used
256 names selected from the full CPU reference, passed with `--subset`.

| Mode | Loaded samples | Block limit | Coverage | Result |
| --- | ---: | ---: | --- | --- |
| Jaccard | 1000 | 10000 (default) | Single tile | Exact CPU match |
| Jaccard | 1000 | 256 | Diagonal and rectangular tiles; final block 232 | Exact CPU match |
| Jaccard | 1000 | 250 | Exact division into four blocks | Exact CPU match |
| Jaccard | 1000 | 999 | One-sample remainder; empty diagonal skipped | Exact CPU match |
| ANI | 256 | 10000 (default) | Minimum-size dataset; single tile | Exact CPU match |
| ANI | 256 | 128 | Diagonal/rectangular tiles; blocks below dataset minimum | Exact CPU match |
| Jaccard rejection | 100 | 10000 (default) | Dataset below minimum | Expected exit 1 and minimum-size error |

All six CPU/GPU comparisons passed the agreed structural and numeric rules:
matching sample and neighbour sets, ten distinct non-self neighbours per row,
finite values, correct sorting, maximum absolute distance error **0**, and
**0** cutoff-tie substitutions. The rejection message was
`GPU kNN requires at least 256 samples`.

CPU reference outputs are `1k_tmp_cpu` and `256_tmp_cpu_ani`; GPU outputs are
`1k_tmp`, `1k_tmp_tiles`, `1k_tmp_exact`, `1k_tmp_singleton`, `256_tmp_gpu_ani`, and
`256_tmp_gpu_ani_tiles`. They remain in the remote checkout. The temporary name
subset is `/tmp/sketchlib-cuda-256.txt`.

Local command logs include `/private/tmp/sketchlib-cuda-build-index-fix.log`,
`/private/tmp/sketchlib-cuda-compare-jaccard.log`,
`/private/tmp/sketchlib-cuda-compare-ani.log`, and
`/private/tmp/sketchlib-cuda-below-minimum.log`.

## References

- User-supplied source example:
  https://github.com/NVIDIA/cuda-rust/tree/main/cuda-oxide/crates/rustc-codegen-cuda/examples/async_mlp
- Documentation linked in the GPU source:
  https://nvidia.github.io/cuda-rust/cuda-oxide/async-programming/concurrent-execution.html
  (read during this interview; documents `.sync()` for blocking execution and
  `and_then` for ordering dependent stages on one stream).
- Stream scheduling documentation:
  https://nvidia.github.io/cuda-rust/cuda-oxide/async-programming/scheduling-and-streams.html
  (read during this interview; the same operation can be run synchronously or
  scheduled through a stream pool).

The current brief supersedes earlier-session plans for fixed-four-block async
execution. None of that earlier implementation is assumed to be present.
