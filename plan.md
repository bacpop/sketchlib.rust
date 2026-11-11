# CUDA self-distance planning

## Current round: multi-device GPU self-distances

Status: implementation, correctness validation, single-device performance
checks, mixed-device benchmarks and three-device profiling complete (2026-11-11).

### Accepted design and implementation

- One current-thread coordinator drives independent device-local pipelines,
  using each GPU's own default four-stream policy and explicitly scheduled
  futures. A shared tile iterator and device-tagged ready slots give faster
  devices more work without fixed weights, locks or extra GPU worker threads.
- Device state owns its module, context, policy, metadata and CUB shape cache.
  Raw operations bind their captured context and return allocation/transfer
  errors rather than panic. CUB scratch queries explicitly bind their target;
  dispatch binds again even after foreign-context pinned-buffer cleanup.
- Shared pinned sketches are read once; all selected devices must support UVA.
  Tile kNN selection stays on the GPU; one heap per sample and one dense writer
  remain shared on the host. Slots return only after host consumption/writer
  acknowledgement, preserving bounded buffering and cancellation ownership.
- The CLI accepts `--gpu 0,1` and repeated `--gpu` flags. Duplicate/invalid CUDA
  ordinals fail before loading sketch data. The original public single-device
  entry point wraps the device-list entry point. Existing CPU routing is retained.
- The per-device in-flight upper bound defaults to 4. Only multi-device calls
  automatically size memory: budget 80% of free VRAM after setup, halve the
  common block size if one tile cannot fit, then cap each device's slots using
  actual sketch/distance/candidate/workspace requirements. Aggregate slots do
  not exceed nonempty tiles. Single-device block/admission settings are unchanged.
- `.cargo/cuda-oxide.toml` sets the default Rust target to SM75; explicit oxide
  target overrides remain supported. Add SM75 native/PTX CUB coverage while
  retaining SM80 native/PTX and SM120 native images. The source package includes
  the new configuration and host scheduling module. No new C++ kernels or CUB
  algorithm changes were introduced.
- On device/writer failure, stop new work, abort/join jobs and retain submitted
  resources until completion. Errors identify the device and tile ranges;
  there is no silent skipping/retry. Logs report UUIDs, architecture, memory,
  chosen limits, per-device completed tiles/pairs and aggregate host merge time.

### Correctness validation completed

- Local CUDA 13.5 builds passed. The release suite passed 81 tests plus one
  doctest, with the hardware-dependent context test ignored by default.
  Explicitly running that ignored test on all three local GPUs also passed:
  it changes the current context through foreign pinned-buffer cleanup, queries
  CUB on the intended device, and interleaves scheduled CUB operations/data copies
  across the 5090, 3090 and 2080 Ti.
- The existing Turing CLI validator passed its smoke run, 33 CPU/GPU comparisons
  and two writer failures. Report: `/tmp/sketchlib-cuda-async-ms13oa_7/report.json`.
- Twelve additional dense/kNN CPU comparisons passed for mixed/reversed device
  lists, comma/repeated CLI forms, single/two/three GPUs, in-flight limits 1/4/8,
  both k-mers, Jaccard/ANI, completeness correction, reordered subsets, k=1/50/256,
  singleton remainders, a single tile and an 8,384-tile small-block stress case.
  Maximum error remained below 1e-6. Seven invalid-argument cases and both
  three-GPU writer failures (one/four formatter threads) passed.
- Real memory planning with 100k samples and a requested 100k block completed
  after selecting a 25k common block and 4/3/3 slots, covering all 4,999,950,000
  pairs across three devices. The single-device oversized case instead preserved
  its requested block and returned a clean GPU 2/tile allocation error requesting
  19,999,800,000 bytes; no panic or final sparse output occurred.
- Package listing includes `.cargo/cuda-oxide.toml`, the host scheduling module
  and the CUB source. Source/format checks pass. Artifacts and logs are in
  `/tmp/sketchlib-multigpu-implementation`. Earlier remote-testing instructions
  are superseded by the user's local CUDA 13.5 testing instructions.

### Performance acceptance

- Rebuilt original SM80 and architecture-matched SM75 baselines with the current
  lockfile/toolchain. Compared the final executable using five alternating
  measured runs after warm-ups on each GPU, rejecting samples with another
  process on the measured GPU or observed thermal throttling. All three devices
  passed the 3% single-device median regression gate.
- Completed uncontended single-device sets (five valid measured runs per
  executable, after warm-ups; 100k samples, k=50, ANI/k-mer 31, block 5000,
  four slots, 16 threads, default reactor, stdout `/dev/null`):

  | GPU | Original SM80 baseline | Matched SM75 baseline | New SM75 implementation |
  | --- | ---: | ---: | ---: |
  | RTX 5090 / CUDA 0 | 5.689 s | 5.486 s | 5.531 s |
  | RTX 3090 / CUDA 1 | 6.944 s | 6.941 s | 6.891 s |
  | RTX 2080 Ti / CUDA 2 | Unsupported target | 9.451 s | 9.453 s |

  The matched-baseline differences are +0.83%, -0.71% and +0.02%, respectively.
  One thermally throttled warm-up was rejected automatically. Monitoring excludes
  other compute processes on each measured GPU; workloads on other GPUs may still
  affect shared host resources. Raw measurements, rejected-sample diagnostics and
  final medians are in the performance subdirectory recorded by
  `/tmp/sketchlib-multigpu-implementation/latest-performance.txt`.
- After the user announced GPU availability, completed all four multi-device
  combinations with the same workload/settings, one warm-up and three valid
  measured runs per combination. One thermally throttled sample was rejected
  and replaced. Each run covered 210 tiles and all 4,999,950,000 pairs.

  | CUDA devices | GPUs | Median wall time |
  | --- | --- | ---: |
  | 0,1 | RTX 5090 + RTX 3090 | 4.485 s |
  | 0,2 | RTX 5090 + RTX 2080 Ti | 4.835 s |
  | 1,2 | RTX 3090 + RTX 2080 Ti | 5.338 s |
  | 0,1,2 | All three | 4.381 s |

  Compared with the new 5090-only median of 5.531 s, the 5090/3090 pair reduces
  total time by about 19%; all three reduce it by about 21%. The third device
  provides a small further gain on this workload. Raw measurements, monitoring
  and summaries are in the directory recorded by
  `/tmp/sketchlib-multigpu-implementation/latest-scaling.txt`.
- The two-device trace (`pair.nsys-rep`, `pair.sqlite`, `pair-summary.json` under
  the implementation artifact directory) confirms 210 distance tiles split
  113/97 across CUDA 0/1, four tile streams on each, and 1.238 s with distance
  kernels executing simultaneously on the two devices. All four explicit stream
  synchronisations precede the first tile. D2H volume remains 800 MB, and 3,852
  traced context switches total about 0.79 ms of API time. This is overlap/context
  validation, not a clean scaling comparison: another workload occupied CUDA 2
  during capture.
- The final three-device trace (`triple.nsys-rep`, `triple.sqlite`,
  `triple-summary.json` in the same artifact directory) splits the 210 tiles
  87/69/54 across CUDA 0/1/2, using four streams per device. Distance kernels
  overlap across at least two devices for 1.569 s and across all three for
  0.778 s. The GPU phase spans 1.991 s, with only 1.116 ms of aggregate time
  when no device runs a kernel or memory operation; the largest such gap is
  0.370 ms. All six explicit stream synchronisations precede the first distance
  kernel. There are no repeated tile-group synchronisation barriers.
- That trace records one shared 204.8 MB pinned sketch allocation, about 612 MB
  peak tracked device allocations per GPU, and 800 MB total D2H output. Its
  5,694 context-binding calls total about 1.10 ms of CPU API time. Host global
  heap merging takes 1.390 s in this profiled run. Profiling includes memory
  tracking, so benchmark wall times above come from unprofiled runs. No benchmark
  process remains paused; all requested validation is complete.

## Completed round: CUB tile-level kNN with host global selection

Status: implementation, local validation and paired profiling complete.

### Implementation and build

- A private CUDA C++ wrapper dispatches CUB through two C ABI functions: a
  workspace-query/execution entry point and runtime error-string lookup. Cargo's
  GPU-feature-gated build script uses `cc`/`nvcc` and the static CUDA runtime.
  The wrapper never allocates/frees GPU memory or synchronizes streams.
- GPU builds explicitly select `CUDA_TOOLKIT_PATH` (first) or `CUDA_HOME` and
  require CUDA 13.5+ with CCCL/CUB 3.5+. The build script reports missing/older
  toolkits clearly. The wrapper and build script are included in source packages.
- Native SM80 and SM120 code plus compute_80 PTX support the local Ampere and
  Blackwell devices. The extra native SM120 image avoids requiring newer-PTX JIT
  support from a driver operating under CUDA minor-version compatibility.
- CUB `DeviceBatchedTopK::MinPairs` reads logical reference rows, query columns
  and self-excluding condensed diagonal rows through device iterators, without
  another full matrix allocation. It uses float distances and u32 block-local IDs.
- Effective k is clamped per row. For k <= 2,048, 4,096-item chunks are reduced
  repeatedly on the GPU until each tile row has at most k candidates. Larger
  effective k uses `DeviceTopK::MinPairs` per row with reusable scratch. Output
  is unsorted and may choose different valid cutoff ties, matching the existing
  validation contract. Existing distance/inverse-ANI ranking is preserved.
- Selection and compact candidate D2H copies run on the tile's stream after
  `dists()`. Raw CUDA accesses retain scratch, device results, pinned transfer
  buffers and source owners in submission leases through completion/cancellation.
  Pinned outputs remain reusable per-slot allocations.
- The GPU kNN path no longer scans full matrices or builds tile-local CPU heaps.
  It merges candidates directly into the preallocated global host heaps, logs
  aggregate host merge time at verbose level, and preserves final sorting/output.
  Global merging stays on the dispatch thread; tile independence, stream policy,
  admission limits and dense-distance behaviour are preserved. Zero GPU k is
  rejected explicitly; public signatures and CLI arguments are unchanged.

### Validation

- `CUDA_HOME=/usr/local/cuda-13.5 cargo oxide build` passed, followed by the full
  existing suite via `cargo oxide test -- --release`: 76 tests and one doctest
  passed. No Cargo harness tests were added.
- The existing CLI runner passed its smoke run, 33 CPU/GPU comparisons and two
  writer-failure cases. These cover dense/kNN, ANI, completeness correction,
  reordered subsets, multiple layouts and worker counts. Maximum distance error
  remained approximately 1e-7 within the existing 1e-6 tolerance. Report:
  `/tmp/sketchlib-cuda-async-lx932h47/report.json`.
- Four additional end-to-end comparisons passed for both stored k-mers, k=1/50/256,
  reversed subsets, singleton remainders and in-flight limits 1/4/8. Both local
  CUDA ordinals 0 (RTX 5090) and 1 (RTX 3090) were exercised. Zero k rejection passed.
- Temporary ABI checks passed for rectangular rows/columns, condensed diagonals,
  ties, singleton endpoints, segment lengths 4,096/4,097, multi-row recursive
  reduction (3 x 20,000, k=2,048), and the large-k fallback (3 x 5,003, k=2,500).
  Probe: `/tmp/sketchlib-cub-abi.py`; packaged source listing contains all new files.
- Build-script checks passed for missing toolkit selection and precedence/version
  rejection using a synthetic older-version header. Manual ownership/diff review,
  Rust formatting and `git diff --check` passed.

### Paired performance profile

- The baseline was saved after global heap preallocation. Both runs use the
  RTX 3090 UUID `4b643a37-e171-2c2f-5e65-4abde7aca8c8` (CUDA ordinal 1), verified
  from trace metadata, with `100k_subset`, k-mer 31, 5,000-sample blocks, k=50,
  ANI, 16 CPU threads, four in-flight tiles and the default completion reactor.
- Logged application time: **10.051 s before, 6.606 s after**. Tile-phase duration
  (first `dists()` start to final D2H completion): **8.065 s to 4.946 s**.
- D2H volume fell from **19,999,800,000 bytes to 800,000,000 bytes**, approximately
  25-fold. Each off-diagonal tile transfers two 2 MB arrays instead of 100 MB of
  distances; each diagonal tile transfers two 1 MB arrays instead of 49.99 MB.
- Time with no traced kernels/copies/memsets during the tile phase decreased
  from **2.895 s to 0.000461 s**; the longest such gap fell from 85.615 ms to
  0.345 ms. The 800 CUB selection kernels totalled 0.370 s of kernel duration;
  `dists()` kernel durations totalled 4.556 s versus 4.279 s before. These are
  summed kernel durations, not isolated exclusive GPU times. Host global merging
  was 1.349 s and remains a potential future optimisation.
- Both traces contain only the two startup `cuStreamSynchronize` calls. CPU
  sampling remains unavailable under current user permissions, so profiles use
  CUDA/NVTX/OSRT tracing without sampling or CPU context-switch tracing. These
  measurements are one paired run, not a repeated benchmark or a hardware SM
  utilisation measurement.
- Baseline executable, logs, traces, exported SQLite, command measurements and
  metric summary are in `/tmp/sketchlib-cub-p02qvhml`.

## Completed round: preallocate global GPU kNN heaps

Status: implementation, correctness checks and paired profiling complete.

- Construct all global per-sample heaps individually using
  `(0..n).map(|_| BinaryHeap::with_capacity(knn + 1)).collect::<Vec<_>>()`.
  Repeating an empty heap with `vec![heap; n]` clones it without preserving its
  reserved capacity. The extra element accommodates the temporary insertion in
  `push_heap()` before the worst candidate is removed.
- This reserves storage before tile processing and removes global heap-growth
  allocations from serial merging. For n=100,000 and k=50, candidate capacity is
  about 81.6 MB on this 64-bit machine, excluding heap/allocator metadata.
- Tile-local heaps, ranking, scheduling, public APIs and dense output are unchanged.
  CUB tile selection remains separate future work.
- Local CUDA 13.5 builds and all 22 existing distance integration tests passed.
  Four CPU/GPU kNN comparisons (Jaccard/ANI, k=1/50) across multiple tiles passed
  with zero distance difference, using the existing cutoff-tie validation policy.
- CPU sampling is unavailable to this user: `nsys status --environment` reports
  `perf_event_open` denied, and `sudo -n` requires a password. Paired profiles use
  CUDA/NVTX/OSRT tracing with sampling and CPU context-switch tracing disabled.
  This permits GPU timeline and elapsed-time comparison, but not a new CPU
  allocation breakdown comparable to report4's privileged sampling.
- The pre-change executable and validation/profile artifacts are preserved in
  `/tmp/sketchlib-global-prealloc-4o9j43vn`. Both profile runs use CUDA device 1
  (RTX 3090),
  `100k_subset`, k-mer 31, block size 5,000, four tiles in flight, 16 CPU threads,
  k=50, ANI, stdout redirected to `/dev/null`, and the default completion reactor.
- GPU identity was verified from the CUDA-device mapping in report4 and both
  paired traces: all use UUID `4b643a37-e171-2c2f-5e65-4abde7aca8c8` (RTX 3090).
  The earlier RTX 5090 label was incorrect: CUDA ordinal 1 maps to GPU inventory
  ID / `nvidia-smi` index 0 here, while CUDA ordinal 0 maps to the RTX 5090.
  Compare physical UUIDs rather than assuming these numeric IDs are interchangeable.
- One paired profile measured logged application time of 10.269 s before versus
  9.563 s after. The tile phase (first `dists()` start to last D2H completion)
  decreased from 8.253 s to 7.484 s. Intervals without traced kernels, copies or
  memsets within that phase totalled 3.073 s before versus 2.372 s after; the
  longest such gap decreased from 92.209 ms to 63.921 ms.
- The CPU setup interval between the last initialization stream synchronization
  and the first output-buffer allocation increased from 2.561 ms to 40.381 ms.
  That interval includes global heap construction and sample-name preparation,
  not an isolated measurement of heap allocation alone.
- Profiles are `before.nsys-rep` and `after.nsys-rep` in the artifact directory,
  with exported SQLite traces, correctness logs, `measurements.json` (commands
  and profiler wall times), and `summary.json` (application/timeline metrics).
  Profiler wall times include collection/export overhead. The logged application
  interval was used above. These single-run results indicate improvement but
  are not a repeated benchmark or evidence that all host-side gaps are resolved.
- Manual review confirmed only the global heap constructor changed; final
  `git diff --check` passed. No Cargo harness tests were added.

## Completed round: configurable in-flight tile limit

Status: implementation and local CUDA validation complete.

- Add `dist --max-in-flight-tiles <N>`, defaulting to 4 and accepting positive
  integers. CPU paths accept the option without effect, like `--max-sketch-batch`.
- Append `max_in_flight_tiles: usize` to the public GPU `self_dists_stream()`
  arguments. Direct callers must supply the limit; zero returns an error before
  CUDA initialization. This supersedes the historical fixed-four-tile decision.
- Use the argument for the initial tile window. Allocate no more output slots
  than the number of nonempty tiles; retain slots through CPU consumption and
  preserve existing cancellation cleanup. Larger limits may increase pinned host
  and device memory use.
- Keep the separate CUDA stream pool at four streams with round-robin assignment.
  Increasing the tile limit does not create additional streams.
- Validate locally using CUDA 13.5: CLI defaults and invalid values, dense/kNN
  results for omitted/1/4/8 limits, window refills, fewer tiles than slots, writer
  failure at a nondefault limit, and direct API zero rejection.

Validation results:

- `CUDA_HOME=/usr/local/cuda-13.5 cargo oxide build` passed. Help displays the
  default of 4; zero and negative CLI values are rejected during parsing.
- `CUDA_HOME=/usr/local/cuda-13.5 cargo oxide test -- --release --test distance`
  passed all 22 existing tests. No Cargo harness tests were added.
- On the 257-sample two-k-mer fixture, ten dense/kNN GPU comparisons against CPU
  output passed: omitted/1/4/8 limits with block size 64, and limit 8 with block
  size 128 (fewer nonempty tiles than slots). Existing cutoff-tie rules were used.
- Writer failure with limit 8 exited correctly. An external API probe confirmed
  zero is rejected before CUDA initialization, even with an invalid device ID.
  Logs and the temporary probe are in `/tmp/sketchlib-tile-limit-gt4oq3ss`.
- Manual diff review and `git diff --check` passed. The stream policy and tile
  resource ownership mechanisms were not changed.

## Completed round: read directly into final host sketch storage

Status: implementation, local CUDA validation and independent reviews complete.

### Current validation policy

The user explicitly superseded the historical remote-only instructions below:
build and test locally on this machine using CUDA 13.5, for example
`CUDA_HOME=/usr/local/cuda-13.5 cargo oxide build`. Use
`CUDA_HOME=/usr/local/cuda-13.5 cargo oxide test -- --release` for the existing
test suite and run the CLI validation runner locally.
The earlier remote host, account and per-command approval requirements are
historical, not instructions for this round. No new Cargo harness tests are
required for this refactor.

### Agreed implementation

- Select one final host allocation before loading: aligned `SketchVec` for CPU
  use, or `Arc<PinnedHostBuffer<u64>>` for GPU self-distances. A private storage
  enum preserves the existing borrowed slice accessors and metadata file format.
- Keep existing CPU loading APIs. Add fallible `read_sketch_data_pinned(prefix,
  device_id)` and `read_sketch_data_block_pinned(prefix, names, device_id)` APIs.
  Both allocate and fill pinned memory before sharing it. The CLI selects these
  automatically for GPU self-distances, including subsets.
- Bulk-read full `.skd` files directly into initialized destination slices;
  subset loading copies sample byte ranges from the mmap into the final buffer.
  Decode little-endian words in place on big-endian hosts. Validate full-file
  length and subset ranges; propagate failures from the new pinned APIs.
- The public GPU distance API requires pinned-loaded sketches and returns an
  explanatory error for ordinary storage. There is no fallback host copy.
  Existing CPU APIs and CLI arguments remain available as before.
- GPU tiles clone the existing pinned owner and retain it through completion
  and cancellation. This supersedes the previous full host-copy decision.
  Kernel code, transposition, device allocations and tile scheduling are unchanged.
- Keep initialized allocations and their initial zero-fill. Add no allocator
  framework, uninitialized-memory API, registration bindings or new dependency.

### Validation record

- Local release build passed with `CUDA_HOME=/usr/local/cuda-13.5 cargo oxide build`.
  Existing compiler/documentation warnings remain; no build errors occurred.
- Focused existing reader and sketch integration tests passed, followed by the
  full existing suite via `CUDA_HOME=/usr/local/cuda-13.5 cargo oxide test -- --release`:
  76 tests and one doctest passed. No Cargo harness tests were added.
- On local GPU 0 (RTX 3090), `scripts/validate_gpu_distances.py --dataset
  /tmp/sketchlib-pinned-fixture-mbbd2_1a/small --kmer 31 --gpu 0` passed its
  full-loading smoke run, 33 CPU/GPU comparisons and two writer-failure cases.
  This covers reordered subsets, dense/kNN output, ANI, completeness correction,
  multiple tile layouts, singleton remainder and concurrent CPU consumption.
  The largest measured distance difference was approximately 1e-7 (within the
  existing 1e-6 tolerance). Report: `/tmp/sketchlib-cuda-async-je9wbdol/report.json`.
- A temporary 257-sample fixture containing both k=21 and k=31 passed eight
  additional dense/kNN CPU/GPU comparisons across full and reordered subset
  loading. Distances matched exactly; cutoff ties used the existing validation
  policy. Truncated/oversized full files were rejected on CPU and GPU, and a
  truncated selected sample returned a range error.
- A temporary external Rust API probe passed full/subset content and name-order
  equality, failed-load state preservation, full reload after subsetting, and
  rejection of ordinary host storage by the GPU API. Temporary fixtures/probe
  remain under `/tmp`; no benchmark or speedup claim is made.
- Independent standards and spec reviews found no issues. Formatting checks for
  edited reader/storage files, Python syntax and `git diff --check` passed.

## Memory optimisation backlog

The host-loading changes above are complete. The remaining items below are
proposals for later work, prioritising small changes; their speed benefits have
not been benchmarked.

- [x] **Bulk-read sketch data.** Replace per-word reads/pushes with direct bulk
  reads into the final destination; copy subset samples as mmap byte ranges.
- [x] **Eliminate the full host-to-host pinned copy.** Read directly into pinned
  storage for GPU use, while retaining aligned ordinary storage for CPU use.
  Registering the existing allocation with CUDA was considered as an alternative;
  direct pinned loading was chosen instead.
- [ ] **Remove unnecessary device zero-fills.** The transpose and distance
  kernels overwrite every valid output element, so their destination buffers
  can be allocated without a preceding device memset. Preserve initialized host
  buffers and verify that every device element is written before it is read.
  This is the smallest suggested next change.
- [ ] **Reuse diagonal-tile sketch storage.** Upload and transpose a diagonal
  tile's sketch block once, then share that read-only device allocation between
  the reference and query inputs. This avoids preparing the same block twice.
- [ ] **Transfer only the selected k-mer.** The implemented GPU distance path
  uses one k-mer, but currently transfers and transposes all stored k-mers.
  Packing only the selected k-mer would reduce pinned storage, transfer volume
  and transpose work for databases containing multiple k-mers.
- [ ] **Cache uploaded/transposed sketch blocks.** Reuse blocks across distance
  tiles to avoid repeated transfers and transpositions. This needs more ownership
  and scheduling changes than the preceding proposals; keep any cache bounded.
- [ ] **Consider bounded pinned staging buffers.** Use a few reusable pinned
  buffers instead of pinning all loaded sketches when host memory is limiting.
  This reintroduces repeated host copies for reused blocks, so it is primarily
  a memory-footprint tradeoff whose speed impact needs measurement.
- [ ] **Measure and optimise dense host formatting/writing.** On an H200, the
  user reports GPU transfers and kernels finishing within about 4.1 s with
  `--max-sketch-batch 100000 --max-in-flight-tiles 1`, while total execution
  exceeds 30 s. Smaller tiles and more in-flight tiles distribute the CPU work
  across execution but still leave substantial accumulated host time. These
  are user observations; the split between formatting, channel overhead and
  output I/O has not yet been measured. Dense consumption formats every pair
  into 1000-record text chunks using Rayon, feeds a bounded channel to one
  writer, and retains each tile's slot until writing is acknowledged. With one
  slot, successive tiles cannot overlap GPU work with host consumption.
  Compare the normal output destination with `/dev/null` (which still performs
  all formatting) and time host consumption/writing before choosing changes.
  Start with small improvements to formatting, allocation and chunk/channel
  overhead; preserve bounded memory and the existing output contract. This
  host bottleneck predates CUB selection; dense mode does not run CUB or reduce
  the number of output records.

## Completed round: finish GPU distance progress before sparse output

Status: implementation, remote release build and GPU kNN smoke run complete.
The user approved the minimal fix with "Great, implement that minimal fix".
Previous rounds and their existing working-tree changes have been preserved.

Reported problem and source findings:

- The user observes the GPU progress bar left at 99% when distance calculation
  ends and requests completion/closure before the sparse-writing log message.
- Before this fix, `self_dists_stream()` finalised and wrote the sparse matrix in an
  `and_then` closure, emits "Writing out in sparse matrix form", and only then
  calls `progress.finish_and_clear()`. This completion is after the requested
  distance-phase boundary and allows the writing message to interleave with an
  unfinished progress display.
- The inspected indicatif 0.17.11 implementation of `finish_and_clear()` sets
  position to the configured length, marks the display finished/hidden, and
  forces a final draw. No additional tick, progress counter, or manual position
  adjustment is needed.
- The user's example shows a 99% bar immediately before the sparse-writing
  log line, whereas CPU progress disappears before output. Further source
  inspection confirms that our `simple_logger` writes directly to stderr via
  `eprintln!`, while indicatif clears its previously drawn lines relative to
  the current terminal cursor. The logger does not coordinate with indicatif.
  Inference: logging first can move the cursor away from the bar, so the later
  clear may leave its last rendered 99% line visible. The method can finish its
  internal state without removing that stale earlier line. Rate-limited drawing
  also means the last rendered percentage need not equal the final counter.
  This is a source-supported explanation of the reported output, not a terminal
  reproduction; automated progress/quiet checks remain excluded.

Approved minimal fix:

- Finish and clear the existing bar immediately after successful GPU dispatch,
  CPU consumption, and heap merging have completed, before sparse finalisation
  or its writing message. Preserve the existing clear-on-completion behaviour.
- On distance-phase failure, abandon the bar and return the error. A subsequent
  sparse-output failure still returns an error; the completed distance bar does
  not need to be reopened or abandoned again. Dense completion remains after
  its producer/writer pipeline has completed.
- Keep quiet handling, progress increments, kernels, canonical transfer helpers,
  thread budgets, tile ownership and numerical/output behaviour unchanged. No
  new wrapper or additional progress mechanism is needed.
- Validate source ordering and formatting locally without compilation/runtime
  tests. Any remote copy/build/smoke command requires its own approval and must
  run sequentially as `jlees`. Prior exclusions of automated progress/quiet
  checks and agent-run scaling benchmarks remain in force; do not add new tests
  for this lifecycle-only change.

Implementation detail verified in the installed dependency:
`indicatif-0.17.11/src/state.rs::finish_using_style()` handles `AndClear` by
setting the position to the length and calling `draw(true, now)`.

Implementation and validation record:

- Moved progress completion/abandonment onto the distance-phase Result, before
  the existing sparse-finalisation closure. Successful completion now clears
  the bar before the writing log; sparse-output errors do not reopen that bar.
- The source change only relocates these lifecycle calls and adds an ordering
  comment. Local rustfmt and `git diff --check` passed. No local build or
  runtime tests ran.
- With individual approval, copied the source and plan to the remote checkout,
  then ran the release build and one existing GPU kNN smoke command sequentially
  as `jlees`, with a non-root guard. The build succeeded in 20.30 seconds; the
  smoke run (`-k 21 --knn 10 -v --gpu 0 1k_subset`) reached "Complete".
- Remote smoke output:
  `/tmp/sketchlib-gpu-progress-smoke-1793116255615.tsv`.
  Local captured build/smoke log:
  `/private/tmp/sketchlib-gpu-progress-finish-validation.log`.
- No automated progress/quiet checks or scaling benchmarks ran. The captured
  noninteractive log does not verify rendered terminal clearing; source ordering
  now guarantees the finish call precedes the sparse-writing log message.

## Completed round: parallel CPU consumption within each distance tile

Status: implementation, remote release build, smoke test, 33 CPU/GPU correctness
comparisons and both writer-failure checks complete. The user authorised this
round with "Implement it". Existing working-tree changes have been preserved.
Scaling validation remains reserved for the user's manual testing.

Confirmed scope and constraints:

- Completed tiles must expose enough CPU work to use the provided `--threads`
  budget, including when only one tile is available. Remove the effective
  four-consumer/formatter limit imposed by one sequential CPU job per tile.
- Use the single existing shared Rayon pool across every tile. kNN has
  `--threads` workers; dense mode keeps its existing `--threads + 1` pool,
  including the one writer worker. Do not create per-tile pools, reserve extra
  workers, or multiply the worker count by the four concurrent tiles.
- Keep the four-tile admission/memory bound through consumption and the current
  pinned-buffer ownership/reuse. Preserve the single dense writer, 1000-record
  formatted chunks, bounded text queue, progress semantics, prompt failure
  behaviour, output format, unordered dense-pair contract, and cutoff-tie rule.
- Keep distance kernels and canonical H2D/D2H/zeros helpers outside this change.
  Avoid full distance copies, full-matrix expansion, redundant parameter checks,
  additional wrapper structures, and unrelated algorithm changes.
- Compilation and runtime checks remain remote-only, sequentially as `jlees`,
  with approval before every remote command. Reuse the repository's existing
  Python validation runner; do not use the Cargo test harness. Continue to
  exclude progress/quiet checks and f32/f64 cutoff precision-boundary tests.
- The user will validate scaling manually. Do not run scaling/performance
  benchmarks, measure CPU utilisation, or add worker-count sweeps for scaling.

Source findings at the start of this round:

- `DistTile::consume()` scanned a tile sequentially and updated heaps
  for both endpoints of each unordered pair. The four-tile admission bound thus
  permitted at most four CPU consumers, even with a larger configured Rayon pool.
- `n_formatters = min(MAX_IN_FLIGHT_TILES, pool_size - 1)` sized the dense text
  channel; changing that value alone would not parallelise consumption.
- kNN tile heaps are merged into global heaps serially on the dispatcher. The
  user decided to keep this merge unchanged, as recorded below.
- Existing compact-index arithmetic (`square_to_condensed` for diagonal tiles,
  row-major indices for rectangles) permits accessing a sample's neighbours
  without allocating a square matrix or recomputing distances.
- Rayon parallel iterators and nested tasks reuse their current pool. This
  allows all tiles to share the configured workers without creating four pools.

References:

- [Rayon shared pool execution](https://docs.rs/rayon/latest/rayon/struct.ThreadPool.html#method.install)
- [Parallel mutable slice partitions](https://docs.rs/rayon/latest/rayon/slice/trait.ParallelSliceMut.html#method.par_chunks_mut)

Agreed kNN work partitioning:

- Process sample rows (or groups of rows) in parallel, giving
  each row exclusive access to its existing tile-local heap. Each worker reads
  the stored distances for that row and updates only that heap. A pair's stored
  value is read once for each endpoint, but calculated on the GPU only once.
  This keeps heap memory independent of worker count and needs no per-distance
  locks or additional heap copies. Compact triangular/rectangular storage stays
  intact. The user accepted the extra read and selected parallel sample rows.
- At the user's request, a comment beside `DistTile::consume()` records this
  decision and a future option: benchmark locks or a suitable atomic heap design
  if work partitioning keeps collisions rare. Such synchronization may permit a
  single pass over the distance buffer; it is not part of this round.
  The comment now describes the implemented approach.
- Splitting pair ranges into private partial heaps was considered but not
  selected because heap storage and merging grow with parallel partitions.

Global heap merging:

- Keep the existing serial merge and unsorted iteration through the tile-local
  heap unchanged. The user considered a sorted-vector merge with early stopping,
  but explicitly rejected the extra complexity if sorting costs O(k log k).
- Rust's BinaryHeap source explicitly documents `into_sorted_vec()` as in-place
  heapsort with O(k log k) complexity. Ordinary iteration costs O(k) but has no
  sorted-order guarantee, so it cannot use the proposed prefix early stopping.
  The condition for leaving the merge unchanged is therefore satisfied. Do not
  add sorted conversion, a new merge loop, an empty-heap fast path, or parallel
  global updates in this round.

Dense formatting within a tile:

- Distribute the existing 1000-record chunks across parallel tasks in the same
  shared Rayon pool. Keep one writer worker and the existing formatted-text
  channel; size its bound at twice the formatter budget (`pool_size - 1`),
  removing the four-tile cap from that calculation. Tasks across all tiles share
  that budget; there is no new pool or per-tile worker reservation.
- Wait for a tile's formatting tasks to finish sending before sending its final
  completion marker. The writer acknowledges the marker after all that tile's
  chunks have been written, retaining the current buffer/slot release contract.
  Keep allocation bounded by active tasks and the channel instead of collecting
  all formatted chunks into a vector. Preserve prompt cancellation on failure.

Validation scope:

- The user explicitly reserved scaling validation for their own manual testing.
  No representative benchmark command or dataset is needed from them, and no
  scaling claims should be inferred from correctness tests.
- After implementation, reuse the existing GPU smoke run, 33 CPU/GPU numerical
  comparisons and two writer-failure checks. Keep their existing one/four-worker
  coverage and comparison/tie rules; do not add scaling sweeps or benchmarks.
  All execution remains remote-only and individually approved, as before.
- Inspect the shared-pool scheduling and removal of the four-formatter cap in
  source. Confirm that no new pool or writer worker is created and that all
  work within and across tiles uses the configured worker budget. The four-tile
  admission/buffer bound remains independent of CPU parallelism.

Completed implementation sequence:

1. Change kNN consumption to parallel sample rows, with exclusive access to
   each existing tile-local heap. Address the compact triangular/rectangular
   buffer directly; perform no distance recomputation, copies, or matrix
   expansion. Update the requested comment to describe the implemented approach
   and retain its future low-contention locks/atomics option.
2. Change dense consumption to parallel 1000-record chunks within the same
   shared pool, sending text to the existing writer. Decode each chunk's first
   pair and advance the indices within the chunk, matching CPU streaming's
   traversal. Wait for all sends before the tile completion marker and retain
   prompt error cancellation.
3. Size the text queue from the full formatter budget, without the four-tile
   cap. Ensure the single reserved writer is running before formatting tasks
   can block on the bounded queue. Do not create another pool or change the
   `--threads`/`--threads + 1` convention. Preserve buffer lifetime/reuse, tile
   admission, progress and the serial unsorted global heap merge.
4. Format and inspect the diff locally without compilation or runtime tests.
   Request approval for remote copies, build, and the existing validation
   runner. Record actual outcomes here; leave scaling validation to the user.

Implementation and validation record:

- kNN now uses parallel mutable iteration over the existing tile-local row
  heaps. Each row reads its neighbours directly from compact diagonal or
  rectangular storage; pair values are read once for each endpoint. The requested
  comment describes this choice and preserves the future locks/atomics option.
- Dense mode now formats parallel 1000-record chunks. Its final marker is sent
  only after those tasks join, retaining the existing writer acknowledgement and
  pinned-output lifetime. The text queue uses twice the full consumer budget,
  without a four-tile cap, and a startup acknowledgement reserves the existing
  writer worker before producers can block on the queue.
- Source inspection confirms no new pools/workers, distance copies, kernel or
  canonical-helper changes, global-merge changes, or validation-script changes.
  Local rustfmt and `git diff --check` passed; no local build or runtime tests
  have run.
- After the approved copy, the sequential remote build and unchanged validation
  runner both passed as `jlees`. The release build used
  `CUDA_HOME=/usr/local/cuda-13.4 /home/jlees/.cargo/bin/cargo oxide build` and
  completed in 20.56s. No build fixes, dependency pins or toolchain changes were
  needed in this round.
- The GPU kNN smoke run on `1k_subset` passed. All 33 CPU/GPU comparisons passed
  at absolute tolerance 1e-6, covering dense/kNN, Jaccard/ANI, absent/mixed
  completeness, one/two/four-block layouts, and the existing one/four-worker
  correctness cases. Dense outputs contained all 32640 unordered pairs, with
  no self-pairs or duplicates; kNN outputs had 256 rows with ten distinct
  neighbours in score order. Alternative neighbours satisfied the previously
  agreed CPU-verified cutoff-tie rule.

| Completeness | Metric | Dense maximum absolute error | kNN maximum absolute error |
| --- | --- | ---: | ---: |
| Absent | Jaccard | 0 | 0 |
| Absent | ANI | 0 | 0 |
| Mixed | Jaccard | 1.4e-7 | 1.3e-7 |
| Mixed | ANI | 1e-7 | 1e-7 |

- The singleton-remainder dense case also passed, with maximum error 1.4e-7.
  Mixed-completeness checks retained their existing above/below-cutoff coverage;
  no precision-boundary tests were added.
- Both existing `/dev/full` runs produced the expected writer error and exited
  without hanging, with `--threads 1` and 4. Progress/quiet checks remain
  excluded. No scaling benchmarks, CPU utilisation measurements, extra worker
  sweeps or performance claims were made; the user will validate scaling.
- The combined build/validation log is
  `/private/tmp/sketchlib-cuda-parallel-rows-validation.log` locally. Outputs,
  fixtures, per-command logs and `report.json` remain remotely in
  `/tmp/sketchlib-cuda-async-n38oj4e6/`. The unchanged reusable runner is
  `scripts/validate_gpu_distances.py`.
- Final source review and `git diff --check` passed. The four-tile admission
  bound is now independent of parallel tasks within those tiles, while all CPU
  row/chunk tasks share the existing consumer pool and dense mode's single
  reserved writer. The global heap merge remains serial and unsorted.

Merge reference:

- [Rust BinaryHeap sorted conversion](https://doc.rust-lang.org/std/collections/struct.BinaryHeap.html#method.into_sorted_vec)
- [Rust BinaryHeap source and heapsort complexity](https://doc.rust-lang.org/src/alloc/collections/binary_heap/mod.rs.html)

## Completed round: asynchronous tiles, pinned sketches, and parallel consumers

Status: implementation, remote release build, smoke test, 33 CPU/GPU comparisons,
and both writer-failure checks complete. The user confirmed shared understanding
and authorised implementation/testing. Every remote command was approved and
run sequentially as `jlees`; no local compilation or runtime tests ran.

Requested scope:

- Schedule independent `DistTile` pipelines concurrently through cuda-oxide's
  asynchronous operations. Keep the public `self_dists_stream()` entry point
  synchronous, using a Tokio runtime inside `run_gpu_batches()` to execute and
  join an async block. The existing tile pipeline already composes H2D,
  interleaving, distances, and D2H and retains its buffer owners.
- Use `PinnedHostBuffer` for host sketch transfers and add CPU parallelism for
  kNN consumption. Resolve ownership, heap merging, and bounded buffering before
  implementation. Preserve dense streaming through the existing buffered writer
  and its unordered-pair output contract.
- Add progress reporting that observes `quiet`, using the agreed host-side
  consumed-pair count. The user accepted pauses during GPU computation in
  exchange for avoiding device counters and polling.
- Preserve Jaccard/ANI and f32 completeness behaviour from the validated previous
  round, the user-controlled block limit (default 10000), and the total-dataset
  GPU minimum of 256. The original `self_dists_knn_generic()` name in the brief
  now corresponds to the unified GPU `self_dists_stream()` entry point.
- Preserve canonical helpers except the explicitly approved pinned-buffer
  `d2h()` change and the permitted generic H2D helper, retaining the original
  H2D implementation with an `_u64` suffix. Type-only variants remain permitted;
  unrelated changes to canonical helpers require agreement. Keep `zeros()` intact.
  Avoid redundant parameter checks and unnecessary wrapper abstractions.
- Prefiltering, shared-memory query caching, reference/query comparisons, and
  core/accessory regression remain outside this round.
- Compile and test only on the remote machine, sequentially as `jlees`, with
  explicit approval before every remote command, including reads and copies.
  Use `CUDA_HOME=/usr/local/cuda-13.4 cargo oxide build`. Keep a reusable Python
  smoke/CPU-comparison runner in this repository and document its invocation;
  do not use the Cargo test harness.

Source findings at the start of this round:

- `DistTile::pipeline()` already built one stream-ordered operation; its caller
  ran `.sync()` and consumed each result before constructing the next tile.
  `run_gpu_batches()` was a deferred stub.
- The documented default `StreamPoolRoundRobin` has four streams. It assigns
  operations to streams; the application must separately bound submitted tiles
  and completed results waiting for CPU consumption.
- Existing dense consumption wrote directly to the supplied writer, which is
  already buffered by `set_ostream()`. The prior redundant text batching was
  removed from that synchronous consumer; this round's requested parallel text
  producers give bounded formatted chunks a distinct inter-thread purpose.

Documentation:

- [Concurrent execution](https://nvidia.github.io/cuda-rust/cuda-oxide/async-programming/concurrent-execution.html)
- [Scheduling and streams](https://nvidia.github.io/cuda-rust/cuda-oxide/async-programming/scheduling-and-streams.html)
- [Overlapping transfers and compute](https://nvidia.github.io/cuda-rust/cuda-oxide/async-programming/overlapping-transfers-and-compute.html)
- [Async MLP pipeline](https://nvidia.github.io/cuda-rust/cuda-oxide/projects/async-mlp-pipeline.html)
- [Rayon in-place scope](https://docs.rs/rayon/latest/rayon/fn.in_place_scope.html)

Resolved decisions for this round:

- Admit at most four tiles across GPU execution and CPU consumption. A tile
  retains its slot until consumption finishes, so slow consumers also bound
  completed host results. Implement the limit as a constant, proposed name
  `MAX_IN_FLIGHT_TILES: usize = 4`, so it can be changed later without adding a
  command-line option in this phase. The user explicitly accepted this limit
  and requested a constant.

Agreed pinned-memory design and source findings:

- The inspected `PinnedHostBuffer::from_slice()` allocates new page-locked
  storage and copies the input; it does not pin the existing `Arc<SketchVec>`
  allocation in place. It exposes immutable slices and implements Send/Sync for
  u64, so an `Arc<PinnedHostBuffer<u64>>` can safely supply overlapping read-only
  ranges to multiple tile H2D transfers, provided ownership survives completion.
- The user approved copying all loaded sketch data once when creating pinned
  memory, then sharing ranges of that allocation among tiles. Implement this:
  create one GPU-call-owned pinned allocation from `get_sketch_all()` and share
  it with Arc. Retain the pinned owner through each operation's completion,
  including error cleanup, rather than relying only on a borrowed host slice.
  This replaces the unaccepted per-slot staging proposal; no repeated host
  staging copies or per-slot sketch buffers would be needed.
- `get_sketch_all_arc()` is used only by `DistTile::new()`. Once tiles share the
  pinned allocation instead, remove that accessor, the Arc storage wrappers in
  `MultiSketch`, and its Arc-specific empty-storage helper. Restore ordinary
  aligned `SketchVec` ownership there and use the existing `get_sketch_all()`
  slice accessor for the one-time copy. No CUDA-specific storage variant is
  needed in `MultiSketch`; CPU storage and on-disk layout remain independent of
  CUDA.
- Accepted memory tradeoff: the original aligned sketch allocation remains in
  `MultiSketch` while one full pinned duplicate exists for the GPU call. Thus
  pinned sketch memory scales with the loaded dataset, whereas the four-tile
  limit bounds device allocations and pending distance results. Add a comment
  at the full pinned allocation explaining that a reusable copy buffer might
  reduce this memory cost in future, as explicitly requested by the user.
  Allocation failures should propagate; do not silently switch strategies.

Agreed CPU consumer design:

- Use the existing Rayon pool (and existing `--threads`
  setting) for kNN tile consumers. Each consumer builds local top-k heaps only
  for the sample rows touched by its tile, updating both endpoints with the
  shared pair decoder. A single owner merges those retained candidates into
  the global heaps as consumers complete, then releases each tile's slot.
  Keeping each tile's top k per sample suffices for the global top k; tie
  identities retain the existing unspecified cutoff-tie behaviour.
- This avoids locking global heaps for individual distance updates and
  avoids allocating a full-dataset heap array for every consumer. Merging on
  completion also avoids retaining per-tile heaps until the end of the run.
  The user accepted this approach with an explicit constraint against excess
  thread creation.
- All consumers share that one Rayon pool. Submit one sequential heap-building
  job per tile; do not create per-tile pools or run nested parallel loops within
  tile consumers. At most `min(MAX_IN_FLIGHT_TILES, Rayon pool size)` such jobs
  execute simultaneously, and queued consumers still hold their tile slots.

Agreed runtime arrangement:

- Implement `run_gpu_batches()` using a Tokio current-thread
  runtime on the calling host thread, rather than adding a Tokio worker pool.
  GPU tile futures can overlap across the four CUDA streams while this single
  thread submits/polls them. Offload substantial CPU consumption to the shared
  Rayon pool so it does not block GPU dispatch. This also keeps cuda-oxide's
  thread-local device state and stream pool on one host thread. Release completed
  device-buffer owners on that thread before handing host values to consumers.
  The user accepted this arrangement and reaffirmed the helper's purpose:
  it accepts an async block and drives it to completion while the public
  distance function and its callers remain synchronous.

Agreed dense producer/writer topology and chunk settings:

- The user requested the CPU streaming pattern: multiple Rayon producers
  receive completed distance buffers and format text, while one consumer alone
  writes that text to the existing `BufWriter`. This replaces the unaccepted
  proposal to perform formatting and writing in a single consumer.
- Keep one sequential formatting producer per tile in the same shared Rayon
  pool, without nested parallelism or additional pools. Send bounded text chunks
  through a multi-producer/single-consumer channel to one writer job. The async
  GPU dispatcher must not perform bulk formatting, wait on a blocking channel,
  or perform output I/O.
- Source inspection confirms the dense CLI path already calls
  `set_threads_with_writer(threads)`, configuring one shared pool with
  `threads + 1` workers. Use that existing writer reservation: at most
  `min(MAX_IN_FLIGHT_TILES, threads)` formatting producers plus the one writer
  execute concurrently. With `--threads 1`, one formatter and one writer can
  still make progress. Do not add another reservation or create per-tile pools.
  The user explicitly confirmed the CPU convention of `--threads + 1` as the
  total dense-mode pool budget. kNN keeps the `--threads` pool budget. The
  current-thread Tokio dispatcher adds no Tokio worker pool.
- The four-tile admission limit includes formatting, queued chunks, and writing.
  Release a dense tile's slot only after all its chunks have been written, not
  just after formatting or enqueueing. Output order, including chunks from
  different tiles, is unspecified. Slow output must bound both raw results and
  queued formatted text; formatting an entire tile into one String is excluded.
- Scoped Rayon jobs can borrow the existing writer and sample names without
  cloning them or requiring a new public `'static` ownership contract. Use an
  in-place scope so the dispatcher stays on its original host thread.
- Use the CPU streaming chunk size of 1000 records
  and bound the text channel at twice the number of possible active formatting
  producers. These chunks transfer formatting work to the writer; retain the
  standard `BufWriter` without introducing another byte-buffering wrapper.
  The user accepted these chunk and channel settings.

Agreed progress reporting:

- Count fully consumed unordered sample pairs
  out of `n * (n - 1) / 2`, weighted by comparisons rather than tile count.
  Dense mode advances after each text chunk has been written; kNN advances
  after each tile's local heaps have been merged into the global heaps.
- Use the existing host progress-bar helper, with `quiet` hiding the bar.
  The single writer worker advances dense progress only after successful chunk
  writes; the host coordinator advances kNN progress after merging each tile's
  local heaps. There are no per-pair progress updates from heap consumers.
  No polling thread/task, device atomic, unified-memory allocation, or kernel
  progress argument is needed. The user explicitly accepted pauses during GPU
  computation and preferred this simpler design.

Approved D2H helper change:

- Before this round, canonical `d2h()` allocated an ordinary `Vec<f32>` destination.
  The inspected `memcpy_dtoh_async()` documentation requires pinned host memory
  for guaranteed asynchronous behaviour; CUDA documents that transfers involving
  pageable host memory may block the host thread even with an Async API.
- Pass a reusable pinned output allocation into `d2h()` and
  return `(PinnedHostBuffer<f32>, DeviceBox<[f32]>)`, retaining the current device
  owner lifetime. Copy `dev.len()` values into that destination; consumers use
  the tile's logical distance range, without copying back into a Vec.
  Allocate the bounded output buffers before concurrent tile execution and reuse
  them after consumption, avoiding pinned allocations in the submission loop.
- The user explicitly approved this canonical-helper change and required that
  it avoid extra host copies. Keep the stream-ordered memcpy behaviour and
  device-buffer retention; `zeros()` stays unchanged. This approval does not
  authorise unrelated changes to canonical helpers.

Agreed failure behaviour:

- The user rejected treating partial results as useful output. On the first
  failure, fail the whole run and return an Error immediately or as soon as
  graceful cleanup permits. A full rerun is required; do not attempt to complete,
  salvage, resume, merge, or keep emitting queued results after detecting failure.
- Stop admitting tiles, signal CPU consumers to stop promptly, and close/drop
  channels as appropriate to unblock producers. Discard queued distances and
  formatted text. Join CPU jobs as they exit to protect borrowed data, rather
  than letting them finish processing their entire remaining tiles. Preserve
  the original substantive error over secondary channel-closure errors.
- Aborting host futures does not cancel already-submitted CUDA work. The
  inspected cuda-oxide `DeviceFuture` implementation defers reclamation of its
  owned output after cancellation, retaining resources behind a completion event.
  Retain all pinned/device owners in the operation output so this mechanism can
  protect them. Prefer prompt future cancellation and SDK-managed safe cleanup;
  wait for GPU completion only if required for safe cleanup with the dependency
  actually compiled remotely. Do not keep normal result consumption running just
  to drain work, or force a blocking SDK reclamation drain without necessity.
- Stop/clear progress on failure without reporting successful completion and
  never write an incomplete sparse matrix. Any output emitted before an error
  is invalid and requires a complete rerun. No output-salvage or file-transaction
  feature is requested.

Thread configuration for callers:

- Reuse the already configured Rayon pool; the GPU helper must not create or
  resize another one. Keep the CLI's existing dense writer reservation and kNN
  thread count. Document the same setup requirement for direct library callers:
  dense mode needs producer capacity plus its writer worker. Scheduling must not
  deadlock if that precondition is missing; use a clear host-side configuration
  error if needed, rather than silently adding threads.

Agreed remote validation:

- Keep the reusable stdlib Python runner in
  `scripts/validate_gpu_distances.py` (not the Cargo test harness), adapting the
  existing one-off runner and documenting its command in this plan. Only run it
  remotely after explicit approval; its subprocesses execute sequentially.
- Retain the supplied 1000-sample GPU kNN smoke command and the previous CPU/GPU
  comparison policy: absolute tolerance 1e-6, exact unordered-pair coverage for
  dense results, kNN structure/order, and CPU-verified alternatives at cutoff
  ties. Cover Jaccard/ANI and absent/mixed completeness away from the excluded
  f32/f64 cutoff precision boundary.
- Retain single-tile and singleton-remainder coverage. Add forced four-block
  integration runs (256 samples, block limit 64, ten tiles) so the four-tile
  window must refill. Exercise both dense and kNN with `--threads 1` and with
  multiple workers, checking results against the same CPU references.
- Exercise a dense writer failure using `/dev/full`, with a subprocess timeout
  and a required nonzero exit, to check prompt error propagation and absence of
  channel/pool deadlocks. The user explicitly excluded progress/quiet checks:
  do not add terminal capture, progress assertions, or quiet-mode test runs.
  Progress must still observe `quiet` in the implementation.
- Numerical validation establishes correctness, not a guaranteed speedup or
  hardware-level overlap for kernels that saturate the GPU. Inspect dispatch
  ownership and the bounded scheduling path alongside runtime comparisons.

Completed implementation sequence:

1. Restore aligned `SketchVec` ownership in `MultiSketch`, remove its now-unused
   Arc accessor/setup, and create one `Arc<PinnedHostBuffer<u64>>` inside the GPU
   call from the existing slice accessor. Add the requested future-copy-buffer
   comment. Preserve the metadata/data file format. Consolidate typed H2D
   helpers using the permitted generic implementation and retain `h2d_u64()`.
2. Allocate at most four reusable pinned output buffers before launching tiles.
   Apply the approved `d2h()` change and adapt `DistTile::pipeline()` to retain
   the pinned sketch owner and all device owners through async completion and
   cancellation. Consumers read only the tile's logical output range, directly
   from pinned memory. Preserve existing pair decoding, kernels, completeness
   arithmetic, and empty-diagonal handling; keep `zeros()` unchanged.
3. Implement `run_gpu_batches()` as the synchronous wrapper around a Tokio
   current-thread runtime accepting an async block. Dispatch independent tile
   operations through the existing stream policy, retaining a constant four-tile
   admission bound through CPU consumption. Use completion notifications to
   refill slots as they become available; no per-tile `.sync()`. Add only the
   dependency features needed for these channels/tasks; do not pin cuda-oxide.
4. Use the single existing Rayon pool through an in-place scope. kNN consumers
   build local heaps for each tile's participating rows, then the coordinator
   merges them. Dense consumers format 1000-record chunks and feed the single
   writer through the agreed bounded channel. Return/reuse output buffers only
   when the tile's slot is released. Preserve the dense and sparse output formats
   and ordering contracts, and the agreed thread budgets.
5. Add host progress at successful writes/heap merges, respecting `quiet` without
   another thread. Wire the agreed prompt failure path through dispatch,
   consumers, and writer: stop further work, discard pending results, preserve
   the original error, and perform only required safe cleanup. Keep additional
   abstractions and host scheduling checks to what this ownership requires.
6. Add the reusable runner to `scripts/validate_gpu_distances.py` and document
   its invocation/results here. Format and inspect changes locally without
   building or running tests. Obtain approval separately for every remote copy,
   build, read, or validation invocation; run commands sequentially as `jlees`.
   Fix in-scope failures and record actual validation outcomes and limits.

Reusable validation invocation from the remote repository:

```sh
python3 scripts/validate_gpu_distances.py --repo /home/jlees/sketchlib.rust --dataset 1k_subset --gpu 0
```

The build remains a separate authorised command unless explicitly included in
an approved sequential build-and-validation invocation. The runner does not
authorise itself to connect remotely or run commands without user approval.

Implementation and validation record:

- The implementation uses two task sets on a current-thread runtime:
  GPU completions hand pinned outputs to scoped Rayon jobs; consumer completions
  return them only after local heap merging or the writer's final-chunk
  acknowledgement. At most four output buffers exist, bounding the entire
  tile lifetime. Pinned input/output owners survive async cancellation.
- Dense mode uses the existing shared pool with its writer reservation and a
  bounded formatted-text channel. kNN uses tile-local heaps without nested
  parallel loops. `MultiSketch`'s Arc-only accessor/setup have been removed.
- The reusable runner is now in `scripts/validate_gpu_distances.py`. Its matrix
  has 33 CPU/GPU comparisons, the GPU smoke run, and two writer-failure runs
  (`--threads 1` and 4). Progress/quiet checks are excluded as requested.
- Local rustfmt and `git diff --check` passed. No local compilation or runtime
  tests have run. The approved first remote build found a missing `Send` bound
  on the future behind `pipeline()`'s opaque return type. Added that interface
  bound; no runtime tests ran in the failed build-and-validation invocation.
  The initial log is `/private/tmp/sketchlib-cuda-async-validation.log` locally.
- The next approved remote release build passed in 20.13s using
  `CUDA_HOME=/usr/local/cuda-13.4 /home/jlees/.cargo/bin/cargo oxide build`.
  The resolved remote cuda-async source is version 0.5.0 at existing revision
  `9df9e7cef561e03a066230b3bb433b3d8285da81`. An approved source read verified
  that its `DeviceFuture::drop()` invokes deferred reclamation, matching the
  cancellation design. No cuda-oxide revision pin or toolchain change was added.
  Build/source-audit output is in
  `/private/tmp/sketchlib-cuda-async-rebuild.log` locally.
- The approved runner completed successfully against the new release executable
  on GPU 0. The `1k_subset` kNN smoke run passed. All 33 CPU/GPU comparisons
  passed at absolute tolerance 1e-6. The 256-sample integration subset spans the
  source collection and reverses sample order, checking logical indices as well
  as arithmetic. Mixed completeness cycles through 0.95, 0.85, 0.50 and 1.0.

| Completeness | Metric | Dense maximum absolute error | kNN maximum absolute error |
| --- | --- | ---: | ---: |
| Absent | Jaccard | 0 | 0 |
| Absent | ANI | 0 | 0 |
| Mixed | Jaccard | 1.4e-7 | 1.3e-7 |
| Mixed | ANI | 1e-7 | 1e-7 |

- Every table row passed for single-tile, two-block, and four-block layouts.
  Four-block runs produce ten tiles, refilling the four-tile window, and passed
  with both `--threads 1` and 4. Every dense output has exactly 32640 unordered
  pairs without duplicates or self-pairs. Every kNN output has 256 rows with
  ten distinct non-self neighbours, in score order. Changes in selected tied
  neighbours were verified against the CPU dense oracle under the agreed cutoff
  rule, including at different thread counts.
- The singleton-remainder dense run passed with maximum error 1.4e-7. Mixed
  completeness changes eligible pairs and leaves below-cutoff pairs unchanged
  within tolerance. No precision-boundary tests were added or run.
- Both `/dev/full` runs, with `--threads 1` and 4, produced the expected writer
  error and exited without hanging, each taking approximately 0.57s including
  process/GPU setup. No progress/quiet checks were added or run.
- The full successful runner log is
  `/private/tmp/sketchlib-cuda-async-tests.log` locally. Per-command logs,
  outputs, fixtures, and `report.json` remain remotely in
  `/tmp/sketchlib-cuda-async-xqn3gaxo/`. The runner remains in the repository for
  reuse and outside the Cargo test harness.
- Final source inspection and `git diff --check` passed. Runtime comparisons
  validate results and the tested failure path; no GPU overlap profiler or
  performance benchmark was run. The accepted full pinned duplicate remains a
  host-memory cost, with the requested copy-buffer comment at its allocation.

Additional transfer reference:

- [CUDA API synchronization behaviour](https://docs.nvidia.com/cuda/cuda-runtime-api/api-sync-behavior.html)

## Completed round: shared dense/kNN path and completeness

Status: implementation, remote release build, smoke test, and agreed integration
validation complete, including the follow-up removal of redundant dense text
batching and its remote rebuild and validation.
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
   buffered writer before moving to the next tile. kNN still adds each unordered
   pair to both heaps and uses the existing shared finalisation before writing.
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

Follow-up: remove redundant dense text batching

- At the user's request, removed `consume()`'s `String`, 64 KiB threshold,
  manual chunk writes, and formatting-trait import. Dense records now use
  `writeln!(writer, ...)` directly, retaining the existing writer-error context
  and relying on the caller's `std::io::BufWriter`.
- Local formatting and `git diff --check` passed. The approved remote release
  rebuild and existing smoke/integration runner both passed. All 17 CPU/GPU
  comparisons still passed, with maximum absolute error 1.4e-7 and the same
  CPU-verified cutoff-tie policy. No tests were added.
- The combined local build/validation log is
  `/private/tmp/sketchlib-cuda-direct-write-validation.log`. New runtime
  outputs, fixtures, per-run logs, and `report.json` remain remotely in
  `/tmp/sketchlib-cuda-phase2-y1uacegn/`.

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

## Research: multi-GPU runtime foundations

Research scope: runtime feasibility and proposed application integration, checked
2026-11-11 against the restored CUB tile-selection + host-global-heaps baseline.
No application implementation changes or multi-GPU distance runs were made;
limited driver/CUB probes are described below. `Cargo.lock` resolves cuda-oxide to
`a0a4dc78a20e06ec90827f597e950bdc49f723e2`; the inspected checkout is
`/home/marussel/.cargo/git/checkouts/cuda-oxide-6d394bb007f5e114/a0a4dc7`.
Earlier research against `8785190` does not establish the current API.

### Explicit scheduling supports one shared executor

- There is one process-wide primary `Context` per device ordinal; scheduling
  policies, default device, deallocator streams and optional pools are per OS
  thread. `init_device_contexts(default, num_devices)` sets the default and map
  capacity; it does **not** enumerate or initialize that many devices. It errors
  if the calling thread already initialized its map. Device contexts are created
  lazily. Sources: [device_context.rs:92–115,132–205,311–363](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-async/src/device_context.rs#L92).
- `global_policy(device_id)` returns the selected device's policy.
  `operation.schedule(policy.as_ref())?` is the explicit safe async path: it
  chooses that policy's stream and stores an `ExecutionContext` with the stream,
  device and pool in the returned `DeviceFuture`. Await this already-scheduled
  future; avoid using bare operation `.await`/`.sync()` to select devices inside
  concurrent tasks because those consult the calling thread's default device.
  Sources: [device_context.rs:357–363,388–405,479–490](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-async/src/device_context.rs#L357),
  [device_operation.rs:126–147,366–375,641–652](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-async/src/device_operation.rs#L126).
- **Dependency limitation:** selecting a device is not the same as making its
  CUDA context current. `with_cuda_context`/`with_device` only retrieves a context
  and invokes a closure. `ExecutionContext::new`, round-robin `next_stream`, and
  `DeviceFuture::execute` do not bind the captured context before executing user
  operations. Sources: [device_context.rs:374–383](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-async/src/device_context.rs#L374),
  [scheduling_policies.rs:60–66](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-async/src/scheduling_policies.rs#L60),
  [device_future.rs:238–249](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-async/src/device_future.rs#L238).
- Consequently one current-thread Tokio executor can drive all GPUs, provided
  raw CUDA boundaries explicitly call the intended `context.bind_to_thread()?`
  immediately before their synchronous host calls. There must be no `.await`
  between binding and those calls. This is an inference from the captured stream
  API and synchronous operation execution: another task cannot change the host
  thread's current context during that uninterrupted sequence. CUDA itself lists
  both single-host-thread and thread-per-GPU operation as supported approaches.
  A thread per GPU is an alternative if host submission becomes a measured
  bottleneck, not a prerequisite for CUDA overlap. Sources:
  [NVIDIA multi-GPU guide](https://docs.nvidia.com/cuda/cuda-programming-guide/03-advanced/multi-gpu-systems.html),
  [runtime.rs:289–298](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-core/src/runtime.rs#L289).

### Device-affine resources and the CUB boundary

- Load `kernels::load_async(id)` separately for each device; keep module,
  prepared launches, streams, metadata buffers, tile allocations and CUB scratch
  with that device. `DeviceBox` records its owning device from its stream.
  `PreparedLaunch` queries the function's own context and rejects a stream from a
  different context. Sharing a host-side owner across threads does not transfer
  its device allocation to another GPU. Sources:
  [module loader generation:237–257](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-oxide/crates/cuda-macros/src/cuda_module/mod.rs#L237),
  [device_box.rs:49–83,108–110](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-async/src/device_box.rs#L49),
  [launch.rs:899–910,1059–1071](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-core/src/launch.rs#L899).
- Prepared Rust kernel queries/launch wrappers bind their owning context through
  their typed CUDA methods. Raw CUB calls have no such automatic wrapper.
  The runtime uses the Driver API context current on the calling host thread;
  cuda-oxide retains each device's primary context, so explicitly binding the
  matching primary context suffices for runtime/driver interoperability. Sources:
  [module.rs:238–245](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-core/src/module.rs#L238),
  [launch_raw.rs:36–47](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-core/src/launch_raw.rs#L36),
  [runtime.rs:176–190](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-core/src/runtime.rs#L176),
  [NVIDIA runtime/driver interoperability](https://docs.nvidia.com/cuda/cuda-runtime-api/cuda_runtime_api/group__CUDART__DRIVER.html).
- Current `src/gpu/topk.rs:90–105` queries scratch with null pointers and a null
  stream; it carries no device identity. Query scratch under each intended
  device's binding, then use that device's result for its allocations/dispatch.
  Bind again in `Finish::execute` before CUB runtime dispatch and raw stream-only
  memory calls (`src/gpu/topk.rs:191–243`). Do not infer the binding from a prior
  allocation or kernel launch: unrelated resource cleanup can change it.
  This is an implementation requirement inferred from the source and NVIDIA's
  current-context rule, not an assertion that existing single-GPU code fails.

### Pinned host sketches can be shared with a qualification

- `PinnedHostBuffer` owns an allocation made by `cuMemAllocHost`; it exposes no
  portable-allocation flag or public raw-allocation constructor in this revision.
  It is `Send`/`Sync` when its element type permits. Its destructor binds the
  allocation context before freeing; it does not restore the previous binding.
  Sources: [pinned_host_buffer.rs:41–71,129–166](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-core/src/pinned_host_buffer.rs#L41),
  [memory.rs:209–214](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-core/src/memory.rs#L209).
- NVIDIA documents `cuMemAllocHost` allocations as accessible to every context on
  devices supporting unified addressing. Thus sharing the existing immutable
  pinned sketches is reasonable if every selected GPU's unified-addressing
  support is checked/required. Explicit `cuMemHostAlloc` with
  `CU_MEMHOSTALLOC_PORTABLE` guarantees pinned status across contexts; adding
  this as an explicit allocation option requires a dependency extension or a
  separate owning buffer type. Do not assume `PinnedHostBuffer` already requests
  that flag. Source: [NVIDIA memory API, cuMemAllocHost/cuMemHostAlloc](https://docs.nvidia.com/cuda/cuda-driver-api/cuda_driver_api/group__CUDA__MEM.html).
- Keep output buffers owned by their assigned device's pipeline and retain all
  input/output owners until transfer completion. Creating/freeing a pinned
  output through its own context is already handled by the current buffer type;
  CPU-side free is safe under that destructor's ownership contract but can
  change the freeing thread's CUDA binding. No per-GPU duplication of the entire
  host sketch allocation is inherently required by the UVA case above.

### Mixed generations and remaining validation

- `build.rs:53–75` requires Toolkit 13.5+ and CCCL/CUB 3.5+, and builds the CUB
  wrapper with native `sm_80`, native `sm_120`, and forward-compatible
  `compute_80` PTX. That makes the CUB code path potentially usable on mixed
  Ampere/Ada/Hopper/Blackwell devices with compute capability >= 8.0: same-major
  compatible native code is preferred; otherwise the driver can JIT PTX for the
  selected GPU. It cannot cover GPUs below CC 8.0. This is a code-compatibility
  inference, not hardware validation. Source:
  [NVIDIA binary/PTX compatibility](https://docs.nvidia.com/cuda/cuda-programming-guide/01-introduction/cuda-platform.html#binary-compatibility).
- These CUB build flags do not establish the artifact targets of the Rust
  distance kernels. Each selected device must also successfully load its own
  embedded Rust kernel bundle and satisfy its prepared launch contract. The
  checkout's loader distinguishes target-specific NVVM/LTOIR inputs and rejects
  incompatible targets; do not promise whole-application mixed-generation
  support based only on the CUB fatbin. Source:
  [cuda-host/ltoir.rs:104–143,246–255](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-oxide/crates/cuda-host/src/ltoir.rs#L104).
- Driver support must cover every chosen GPU and the emitted PTX/toolkit.
  NVIDIA lists >=580 as the CUDA 13.x minor-compatibility floor, with feature
  limits; a driver too old for the toolkit's PTX can still fail JIT. Native
  images reduce that dependence where available but do not remove all driver
  feature requirements. Source:
  [NVIDIA minor-version compatibility](https://docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html).

Recommended foundation: retain the shared current-thread executor and serial
host heap coordinator; explicitly schedule each complete tile operation through
its device policy, keep GPU resources per device, and bind at every raw CUDA/CUB
boundary. Before implementation is accepted, validate two actual devices
(including mixed generations if that is required), each device's scratch query,
UVA support, cancellation/cleanup, and bundle loading. Peer access, NCCL, and
GPU-side global selection are not prerequisites for independent tile dispatch.

### Limited local probes

- `/tmp/sketchlib-multigpu-research/driver_probe.py` passed eight H2D/D2H round
  trips using four non-blocking streams per device on the RTX 5090 (CUDA 0) and
  RTX 3090 (CUDA 1), alternately submitted by one host thread with explicit
  context binding. Both devices report UVA. All streams read one default
  `cuMemAllocHost` source allocation made under the first GPU's primary context;
  no portable-allocation flag or peer-access enablement was used. Results are in
  the adjacent `driver_probe.json`. This corroborates pinned sharing and basic
  stream/context setup, not Rust future scheduling or concurrent distance kernels.
- CUDA 13.5 `nvcc --list-gpu-arch` includes `compute_75`. Compiling the unchanged
  CUB wrapper into a temporary SM75 shared library succeeded. On the local RTX
  2080 Ti (CUDA 2), the existing ABI oracle passed rectangular (17 x 23, k=5),
  diagonal (257 x 257, k=50), and large-k fallback (3 x 5003, k=2500) cases.
  Probe/library: `/tmp/sketchlib-multigpu-research/turing_cub_probe.py` and
  `cub-sm75.so`. This is evidence that CUB/Toolkit 13.5 do not inherently rule
  out this older GPU; the current application fatbin still does not cover it.
- The locked Rust backend's `arch_satisfies_feature` allows basic kernels on
  SM70+, while cross-compilation defaults to SM80. Explicit SM75 is therefore
  plausible for these kernels, but no SM75 build of the complete application was
  tested. Add SM75 CUB coverage and verify a portable Rust target/bundle before
  claiming whole-application Turing support. Do not drop existing native/PTX
  coverage for newer devices just to add the older one. Source:
  [architecture validation:23–64](https://github.com/NVlabs/cuda-oxide/blob/a0a4dc78a20e06ec90827f597e950bdc49f723e2/cuda-oxide/crates/cuda-oxide-codegen/src/target/arch.rs#L23).

### Recommended application design

This is a proposed direction, not an implemented or benchmarked multi-GPU path.
The current application evidence is `src/gpu/mod.rs` (`DistTile::pipeline`,
`self_dists_stream`, host merge/writer), `src/gpu/topk.rs` (scratch query and
`Finish`), and `src/lib.rs`/`src/cli.rs` (single-device loading and selection).

1. **Keep one current-thread coordinator.** Add a small per-device state holding
   device ID/context, scheduling policy, loaded Rust module, device k-mer and
   completeness arrays, and output-slot pool. Initialise the context map once;
   obtain each selected GPU's policy with `global_policy(id)`. Each policy has
   its own four-stream round-robin pool by default. Stream selection within a
   GPU and load balancing between GPUs are separate decisions.
2. **Choose the GPU before constructing/submitting a tile operation.** Use its
   module/constants and `tile.pipeline(...).schedule(policy.as_ref())?`, then
   await the scheduled future. One complete H2D/interleave/dists/CUB/D2H chain
   stays on that GPU's selected stream. Tag completed tiles and CPU consumers
   with their device/slot ID. Do not build one mixed-device round-robin policy
   around operations containing another GPU's module or buffers.
3. **Make context binding explicit.** Bind the target before setup/prepare/CUB
   queries, at the beginning of execution of each custom tile operation, and
   immediately before raw CUB and stream-only memory calls. Never carry an
   ambient default-device/context assumption across an await. Query/cache CUB
   workspace per device and geometry; the same tile shape need not have the
   same scratch requirement on different architectures. This requires Rust
   host-side changes, not a new C++ algorithm or new distance kernels.
4. **Use ready slots as the dynamic work queue.** Keep one canonical tile
   iterator and a `VecDeque` of free `(device_id, output_buffer)` slots, initially
   interleaved across devices. Assign one next tile to each available slot.
   Return its slot to the same device after host heap merging or dense writer
   acknowledgement, exactly preserving buffer reuse/backpressure. Faster GPUs
   finish and reclaim slots sooner, so naturally receive more tiles. This fits
   the existing JoinSet/select coordinator without new GPU worker threads,
   blocking channels, locks around heaps, or per-device global heaps.
5. **Keep host result handling shared.** Retain one heap per sample globally,
   merging completed candidates serially as now, and one dense writer with the
   existing shared Rayon formatting pool. Read sketches once into shared pinned
   storage after checking every selected device's UVA support. Copy only the
   current tile's blocks to each assigned device; duplicate only small persistent
   device metadata. No NCCL, NVLink or peer access is needed for this design.

### Memory sizing and heterogeneous performance

- Recommend retaining one common block partition. First bound block size by
  `--max-sketch-batch`; if even one tile cannot fit a selected device, reduce the
  common block size and regenerate the grid before dispatch. Size against actual
  diagonal/rectangular/remainder geometries, not just the distance matrix.
  Different block grids per GPU cannot simply be enumerated independently:
  they overlap or leave gaps. Device-specific subdivision of canonical tiles is
  a later option if uniform tiles prove limiting, not needed for initial support.
- Estimate a tile's device footprint as `2 * (r+q) * S + 4*d + 8*c + W`, where S
  is bytes per stored sketch, d is the number of pair distances, c is the number
  of selected directional candidates, and W is that device's queried CUB
  workspace. Dense mode omits candidate arrays/workspace. The first term covers
  both raw and interleaved sketch buffers, including today's duplicate diagonal
  preparation. Include persistent metadata and runtime/allocator headroom.
- Choose a per-device admission cap from that footprint and currently free
  memory after setup, leaving headroom (80% of remaining free memory is a
  reasonable initial conservative budget, not an OOM guarantee). Keep the
  requested in-flight limit as an upper bound, not a requirement to fill every
  GPU with that many tiles. Log the common block size, estimated tile bytes,
  free-memory budget and effective slot count for each device. If one slot still
  cannot fit at the minimum block size, fail clearly before dispatch.
- Recommend defining `--max-in-flight-tiles` as a **per-device upper bound** for
  multi-GPU mode, default 4, with the single-GPU behaviour unchanged. Document
  the aggregate bound as the sum of device caps; pinned output memory grows with
  this sum. Do not allocate more aggregate slots than there are nonempty tiles,
  and preserve bounded host result/text queues. This semantics change for the
  new multi-GPU case must be explicit in help/docs before implementation.
- Pull scheduling balances observed completion throughput without GPU speed
  weights or a calibration benchmark. It also adapts to PCIe differences and
  varying device load. It is not optimal: a slow GPU can retain a few final
  in-flight tiles while a faster GPU becomes idle. Start with sufficiently small
  tiles and bounded queues; measure that tail before adding speed estimates,
  adaptive admission or migration. More GPUs are not guaranteed to be faster.
- The whole calculation does not have to reside on each GPU, and the device
  memories are not pooled: each needs only its own live tiles and metadata.
  Shared PCIe bandwidth, host RAM/NUMA placement, serial host kNN merging and
  dense formatting/writing can limit scaling. Dense output is already known to
  be host-bound in some runs; multi-GPU dispatch alone cannot remove that work.

### Interface, errors and validation before implementation acceptance

- Suggested CLI: accept comma-separated/repeated IDs, e.g. `--gpu 0,1` or
  `--gpu 0 --gpu 1`; preserve `--gpu 0` and CPU execution when omitted. Require
  unique valid visible CUDA ordinals and log names, UUIDs, compute capability
  and memory. Do not confuse CUDA ordinals with `nvidia-smi` enumeration. Explicit
  selection is preferable initially to an `all` option that silently includes an
  unsupported or unsuitable GPU.
- Preserve the existing public single-device Rust entry point as a wrapper;
  introduce a device-list entry point backed by the same coordinator. Keep the
  initial feature limited to the existing GPU self-distance modes (dense and
  tile-kNN/host-global-kNN). Do not accidentally route cross-query/core-accessory
  modes to unsupported GPU implementations.
- Preflight every selected device's module/bundle and CUB support before tile
  output. Mixed generations need compatible code for both the Rust kernels and
  the CUB wrapper, not matching GPUs. Initially fail clearly if a selected device
  is unsupported. A current SM80+ build can target the 3090/5090 pair; including
  the 2080 Ti requires the artifact changes/validation described above.
- On a device error, stop new submissions to all devices, abort/join outstanding
  futures/consumers and retain resources through actual GPU completion via the
  existing leases. Report the failing device and tile. Do not silently skip a
  device or retry already-consumed tiles: repeated candidates can corrupt the
  heap results, and dense output may already contain partial records.
- Validate exact pair coverage/no duplicate tile assignment; compare single-GPU
  and mixed-GPU dense/kNN output against CPU with existing tolerances and valid
  cutoff ties. Exercise device order reversal, one device, more devices/slots
  than tiles, singleton remainders, in-flight limits 1/4/8, memory caps that give
  different devices different slot counts, invalid/duplicate device IDs,
  unsupported targets and writer/device errors. Explicitly test interleaved
  CUB queries/dispatch with an intentionally wrong prior current context.
- Add scheduler-level tests with synthetic fast/slow completions: every tile
  must be assigned exactly once, slots remain device-owned/bounded, and the fast
  device receives more work without fixed weights. A small multi-GPU end-to-end
  run is still required: the driver probe does not establish application safety.
- Benchmark each GPU alone and selected combinations using the same dataset,
  tile parameters and output destination. Log per-device tile/sample-pair counts
  and execution time, trace both devices' overlapping work, and monitor other
  workloads/thermal state. Separate dense host-output time and global heap merge
  time so a host bottleneck is not mistaken for poor device scheduling. Use local
  CUDA 13.5 builds/tests; the old remote-testing instructions do not apply.
