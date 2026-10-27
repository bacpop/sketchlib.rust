# CUDA self-distance planning

## Current round: finish GPU distance progress before sparse output

Status: implementation in progress. The user approved the minimal fix with
"Great, implement that minimal fix". Preserve the previous rounds and their
existing working-tree changes.

Reported problem and source findings:

- The user observes the GPU progress bar left at 99% when distance calculation
  ends and requests completion/closure before the sparse-writing log message.
- `self_dists_stream()` currently finalises and writes the sparse matrix in an
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

Implementation record (in progress):

- Moved progress completion/abandonment onto the distance-phase Result, before
  the existing sparse-finalisation closure. Successful completion now clears
  the bar before the writing log; sparse-output errors do not reopen that bar.
- The source change only relocates these lifecycle calls and adds an ordering
  comment. Local rustfmt and `git diff --check` passed. No local build or
  runtime tests ran. Remote build and one existing GPU kNN smoke run await
  individual approval; no new automated progress/quiet tests will be added.

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
