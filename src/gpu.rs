
#[link(name = "vector_add", kind = "static")]
extern "C" {
    fn vectorAdd_main();
}

fn call_gpu() {
    unsafe {
        vectorAdd_main();
    }
}


sparse_coo query_db_sparse_cuda(std::vector<Reference> &ref_sketches,
                          const std::vector<size_t> &kmer_lengths,
                          RandomMC &random_match,
                          const int kNN, const size_t dist_col,
                          const int device_id,
                          const unsigned int num_cpu_threads) {
  dist_params params =
    cuda_dists_init(ref_sketches, ref_sketches, kmer_lengths, device_id);
  FlatRandom flat_random =
      random_match.flattened_random(kmer_lengths, ref_sketches[0].seq_length());
  std::vector<uint16_t> ref_random_idx =
      random_match.lookup_array(ref_sketches);

  // Flatten all of the sketches in blocks
  //   (as stride is bins) - return vector<vector<uint64_t>>
  // Use a fixed chunk size of up to 3000 samples per chunk
  const size_t chunks = 1 + params.n_samples / 3000;
  const size_t samples_per_chunk = params.n_samples / chunks;
  const unsigned int num_big_chunks = params.n_samples % chunks;
  std::vector<std::vector<uint64_t>> sample_blocks;
  std::vector<SketchStrides> sample_block_strides(chunks, params.ref_strides);
  size_t start_idx = 0;
  for (size_t chunk_idx = 0; chunk_idx < chunks; ++chunk_idx) {
    size_t end_idx = start_idx + samples_per_chunk +
      (chunk_idx < num_big_chunks ? 1 : 0);
    sample_blocks.push_back(
      flatten_by_samples(ref_sketches, kmer_lengths,
              sample_block_strides[chunk_idx], start_idx, end_idx, num_cpu_threads));
    start_idx = end_idx;
  }

  return sparseDists(params, sample_blocks, sample_block_strides,
                     flat_random, ref_random_idx, kmer_lengths,
                     kNN, dist_col, samples_per_chunk, num_big_chunks,
                     num_cpu_threads);
}

NumpyMatrix query_db_cuda(std::vector<Reference> &ref_sketches,
                          std::vector<Reference> &query_sketches,
                          const std::vector<size_t> &kmer_lengths,
                          RandomMC &random_match, const int device_id,
                          const unsigned int num_cpu_threads) {
  dist_params params =
    cuda_dists_init(ref_sketches, query_sketches, kmer_lengths, device_id);

  double est_size =
      (params.ref_strides.bbits * params.ref_strides.sketchsize64 * kmer_lengths.size() * params.n_samples *
           sizeof(uint64_t) + // Size of sketches
       kmer_lengths.size() * std::pow(random_match.n_clusters(), 2) *
           sizeof(float) + // Size of random matches
       params.n_samples * sizeof(uint16_t) +
       params.dist_rows * 2 * sizeof(float)); // Size of distance matrix
  std::cerr << "Estimated device memory required: " << std::fixed
            << std::setprecision(0) << est_size / (1048576) << "Mb"
            << std::endl;
  std::cerr << "Total device memory: " << std::fixed << std::setprecision(0)
            << params.mem_total / (1048576) << "Mb" << std::endl;
  std::cerr << "Free device memory: " << std::fixed << std::setprecision(0)
            << params.mem_free / (1048576) << "Mb" << std::endl;

  unsigned int chunks = 1;
  if (est_size > params.mem_free * (1 - mem_epsilon)) {
    if (params.self) {
      // To prevent memory being exceeded, total distance matrix is split up into
      // chunks which do fit in memory. The most is needed in the 'corners' where
      // two separate lots of sketches are loaded, hence the factor of two below

      // These are iterated over in the same order as a square distance matrix.
      // The i = j chunks are 'self', i < j can be skipped
      // as they contain only lower triangle values, i > j work as query vs ref
      chunks = floor((est_size * 2) / (params.mem_free * (1 - mem_epsilon))) + 1;
    } else {
      throw std::runtime_error(
          "Using greater than device memory is unsupported for query mode. "
          "Split your input into smaller chunks");
    }
  }

  // Turn the random matches into an array (same for any ref, query or subsample
  // thereof)
  bool missing = random_match.check_present(ref_sketches, true);
  if (missing) {
    std::cerr
        << "Some members of the reference database were not found "
           "in its random match chances. Consider refreshing with addRandom"
        << std::endl;
  }
  FlatRandom flat_random =
      random_match.flattened_random(kmer_lengths, ref_sketches[0].seq_length());
  std::vector<uint16_t> ref_random_idx =
      random_match.lookup_array(ref_sketches);

  // Ready to run dists on device
  SketchSlice sketch_subsample;
  std::vector<float> dist_results(params.dist_rows * 2);
  NumpyMatrix coreSquare, accessorySquare;
  if (params.self) {
    size_t calc_per_chunk = params.n_samples / chunks;
    unsigned int num_big_chunks = params.n_samples % chunks;
    // Only allocate these square matrices if they are needed
    if (chunks > 1) {
      coreSquare.resize(params.n_samples, params.n_samples);
      accessorySquare.resize(params.n_samples, params.n_samples);
    }
    unsigned int total_chunks = (chunks * (chunks + 1)) >> 1;
    unsigned int chunk_count = 0;

    sketch_subsample.query_offset = 0;
    for (unsigned int chunk_i = 0; chunk_i < chunks; chunk_i++) {
      sketch_subsample.query_size = calc_per_chunk;
      if (chunk_i < num_big_chunks) {
        sketch_subsample.query_size++;
      }

      sketch_subsample.ref_offset = sketch_subsample.query_offset;
      for (unsigned int chunk_j = chunk_i; chunk_j < chunks; chunk_j++) {
        if (total_chunks > 1) {
          std::cerr << "Running chunk " << ++chunk_count << " of "
                    << total_chunks << std::endl;
        }

        sketch_subsample.ref_size = calc_per_chunk;
        if (chunk_j < num_big_chunks) {
          sketch_subsample.ref_size++;
        }

        dist_results = dispatchDists(
            ref_sketches, ref_sketches, params.ref_strides, params.query_strides, flat_random,
            ref_random_idx, ref_random_idx, sketch_subsample, kmer_lengths,
            chunk_i == chunk_j, num_cpu_threads, params.shared_size);

        // Read intermediate dists out
        if (chunks > 1) {
          NumpyMatrix blockMat = Eigen::Map<
              Eigen::Matrix<float, Eigen::Dynamic, 2, Eigen::ColMajor>>(
              dist_results.data(), dist_results.size() / 2, 2);

          // Convert each long form column of Nx2 matrix into square distance
          // matrix Add this square matrix into the correct submatrix (block) of
          // the final square matrix
          longToSquareBlock(coreSquare, accessorySquare, sketch_subsample,
                            blockMat, chunk_i == chunk_j, num_cpu_threads);
        }
        sketch_subsample.ref_offset += sketch_subsample.ref_size;
      }
      sketch_subsample.query_offset += sketch_subsample.query_size;
    }
  } else {
    std::vector<uint16_t> query_random_idx =
        random_match.lookup_array(query_sketches);

    sketch_subsample.ref_size = ref_sketches.size();
    sketch_subsample.query_size = query_sketches.size();
    sketch_subsample.ref_offset = 0;
    sketch_subsample.query_offset = 0;

    dist_results = dispatchDists(
        ref_sketches, query_sketches, params.ref_strides, params.query_strides, flat_random,
        ref_random_idx, query_random_idx, sketch_subsample, kmer_lengths, false,
        num_cpu_threads, params.shared_size);
  }

  NumpyMatrix dists_ret_matrix;
  if (params.self && chunks > 1) {
    // Chunked computation yields square matrix, which needs to be converted
    // back to long form
    Eigen::VectorXf core_dists = square_to_long(coreSquare, num_cpu_threads);
    Eigen::VectorXf accessory_dists =
        square_to_long(accessorySquare, num_cpu_threads);

    dists_ret_matrix.resize(samples_to_rows(coreSquare.rows()), 2);
    dists_ret_matrix << core_dists, accessory_dists; // Join columns
  } else {
    dists_ret_matrix =
        Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, 2, Eigen::ColMajor>>(
            dist_results.data(), dist_results.size() / 2, 2);
  }

  return dists_ret_matrix;
}
