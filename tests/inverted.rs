use sketchlib::hashing::{RollHash, HashType};
use sketchlib::multisketch::MultiSketch;
use sketchlib::sketch::Sketch;

const BBITS: u64 = 14;
const SIGN_MOD: u64 = (1 << 61) - 1;

struct TestHasher {
    hashes: Vec<u64>,
    current_index: usize,
    k: usize,
    seq_len: usize,
    initialized: bool,
}

impl TestHasher {
    fn new(hashes: Vec<u64>, seq_len: usize) -> Self {
        println!("Creating TestHasher with original hashes: {:?}", hashes);
        // Ensure hashes are valid and within range
        let scaled_hashes = hashes.into_iter()
            .map(|h| ((h.wrapping_mul(0x517cc1b727220a95)) << BBITS) % SIGN_MOD)
            .filter(|&h| h > 0)
            .collect::<Vec<_>>();
        println!("Scaled hashes: {:?}", scaled_hashes);
        
        Self {
            hashes: scaled_hashes,
            current_index: 0,
            k: 3,
            seq_len,
            initialized: false,
        }
    }
}

impl Iterator for TestHasher {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.initialized {
            self.initialized = true;
            return Some(self.hashes[0]);
        }

        if self.current_index < self.hashes.len() {
            let hash = self.hashes[self.current_index];
            println!("Yielding hash: {}", hash);
            self.current_index += 1;
            Some(hash)
        } else {
            println!("No more hashes to yield");
            None
        }
    }
}

impl RollHash for TestHasher {
    fn set_k(&mut self, k: usize) {
        println!("Setting k to {}", k);
        self.k = k;
    }

    fn curr_hash(&self) -> u64 {
        let hash = if self.current_index > 0 {
            self.hashes[self.current_index - 1]
        } else {
            self.hashes[0]
        };
        println!("Current hash: {}", hash);
        hash
    }

    fn hash_type(&self) -> HashType {
        HashType::DNA
    }

    fn seq_len(&self) -> usize {
        self.seq_len
    }

    fn sketch_data(&self) -> (bool, [usize; 4], usize) {
        (false, [1, 1, 1, 1], 0)
    }
}

#[test]
fn test_simple_inverted_index() {
    println!("Starting test with SIGN_MOD: {}", SIGN_MOD);
    
    // Use larger initial values to ensure they remain positive after scaling
    let mut hasher1 = TestHasher::new(vec![100, 200, 300], 4);
    let mut hasher2 = TestHasher::new(vec![100, 200, 400], 4);

    let kmer_lengths = vec![3];
    let sketch_size = 3;
    
    println!("\nCreating sketch1...");
    let mut sketch1 = Sketch::new(
        &mut hasher1,
        "genome1",
        &kmer_lengths,
        sketch_size,
        false,  // rc
        1,      // min_count
    );
    println!("Sketch1 created successfully");

    println!("\nCreating sketch2...");
    let mut sketch2 = Sketch::new(
        &mut hasher2,
        "genome2",
        &kmer_lengths,
        sketch_size,
        false,  // rc
        1,      // min_count
    );
    println!("Sketch2 created successfully");

    println!("\nCreating MultiSketch...");
    let mut sketches = vec![sketch1, sketch2];
    let multi_sketch = MultiSketch::new(
        &mut sketches,
        sketch_size,
        &kmer_lengths,
        HashType::DNA,
        false,
    );
    println!("MultiSketch created successfully");

    println!("\nCreating inverted index...");
    let inverted_index = MultiSketch::invert_index(&multi_sketch);
    println!("Inverted index created with {} entries", inverted_index.len());

    // Print the contents of the inverted index
    println!("\nInverted index contents:");
    for (hash, genomes) in &inverted_index {
        println!("Hash {}: present in genomes {:?}", hash, genomes);
    }

    assert_eq!(inverted_index.len(), 4, "Expected 4 unique hashes");

    // Get the actual hash values from our hasher instances for verification
    let hasher1_values: Vec<u64> = TestHasher::new(vec![100, 200, 300], 4).hashes;
    let hasher2_values: Vec<u64> = TestHasher::new(vec![100, 200, 400], 4).hashes;

    // Verify the presence of hashes in the inverted index
    let hash1 = hasher1_values[0];
    let hash2 = hasher1_values[1];
    let hash3 = hasher1_values[2];
    let hash4 = hasher2_values[2];

    assert!(inverted_index.get(&hash1).unwrap().contains(&0));
    assert!(inverted_index.get(&hash1).unwrap().contains(&1));
    assert!(inverted_index.get(&hash2).unwrap().contains(&0));
    assert!(inverted_index.get(&hash2).unwrap().contains(&1));
    assert!(inverted_index.get(&hash3).unwrap().contains(&0));
    assert!(!inverted_index.get(&hash3).unwrap().contains(&1));
    assert!(inverted_index.get(&hash4).unwrap().contains(&1));
    assert!(!inverted_index.get(&hash4).unwrap().contains(&0));
}