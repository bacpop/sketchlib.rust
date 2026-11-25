//! Distance-specific device state, tile memory admission and scheduling.
//! Kernel metadata and selection planning here belong to GPU self-distances;
//! this module is not general device infrastructure.
use super::{kernels, operations::h2d, topk};
use anyhow::{ensure, Context, Result};
use cuda_async::{
    device_box::DeviceBox,
    device_context::{global_policy, with_cuda_context},
    device_operation::DeviceOperation,
    scheduling_policies::SchedulingPolicy,
};
use cuda_core::{sys, Context as CudaContext, IntoResult};
use std::{
    collections::{HashSet, VecDeque},
    iter::Peekable,
    sync::Arc,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum TileKind {
    FullDiagonal,
    FullOffDiagonal,
    TailDiagonal,
    EdgeOffDiagonal,
}

impl TileKind {
    const COUNT: usize = 4;

    pub fn new(diagonal: bool, queries: usize, block: usize) -> Self {
        // In the upper-triangular traversal, only the query block can be
        // shorter than `block` on an off-diagonal tile.
        match (diagonal, queries < block) {
            (true, false) => Self::FullDiagonal,
            (false, false) => Self::FullOffDiagonal,
            (true, true) => Self::TailDiagonal,
            (false, true) => Self::EdgeOffDiagonal,
        }
    }

    pub fn is_diagonal(self) -> bool {
        matches!(self, Self::FullDiagonal | Self::TailDiagonal)
    }
}

#[derive(Clone, Copy)]
struct Geometry {
    refs: usize,
    queries: usize,
    kind: TileKind,
}
impl Geometry {
    fn distances(self) -> Result<usize> {
        let count = if self.kind.is_diagonal() {
            self.refs.checked_mul(self.refs - 1).map(|v| v / 2)
        } else {
            self.refs.checked_mul(self.queries)
        };
        count.context("GPU distance count overflow")
    }
}

pub(super) fn validate(ids: &[usize]) -> Result<()> {
    ensure!(!ids.is_empty(), "Select at least one GPU");
    let count = CudaContext::device_count()? as usize;
    let mut seen = HashSet::new();
    for &id in ids {
        ensure!(seen.insert(id), "GPU {id} was selected more than once");
        ensure!(
            id < count,
            "Invalid visible CUDA device {id}; {count} device(s) available"
        );
        with_cuda_context(id, |ctx| -> Result<()> {
            ctx.bind_to_thread()?;
            let mut unified = 0;
            unsafe {
                sys::cuDeviceGetAttribute(
                    &mut unified,
                    sys::CUdevice_attribute_enum_CU_DEVICE_ATTRIBUTE_UNIFIED_ADDRESSING,
                    ctx.cu_device(),
                )
            }
            .result()?;
            ensure!(
                unified != 0,
                "GPU {id} must support unified addressing for shared pinned sketches"
            );
            Ok(())
        })??;
    }
    Ok(())
}

/// Per-device distance kernels, sketch metadata and counters.
pub(super) struct DistanceDevice {
    pub id: usize,
    pub context: Arc<CudaContext>,
    pub policy: Arc<dyn SchedulingPolicy>,
    pub module: kernels::LoadedModule,
    pub kmers: Arc<DeviceBox<[usize]>>,
    pub completeness: Arc<DeviceBox<[f32]>>,
    pub tiles: usize,
    pub pairs: usize,
}
impl DistanceDevice {
    pub fn new(id: usize, kmers: &[usize], completeness: &[f32]) -> Result<Self> {
        let context = with_cuda_context(id, Arc::clone)?;
        context.bind_to_thread()?;
        let policy = global_policy(id)?;
        // SAFETY: this package owns the contracted kernel bundle.
        let module = unsafe { kernels::load_async(id) }.with_context(|| {
            format!("GPU {id}: cannot load distance kernels; rebuild for sm_75 to include Turing")
        })?;
        let kmers = Arc::new(h2d(kmers).sync_on(&policy.next_stream()?)?);
        let completeness = Arc::new(h2d(completeness).sync_on(&policy.next_stream()?)?);
        context.bind_to_thread()?;
        let mut uuid = std::mem::MaybeUninit::<sys::CUuuid>::uninit();
        let uuid: sys::CUuuid = unsafe {
            (
                sys::cuDeviceGetUuid_v2(uuid.as_mut_ptr(), context.cu_device()),
                uuid,
            )
                .result()?
        };
        let hex: String = uuid
            .bytes
            .iter()
            .map(|&b| format!("{:02x}", b as u8))
            .collect();
        let arch = unsafe { cuda_core::get_device_sm_name(context.cu_device())? };
        let (free, total) = memory_info(&context)?;
        log::info!(
            "CUDA device {id}: {}, UUID {}-{}-{}-{}-{}, {arch}, {free}/{total} bytes free",
            context.name()?,
            &hex[..8],
            &hex[8..12],
            &hex[12..16],
            &hex[16..20],
            &hex[20..]
        );
        Ok(Self {
            id,
            context,
            policy,
            module,
            kmers,
            completeness,
            tiles: 0,
            pairs: 0,
        })
    }
}

fn memory_info(context: &CudaContext) -> Result<(usize, usize)> {
    context.bind_to_thread()?;
    let (mut free, mut total) = (0, 0);
    unsafe { sys::cuMemGetInfo_v2(&mut free, &mut total) }.result()?;
    Ok((free, total))
}

fn geometries(n: usize, block: usize) -> Vec<Geometry> {
    let full = n / block;
    let tail = n % block;
    let mut out = Vec::new();
    let mut push = |refs, queries, diagonal| {
        out.push(Geometry {
            refs,
            queries,
            kind: TileKind::new(diagonal, queries, block),
        });
    };
    if full > 0 && block > 1 {
        push(block, block, true);
    }
    if full > 1 {
        push(block, block, false);
    }
    if tail > 1 {
        push(tail, tail, true);
    }
    if full > 0 && tail > 0 {
        push(block, tail, false);
    }
    out
}

fn footprint(geometry: Geometry, sketch_bytes: usize, shape: Option<topk::Shape>) -> Result<usize> {
    let sketches = geometry
        .refs
        .checked_add(geometry.queries)
        .and_then(|n| n.checked_mul(sketch_bytes))
        .and_then(|n| n.checked_mul(2))
        .context("GPU sketch buffer size overflow")?;
    let distances = geometry
        .distances()?
        .checked_mul(4)
        .context("GPU distance byte size overflow")?;
    let selection = match shape {
        Some(s) => s.device_bytes()?,
        None => 0,
    };
    sketches
        .checked_add(distances)
        .and_then(|n| n.checked_add(selection))
        .context("GPU tile byte size overflow")
}

pub(super) struct MemoryPlan {
    pub block: usize,
    pub limits: Vec<usize>,
    pub max_distances: usize,
    pub max_candidates: Option<usize>,
    pub nonempty_tiles: usize,
    selections: Option<Vec<[Option<topk::PreparedSelection>; TileKind::COUNT]>>,
}
impl MemoryPlan {
    pub fn selection(
        &self,
        device_index: usize,
        kind: TileKind,
    ) -> Option<topk::PreparedSelection> {
        self.selections.as_ref().map(|by_device| {
            by_device[device_index][kind as usize].expect("missing top-k plan for a nonempty tile")
        })
    }

    pub fn new(
        n: usize,
        requested_block: usize,
        requested_slots: usize,
        sketch_words: usize,
        knn: Option<usize>,
        devices: &[DistanceDevice],
    ) -> Result<Self> {
        ensure!(
            n >= 2 && requested_block > 0 && requested_slots > 0 && !devices.is_empty(),
            "Invalid GPU memory plan"
        );
        let sketch_bytes = sketch_words
            .checked_mul(8)
            .context("GPU sketch byte size overflow")?;
        let multi = devices.len() > 1;
        let budgets: Vec<usize> = devices
            .iter()
            .map(|d| {
                if multi {
                    memory_info(&d.context).map(|(free, _)| ((free as u128 * 4) / 5) as usize)
                } else {
                    Ok(usize::MAX)
                }
            })
            .collect::<Result<_>>()?;
        let knn = knn.map(|k| (k, topk::Limits::compiled()));
        let mut block = n.min(requested_block);
        loop {
            let geometries = geometries(n, block);
            let mut footprints = Vec::new();
            let mut max_distances = 0;
            let mut max_candidates = 0;
            let mut selections = knn.map(|_| Vec::with_capacity(devices.len()));
            for device in devices {
                let mut worst = 0;
                let mut device_selections = [None; TileKind::COUNT];
                for &geometry in &geometries {
                    let shape = if let Some((k, limits)) = knn {
                        let selection = topk::Shape::new(
                            geometry.refs,
                            geometry.queries,
                            k,
                            geometry.kind.is_diagonal(),
                        )
                        .and_then(|shape| {
                            topk::PreparedSelection::new(shape, limits, &device.context)
                        })
                        .with_context(|| {
                            format!("GPU {}: CUB workspace query failed", device.id)
                        })?;
                        device_selections[geometry.kind as usize] = Some(selection);
                        max_candidates = max_candidates.max(selection.shape.count);
                        Some(selection.shape)
                    } else {
                        None
                    };
                    max_distances = max_distances.max(geometry.distances()?);
                    worst = worst.max(footprint(geometry, sketch_bytes, shape)?);
                }
                footprints.push(worst);
                if let Some(selections) = &mut selections {
                    selections.push(device_selections);
                }
            }
            let decision = if multi {
                admission(
                    block,
                    &budgets,
                    &footprints,
                    requested_slots,
                    &devices.iter().map(|d| d.id).collect::<Vec<_>>(),
                )?
            } else {
                Admission::Fits(vec![requested_slots])
            };
            let caps = match decision {
                Admission::Fits(caps) => caps,
                Admission::Smaller(next) => {
                    block = next;
                    continue;
                }
            };
            let blocks = n.div_ceil(block);
            let tiles = blocks
                .checked_add(1)
                .and_then(|b| blocks.checked_mul(b))
                .map(|v| v / 2)
                .context("GPU tile count overflow")?;
            let empty = if block == 1 {
                blocks
            } else {
                usize::from(n % block == 1)
            };
            let nonempty_tiles = tiles - empty;
            let limits = distribute(&caps, nonempty_tiles);
            for (i, device) in devices.iter().enumerate() {
                let budget = if multi {
                    format!("{} byte memory budget", budgets[i])
                } else {
                    "automatic memory sizing disabled".to_owned()
                };
                log::info!("GPU {}: block size {block}, estimated {} bytes/tile, {} tile slots (requested {requested_slots}; {budget})",
                    device.id, footprints[i], limits[i]);
            }
            return Ok(Self {
                block,
                limits,
                max_distances,
                max_candidates: knn.map(|_| max_candidates),
                nonempty_tiles,
                selections,
            });
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
enum Admission {
    Fits(Vec<usize>),
    Smaller(usize),
}
fn admission(
    block: usize,
    budgets: &[usize],
    footprints: &[usize],
    requested: usize,
    ids: &[usize],
) -> Result<Admission> {
    let caps: Vec<_> = budgets
        .iter()
        .zip(footprints)
        .map(|(&budget, &bytes)| requested.min(budget / bytes.max(1)))
        .collect();
    if let Some(i) = caps.iter().position(|&cap| cap == 0) {
        ensure!(
            block > 1,
            "GPU {}: one minimum-size tile needs {} bytes, above its {} byte memory budget",
            ids[i],
            footprints[i],
            budgets[i]
        );
        Ok(Admission::Smaller((block / 2).max(1)))
    } else {
        Ok(Admission::Fits(caps))
    }
}
fn distribute(caps: &[usize], tiles: usize) -> Vec<usize> {
    let mut limits = vec![0; caps.len()];
    let mut remaining = tiles.min(caps.iter().fold(0usize, |n, &cap| n.saturating_add(cap)));
    while remaining > 0 {
        for (limit, &cap) in limits.iter_mut().zip(caps) {
            if remaining > 0 && *limit < cap {
                *limit += 1;
                remaining -= 1;
            }
        }
    }
    limits
}

/// A returned slot stays with its device. No fixed partition or GPU speed weights.
pub(super) struct TileQueue<I: Iterator, B> {
    tiles: Peekable<I>,
    ready: VecDeque<(usize, B)>,
}
impl<I: Iterator, B> TileQueue<I, B> {
    pub fn new(tiles: I, ready: VecDeque<(usize, B)>) -> Self {
        Self {
            tiles: tiles.peekable(),
            ready,
        }
    }
    pub fn take(&mut self) -> Option<(I::Item, usize, B)> {
        if self.ready.is_empty() || self.tiles.peek().is_none() {
            return None;
        }
        let (device, buffer) = self.ready.pop_front().unwrap();
        Some((self.tiles.next().unwrap(), device, buffer))
    }
    pub fn release(&mut self, device: usize, buffer: B) {
        self.ready.push_back((device, buffer));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn memory_limits_and_block_reduction() {
        assert_eq!(
            admission(64, &[1000, 400], &[100, 200], 4, &[2, 7]).unwrap(),
            Admission::Fits(vec![4, 2])
        );
        assert_eq!(
            admission(64, &[1000, 80], &[100, 200], 4, &[2, 7]).unwrap(),
            Admission::Smaller(32)
        );
        let mut block = 64;
        loop {
            match admission(
                block,
                &[1000, 80],
                &[block * block, block * block],
                4,
                &[2, 7],
            )
            .unwrap()
            {
                Admission::Smaller(next) => block = next,
                Admission::Fits(caps) => {
                    assert_eq!(block, 8);
                    assert_eq!(caps, vec![4, 1]);
                    break;
                }
            }
        }
        let error = admission(1, &[1000, 0], &[100, 200], 4, &[2, 7])
            .unwrap_err()
            .to_string();
        assert!(error.contains("GPU 7") && error.contains("200 bytes"));
    }

    #[test]
    fn slots_are_bounded_by_tiles_and_device_capacity() {
        assert_eq!(distribute(&[4, 2, 1], 20), [4, 2, 1]);
        assert_eq!(distribute(&[4, 2, 1], 4), [2, 1, 1]);
        assert_eq!(distribute(&[4, 4, 4], 1), [1, 0, 0]);
    }

    #[test]
    fn geometries_cover_remainders_and_ignore_empty_diagonals() {
        use TileKind::*;

        let planned_tiles = |n, block| {
            geometries(n, block)
                .into_iter()
                .map(|g| (g.kind, g.refs, g.queries))
                .collect::<Vec<_>>()
        };
        assert_eq!(
            planned_tiles(261, 64),
            [
                (FullDiagonal, 64, 64),
                (FullOffDiagonal, 64, 64),
                (TailDiagonal, 5, 5),
                (EdgeOffDiagonal, 64, 5),
            ]
        );
        assert_eq!(
            planned_tiles(256, 64),
            [(FullDiagonal, 64, 64), (FullOffDiagonal, 64, 64)]
        );
        assert_eq!(
            planned_tiles(257, 64),
            [
                (FullDiagonal, 64, 64),
                (FullOffDiagonal, 64, 64),
                (EdgeOffDiagonal, 64, 1),
            ]
        );
        assert_eq!(planned_tiles(256, 1), [(FullOffDiagonal, 1, 1)]);
        assert_eq!(planned_tiles(256, 256), [(FullDiagonal, 256, 256)]);
    }

    #[test]
    fn tile_memory_accounts_for_both_sketch_copies_and_selection() {
        let g = Geometry {
            refs: 7,
            queries: 5,
            kind: TileKind::EdgeOffDiagonal,
        };
        assert_eq!(footprint(g, 16, None).unwrap(), 2 * 12 * 16 + 35 * 4);
        let shape = topk::Shape::new(7, 5, 3, false).unwrap();
        assert_eq!(
            footprint(g, 16, Some(shape)).unwrap(),
            2 * 12 * 16 + 35 * 4 + 36 * 8
        );
        assert!(footprint(g, usize::MAX, None).is_err());
        assert!(Geometry {
            refs: usize::MAX,
            queries: 2,
            kind: TileKind::EdgeOffDiagonal,
        }
        .distances()
        .is_err());
    }

    #[test]
    fn ready_queue_balances_fast_and_slow_devices_without_duplicates() {
        let ready = VecDeque::from([(0, 10), (1, 20)]);
        let mut queue = TileQueue::new(0..80, ready);
        let mut pending = Vec::new();
        let mut seen = HashSet::new();
        let mut counts = [0; 2];
        let mut time = 0;
        loop {
            while let Some((tile, device, slot)) = queue.take() {
                assert_eq!(slot, if device == 0 { 10 } else { 20 });
                assert!(seen.insert(tile));
                pending.push((time + if device == 0 { 1 } else { 4 }, device, slot));
            }
            if pending.is_empty() {
                break;
            }
            assert!(pending.len() <= 2);
            pending.sort_unstable_by(|a, b| b.cmp(a));
            let (done, device, slot) = pending.pop().unwrap();
            time = done;
            counts[device] += 1;
            queue.release(device, slot);
        }
        assert_eq!(seen.len(), 80);
        assert_eq!(counts.iter().sum::<usize>(), 80);
        assert!(counts[0] > 2 * counts[1]);
        // Exhaustion leaves free buffers retained until queue teardown.
        assert_eq!(queue.ready.len(), 2);
    }
}
