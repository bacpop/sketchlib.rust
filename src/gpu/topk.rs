//! Private CUB selection planning and stream-ordered dispatch.
use super::operations::allocate;
use anyhow::{ensure, Context, Result};
use cuda_async::{device_box::DeviceBox, device_operation::ExecutionContext, error::DeviceError};
use cuda_core::Context as CudaContext;
use std::{
    collections::HashMap,
    ffi::{c_char, c_int, c_void, CStr},
    mem, ptr,
};

extern "C" {
    fn sketchlib_tile_topk(
        workspace: *mut c_void,
        bytes: *mut usize,
        distances: *const f32,
        keys: *mut f32,
        ids: *mut u32,
        nr: usize,
        nq: usize,
        k: usize,
        diagonal: c_int,
        stream: *mut c_void,
    ) -> c_int;
    fn sketchlib_topk_error_string(status: c_int) -> *const c_char;
}

fn check(status: c_int) -> Result<(), DeviceError> {
    if status == 0 {
        return Ok(());
    }
    // CUDA runtime error strings are static, NUL-terminated strings.
    let message = unsafe { CStr::from_ptr(sketchlib_topk_error_string(status)) }.to_string_lossy();
    Err(DeviceError::Internal(format!(
        "CUB tile selection: {message} (runtime error {status})"
    )))
}

/// Prepared selection shapes for one device, reused during planning and dispatch.
#[derive(Default)]
pub(super) struct SelectionCache {
    shapes: HashMap<(usize, usize, usize, bool), Shape>,
}

impl SelectionCache {
    /// Query on a cache miss. Each cache must be used with only one device.
    pub fn prepare(
        &mut self,
        context: &CudaContext,
        nr: usize,
        nq: usize,
        k: usize,
        diagonal: bool,
    ) -> Result<Shape> {
        let key = (nr, nq, k, diagonal);
        if let Some(shape) = self.shapes.get(&key) {
            return Ok(*shape);
        }
        let shape = Shape::new(nr, nq, k, diagonal)?.query_on(context)?;
        self.shapes.insert(key, shape);
        Ok(shape)
    }

    /// Retrieve a shape prepared during memory planning; no CUDA calls occur here.
    pub fn get(&self, nr: usize, nq: usize, k: usize, diagonal: bool) -> Shape {
        self.shapes[&(nr, nq, k, diagonal)]
    }
}

#[derive(Clone, Copy)]
pub(super) struct Shape {
    nr: usize,
    nq: usize,
    k: usize,
    diagonal: bool,
    pub ref_k: usize,
    pub query_k: usize,
    pub count: usize,
    bytes: usize,
}

impl Shape {
    pub fn new(nr: usize, nq: usize, k: usize, diagonal: bool) -> Result<Self> {
        ensure!(
            nr > 0 && nq > 0 && k > 0,
            "CUB tile dimensions and k must be positive"
        );
        ensure!(
            nr <= u32::MAX as usize && nq <= u32::MAX as usize,
            "CUB block-local IDs exceed u32"
        );
        ensure!(!diagonal || nr == nq, "Diagonal CUB tile must be square");
        let ref_k = k.min(if diagonal { nr - 1 } else { nq });
        let query_k = if diagonal { 0 } else { k.min(nr) };
        let count = nr
            .checked_mul(ref_k)
            .and_then(|n| nq.checked_mul(query_k).and_then(|q| n.checked_add(q)))
            .context("CUB candidate count overflow")?;
        ensure!(
            count.checked_mul(mem::size_of::<f32>()).is_some(),
            "CUB candidate byte size overflow"
        );
        Ok(Self {
            nr,
            nq,
            k,
            diagonal,
            ref_k,
            query_k,
            count,
            bytes: 0,
        })
    }

    /// Device storage for candidate keys and IDs plus the queried CUB workspace.
    pub fn device_bytes(self) -> Result<usize> {
        self.count
            .checked_mul(mem::size_of::<f32>() + mem::size_of::<u32>())
            .and_then(|n| n.checked_add(self.bytes))
            .context("GPU selection byte size overflow")
    }

    pub fn query_on(mut self, context: &CudaContext) -> Result<Self> {
        context.bind_to_thread()?;
        // Query mode never dereferences these device pointers or enqueues work.
        check(unsafe {
            sketchlib_tile_topk(
                ptr::null_mut(),
                &mut self.bytes,
                ptr::null(),
                ptr::null_mut(),
                ptr::null_mut(),
                self.nr,
                self.nq,
                self.k,
                self.diagonal.into(),
                ptr::null_mut(),
            )
        })?;
        Ok(self)
    }
}

/// Enqueue selection on the caller's stream and return candidate storage.
/// Allocation leases retain candidates and workspace through submission completion.
///
/// # Safety
/// The execution context must be bound to the calling thread, and the distance
/// storage must already be retained for reading through submission completion.
/// The shape must have been queried on this device for these distances.
pub(super) unsafe fn select(
    ctx: &ExecutionContext,
    distances: &DeviceBox<[f32]>,
    mut shape: Shape,
) -> Result<(DeviceBox<[f32]>, DeviceBox<[u32]>), DeviceError> {
    let device_keys = unsafe { allocate::<f32>(ctx, shape.count)? };
    let device_ids = unsafe { allocate::<u32>(ctx, shape.count)? };
    let workspace = unsafe { allocate::<u8>(ctx, shape.bytes)? };
    check(unsafe {
        sketchlib_tile_topk(
            workspace.cu_deviceptr() as *mut c_void,
            &mut shape.bytes,
            distances.cu_deviceptr() as *const f32,
            device_keys.cu_deviceptr() as *mut f32,
            device_ids.cu_deviceptr() as *mut u32,
            shape.nr,
            shape.nq,
            shape.k,
            shape.diagonal.into(),
            ctx.get_cuda_stream().cu_stream().cast(),
        )
    })?;
    Ok((device_keys, device_ids))
}
