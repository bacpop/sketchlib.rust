//! Stream-ordered allocation, transfers and asynchronous resource ownership.
use cuda_async::{
    device_box::DeviceBox,
    device_context::with_default_device_policy,
    device_future::DeviceFuture,
    device_operation::{DeviceOperation, ExecutionContext},
    error::DeviceError,
};
use cuda_core::{
    device_buffer::DeviceCopy,
    memory::{malloc_async, memcpy_htod_async, memset_d8_async},
};
use std::{future::IntoFuture, marker::PhantomData, mem, sync::Arc};

/// Fallible custom operations avoid panics and bind their captured device at
/// execution time. The closure is monomorphised; borrowed H2D inputs stay borrowed.
struct RawOperation<O, F> {
    run: F,
    output: PhantomData<fn() -> O>,
}
fn raw_operation<O: Send, F: FnOnce(&ExecutionContext) -> Result<O, DeviceError> + Send>(
    run: F,
) -> RawOperation<O, F> {
    RawOperation {
        run,
        output: PhantomData,
    }
}
impl<O: Send, F: FnOnce(&ExecutionContext) -> Result<O, DeviceError> + Send> DeviceOperation
    for RawOperation<O, F>
{
    type Output = O;
    unsafe fn execute(self, ctx: &ExecutionContext) -> Result<O, DeviceError> {
        ctx.context().bind_to_thread()?;
        (self.run)(ctx)
    }
}
impl<O: Send, F: FnOnce(&ExecutionContext) -> Result<O, DeviceError> + Send> IntoFuture
    for RawOperation<O, F>
{
    type Output = Result<O, DeviceError>;
    type IntoFuture = DeviceFuture<O, Self>;
    fn into_future(self) -> Self::IntoFuture {
        match with_default_device_policy(|policy| self.schedule(policy.as_ref())) {
            Ok(Ok(future)) => future,
            Ok(Err(error)) | Err(error) => DeviceFuture::failed(error),
        }
    }
}

pub(super) fn h2d<T: DeviceCopy + Send + Sync + 'static>(
    host_data: &[T],
) -> impl DeviceOperation<Output = DeviceBox<[T]>> + '_ {
    raw_operation(move |ctx| {
        let buffer = unsafe { allocate::<T>(ctx, host_data.len())? };
        unsafe {
            memcpy_htod_async(
                buffer.cu_deviceptr(),
                host_data.as_ptr(),
                std::mem::size_of_val(host_data),
                ctx.get_cuda_stream().cu_stream(),
            )?;
        }
        Ok(buffer)
    })
}
pub(super) fn zeroed<T: DeviceCopy + Send>(
    n: usize,
) -> impl DeviceOperation<Output = DeviceBox<[T]>> {
    raw_operation(move |ctx| {
        let buffer = unsafe { allocate::<T>(ctx, n)? };
        unsafe {
            memset_d8_async(
                buffer.cu_deviceptr(),
                0,
                n * mem::size_of::<T>(),
                ctx.get_cuda_stream().cu_stream(),
            )?;
        }
        Ok(buffer)
    })
}

/// Register an owner before enqueueing raw asynchronous accesses to it.
pub(super) fn retain<T: Send + 'static>(owner: T) -> impl DeviceOperation<Output = ()> {
    raw_operation(move |ctx| ctx.retain(owner))
}

/// Allocate on the operation's stream and retain storage through completion.
///
/// # Safety
/// The execution context must be bound to the calling thread before allocation.
pub(super) unsafe fn allocate<T: DeviceCopy + Send>(
    ctx: &ExecutionContext,
    n: usize,
) -> Result<DeviceBox<[T]>, DeviceError> {
    let bytes = n
        .checked_mul(mem::size_of::<T>())
        .ok_or_else(|| DeviceError::Internal("Device allocation overflow".into()))?;
    let stream = ctx.get_cuda_stream();
    let ptr = unsafe { malloc_async(stream.cu_stream(), bytes) }.map_err(|error| {
        DeviceError::Internal(format!("Allocating {bytes} device bytes: {error}"))
    })?;
    let buffer = unsafe { DeviceBox::from_raw_parts(ptr, n, Arc::clone(stream)) };
    ctx.retain_buffer(buffer.storage().clone(), true)?;
    Ok(buffer)
}
