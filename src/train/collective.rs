use std::{
    ffi::{CStr, c_char, c_void},
    ptr,
};

use anyhow::{Result, bail};
use tch::Tensor;

use super::distributed::Collective;

unsafe extern "C" {
    fn clair3_pg_init(
        output: *mut *mut c_void,
        rank: i32,
        world_size: i32,
        local_rank: i32,
        use_cuda: i32,
        error: *mut c_char,
        error_size: usize,
    ) -> i32;
    fn clair3_pg_allreduce(
        handle: *mut c_void,
        tensor: *mut c_void,
        world_size: i32,
        error: *mut c_char,
        error_size: usize,
    ) -> i32;
    fn clair3_pg_broadcast(
        handle: *mut c_void,
        value: i64,
        output: *mut i64,
        root: i32,
        local_rank: i32,
        use_cuda: i32,
        error: *mut c_char,
        error_size: usize,
    ) -> i32;
    fn clair3_pg_barrier(handle: *mut c_void, error: *mut c_char, error_size: usize) -> i32;
    fn clair3_pg_free(handle: *mut c_void);
}

pub struct C10dCollective {
    handle: *mut c_void,
    rank: usize,
    world_size: usize,
    local_rank: usize,
    use_cuda: bool,
}

fn check(status: i32, error: &[c_char]) -> Result<()> {
    if status == 0 {
        return Ok(());
    }
    let message = unsafe { CStr::from_ptr(error.as_ptr()) };
    bail!("C10d operation failed: {}", message.to_string_lossy())
}

impl C10dCollective {
    pub fn new(rank: usize, world_size: usize, local_rank: usize, use_cuda: bool) -> Result<Self> {
        if world_size == 0 || rank >= world_size || local_rank >= world_size {
            bail!("invalid distributed rank {rank}/{world_size} or local rank {local_rank}");
        }
        let mut handle = ptr::null_mut();
        let mut error = [0 as c_char; 4096];
        let status = unsafe {
            clair3_pg_init(
                &mut handle,
                rank as i32,
                world_size as i32,
                local_rank as i32,
                i32::from(use_cuda),
                error.as_mut_ptr(),
                error.len(),
            )
        };
        check(status, &error)?;
        Ok(Self {
            handle,
            rank,
            world_size,
            local_rank,
            use_cuda,
        })
    }
}

impl Collective for C10dCollective {
    fn rank(&self) -> usize {
        self.rank
    }

    fn world_size(&self) -> usize {
        self.world_size
    }

    fn broadcast_i64(&self, value: i64, root: usize) -> Result<i64> {
        if root >= self.world_size {
            bail!(
                "broadcast root {root} is outside world size {}",
                self.world_size
            );
        }
        let mut output = 0_i64;
        let mut error = [0 as c_char; 4096];
        let status = unsafe {
            clair3_pg_broadcast(
                self.handle,
                value,
                &mut output,
                root as i32,
                self.local_rank as i32,
                i32::from(self.use_cuda),
                error.as_mut_ptr(),
                error.len(),
            )
        };
        check(status, &error)?;
        Ok(output)
    }

    fn all_reduce_mean(&self, tensor: &Tensor) -> Result<Tensor> {
        let mut result = tensor.shallow_clone();
        let mut error = [0 as c_char; 4096];
        let status = unsafe {
            clair3_pg_allreduce(
                self.handle,
                result.as_mut_ptr() as *mut c_void,
                self.world_size as i32,
                error.as_mut_ptr(),
                error.len(),
            )
        };
        check(status, &error)?;
        Ok(result)
    }

    fn barrier(&self) -> Result<()> {
        let mut error = [0 as c_char; 4096];
        let status = unsafe { clair3_pg_barrier(self.handle, error.as_mut_ptr(), error.len()) };
        check(status, &error)
    }
}

impl Drop for C10dCollective {
    fn drop(&mut self) {
        unsafe { clair3_pg_free(self.handle) }
    }
}
