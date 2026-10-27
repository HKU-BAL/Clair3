use anyhow::Result;
use tch::Tensor;

pub trait Collective {
    fn rank(&self) -> usize;
    fn world_size(&self) -> usize;
    fn broadcast_i64(&self, value: i64, root: usize) -> Result<i64>;
    fn all_reduce_mean(&self, tensor: &Tensor) -> Result<Tensor>;
    fn barrier(&self) -> Result<()>;
}

pub struct LocalCollective;

impl Collective for LocalCollective {
    fn rank(&self) -> usize {
        0
    }

    fn world_size(&self) -> usize {
        1
    }

    fn broadcast_i64(&self, value: i64, _root: usize) -> Result<i64> {
        Ok(value)
    }

    fn all_reduce_mean(&self, tensor: &Tensor) -> Result<Tensor> {
        Ok(tensor.shallow_clone())
    }

    fn barrier(&self) -> Result<()> {
        Ok(())
    }
}
