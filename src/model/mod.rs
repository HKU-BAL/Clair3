mod bloom;
mod checkpoint;
mod full_alignment;
mod layers;
mod pickle;
mod pileup;

use std::{fmt, path::Path};

use anyhow::{Result, bail};
use clap::ValueEnum;
use tch::{Device, Tensor, nn};

pub use bloom::BloomModel;
pub use checkpoint::load_checkpoint;
pub use full_alignment::FullAlignmentModel;
pub use pickle::load_python_state_dict;
pub use pileup::PileupModel;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum ModelKind {
    #[value(name = "pileup")]
    Pileup,
    #[value(name = "full_alignment")]
    FullAlignment,
    #[value(name = "full_alignment_bloom")]
    FullAlignmentBloom,
}

impl fmt::Display for ModelKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Pileup => "pileup",
            Self::FullAlignment => "full_alignment",
            Self::FullAlignmentBloom => "full_alignment_bloom",
        })
    }
}

impl ModelKind {
    pub const fn default_input_channels(self) -> i64 {
        match self {
            Self::Pileup => 18,
            Self::FullAlignment | Self::FullAlignmentBloom => 8,
        }
    }
}

pub enum Network {
    Pileup(PileupModel),
    FullAlignment(FullAlignmentModel),
    Bloom(BloomModel),
}

pub struct Model {
    pub var_store: nn::VarStore,
    network: Network,
}

impl Model {
    pub fn new(
        kind: ModelKind,
        device: Device,
        input_channels: Option<i64>,
        add_indel_length: bool,
        checkpoint: Option<&Path>,
    ) -> Result<Self> {
        let mut var_store = nn::VarStore::new(device);
        let input_channels = input_channels.unwrap_or_else(|| kind.default_input_channels());
        let network = match kind {
            ModelKind::Pileup => Network::Pileup(PileupModel::new(
                &var_store.root(),
                add_indel_length,
                input_channels,
            )),
            ModelKind::FullAlignment => Network::FullAlignment(FullAlignmentModel::new(
                &var_store.root(),
                add_indel_length,
                input_channels,
            )),
            ModelKind::FullAlignmentBloom => {
                Network::Bloom(BloomModel::new(&var_store.root(), input_channels))
            }
        };
        if let Some(path) = checkpoint {
            load_checkpoint(&mut var_store, path)?;
        }
        Ok(Self { var_store, network })
    }

    pub fn forward_tasks(&self, input: &Tensor, train: bool) -> Vec<Tensor> {
        match &self.network {
            Network::Pileup(model) => model.forward_t(input, train),
            Network::FullAlignment(model) => model.forward_t(input, train),
            Network::Bloom(model) => model.forward_t(input, train),
        }
    }

    pub fn forward(&self, input: &Tensor, train: bool) -> Tensor {
        Tensor::cat(&self.forward_tasks(input, train), 1)
    }
}

pub fn parse_device(value: &str) -> Result<Device> {
    match value {
        "cpu" => Ok(Device::Cpu),
        "mps" => Ok(Device::Mps),
        "cuda" => Ok(Device::cuda_if_available()),
        value if value.starts_with("cuda:") => {
            let index = value[5..].parse()?;
            Ok(Device::Cuda(index))
        }
        _ => bail!("unsupported device '{value}', expected cpu, cuda, cuda:N, or mps"),
    }
}
