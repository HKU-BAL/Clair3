use std::{fmt, str::FromStr};

use clap::ValueEnum;

pub const FLANKING_BASES: usize = 16;
pub const POSITION_COUNT: usize = FLANKING_BASES * 2 + 1;
pub const PILEUP_CHANNELS: usize = 18;
pub const FULL_ALIGNMENT_CHANNELS: usize = 8;
pub const DWELL_CHANNELS: usize = FULL_ALIGNMENT_CHANNELS + 1;
pub const MAX_DEPTH: usize = 144;
pub const MAX_VARIANT_LENGTH: usize = 50;
pub const MAX_LONG_VARIANT_LENGTH: usize = 100_000;
pub const LABEL_COUNT: usize = 90;
pub const LABEL_COUNT_WITHOUT_INDEL_LENGTH: usize = 24;
pub const PREDICT_BATCH_SIZE: usize = 200;
pub const FULL_ALIGNMENT_PREDICT_BATCH_SIZE: usize = 200;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, ValueEnum)]
pub enum Platform {
    #[value(name = "ont")]
    Ont,
    #[value(name = "hifi")]
    Hifi,
    #[value(name = "ilmn")]
    Ilmn,
}

impl Platform {
    pub const fn matrix_depth(self) -> usize {
        match self {
            Self::Ont => 89,
            Self::Hifi | Self::Ilmn => 55,
        }
    }

    pub const fn min_indel_af(self) -> f64 {
        match self {
            Self::Ont => 0.15,
            Self::Hifi | Self::Ilmn => 0.08,
        }
    }

    pub const fn default_var_pct_full(self) -> f64 {
        match self {
            Self::Ont => 0.7,
            Self::Hifi | Self::Ilmn => 0.3,
        }
    }

    pub const fn default_ref_pct_full(self) -> f64 {
        match self {
            Self::Ont => 0.1,
            Self::Hifi | Self::Ilmn => 0.3,
        }
    }
}

impl fmt::Display for Platform {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Ont => "ont",
            Self::Hifi => "hifi",
            Self::Ilmn => "ilmn",
        })
    }
}

impl FromStr for Platform {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "ont" => Ok(Self::Ont),
            "hifi" => Ok(Self::Hifi),
            "ilmn" => Ok(Self::Ilmn),
            _ => Err(format!(
                "unsupported platform '{value}', expected ont, hifi, or ilmn"
            )),
        }
    }
}
