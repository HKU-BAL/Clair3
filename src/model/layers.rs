use tch::{
    Kind, Tensor,
    nn::{self, Init, Module, ModuleT},
};

pub fn uniform_init(fan_in: i64) -> Init {
    let bound = 1.0 / (fan_in as f64).sqrt();
    Init::Uniform {
        lo: -bound,
        up: bound,
    }
}

pub fn linear(vs: &nn::Path, in_dim: i64, out_dim: i64) -> nn::Linear {
    nn::linear(
        vs,
        in_dim,
        out_dim,
        nn::LinearConfig {
            ws_init: uniform_init(in_dim),
            bs_init: None,
            bias: true,
        },
    )
}

#[allow(clippy::too_many_arguments)]
pub fn conv2d(
    vs: &nn::Path,
    in_channels: i64,
    out_channels: i64,
    kernel: [i64; 2],
    stride: [i64; 2],
    padding: [i64; 2],
    groups: i64,
    bias: bool,
) -> nn::Conv2D {
    let fan_in = in_channels / groups * kernel[0] * kernel[1];
    nn::conv(
        vs,
        in_channels,
        out_channels,
        kernel,
        nn::ConvConfigND::<[i64; 2]> {
            stride,
            padding,
            dilation: [1, 1],
            groups,
            bias,
            ws_init: uniform_init(fan_in),
            bs_init: if bias {
                uniform_init(fan_in)
            } else {
                Init::Const(0.0)
            },
            padding_mode: nn::PaddingMode::Zeros,
        },
    )
}

pub struct BatchNorm2d {
    layer: nn::BatchNorm,
    num_batches_tracked: Tensor,
}

impl BatchNorm2d {
    pub fn new(vs: &nn::Path, channels: i64) -> Self {
        let layer = nn::batch_norm2d(
            vs,
            channels,
            nn::BatchNormConfig {
                eps: 1e-3,
                ws_init: Init::Const(1.0),
                bs_init: Init::Const(0.0),
                ..Default::default()
            },
        );
        let num_batches_tracked = vs.add(
            "num_batches_tracked",
            Tensor::zeros([], (Kind::Int64, vs.device())),
            false,
        );
        Self {
            layer,
            num_batches_tracked,
        }
    }

    pub fn forward_t(&self, input: &Tensor, train: bool) -> Tensor {
        if train {
            tch::no_grad(|| {
                let mut tracked = self.num_batches_tracked.shallow_clone();
                tracked.copy_(&(tracked.shallow_clone() + 1));
            });
        }
        self.layer.forward_t(input, train)
    }
}

pub enum ConvBlock {
    Standard(super::nnpack::NnpConv2D),
    Separable {
        depthwise: nn::Conv2D,
        pointwise: nn::Conv2D,
    },
}

impl ConvBlock {
    fn new(
        vs: &nn::Path,
        in_channels: i64,
        out_channels: i64,
        stride: [i64; 2],
        separable: bool,
    ) -> Self {
        if separable {
            let depthwise = conv2d(
                &(vs / "depthwise"),
                in_channels,
                in_channels,
                [3, 3],
                stride,
                [1, 1],
                in_channels,
                false,
            );
            let pointwise = conv2d(
                &(vs / "pointwise"),
                in_channels,
                out_channels,
                [1, 1],
                [1, 1],
                [0, 0],
                1,
                true,
            );
            Self::Separable {
                depthwise,
                pointwise,
            }
        } else {
            Self::Standard(super::nnpack::NnpConv2D::new(
                conv2d(
                    vs,
                    in_channels,
                    out_channels,
                    [3, 3],
                    stride,
                    [1, 1],
                    1,
                    true,
                ),
                stride,
                [1, 1],
            ))
        }
    }

    pub fn forward_t(&self, input: &Tensor, train: bool) -> Tensor {
        match self {
            Self::Standard(layer) if !train => layer
                .forward(input)
                .unwrap_or_else(|| layer.layer.forward(input)),
            Self::Standard(layer) => {
                layer.invalidate();
                layer.layer.forward(input)
            }
            Self::Separable {
                depthwise,
                pointwise,
            } => pointwise.forward(&depthwise.forward(input)),
        }
    }
}

pub struct BasicConv2d {
    conv: ConvBlock,
    bn: BatchNorm2d,
}

impl BasicConv2d {
    pub fn new(vs: &nn::Path, in_channels: i64, out_channels: i64, separable: bool) -> Self {
        Self {
            conv: ConvBlock::new(&(vs / "conv"), in_channels, out_channels, [2, 2], separable),
            bn: BatchNorm2d::new(&(vs / "bn"), out_channels),
        }
    }

    pub fn forward_t(&self, input: &Tensor, train: bool) -> Tensor {
        self.bn
            .forward_t(&self.conv.forward_t(input, train), train)
            .relu()
    }
}

pub struct ResidualBlock {
    conv1: ConvBlock,
    bn1: BatchNorm2d,
    conv2: ConvBlock,
    bn2: BatchNorm2d,
    downsample: Option<(nn::Conv2D, BatchNorm2d)>,
}

impl ResidualBlock {
    pub fn new(
        vs: &nn::Path,
        in_channels: i64,
        out_channels: i64,
        stride: [i64; 2],
        separable: bool,
    ) -> Self {
        let downsample = if stride != [1, 1] || in_channels != out_channels {
            Some((
                conv2d(
                    &(vs / "downsample" / "0"),
                    in_channels,
                    out_channels,
                    [1, 1],
                    stride,
                    [0, 0],
                    1,
                    true,
                ),
                BatchNorm2d::new(&(vs / "downsample" / "1"), out_channels),
            ))
        } else {
            None
        };
        Self {
            conv1: ConvBlock::new(
                &(vs / "conv1"),
                in_channels,
                out_channels,
                stride,
                separable,
            ),
            bn1: BatchNorm2d::new(&(vs / "bn1"), out_channels),
            conv2: ConvBlock::new(
                &(vs / "conv2"),
                out_channels,
                out_channels,
                [1, 1],
                separable,
            ),
            bn2: BatchNorm2d::new(&(vs / "bn2"), out_channels),
            downsample,
        }
    }

    pub fn forward_t(&self, input: &Tensor, train: bool) -> Tensor {
        let residual = match &self.downsample {
            Some((conv, bn)) => bn.forward_t(&conv.forward(input), train),
            None => input.shallow_clone(),
        };
        let mut output = self
            .bn1
            .forward_t(&self.conv1.forward_t(input, train), train)
            .relu();
        output = self
            .bn2
            .forward_t(&self.conv2.forward_t(&output, train), train);
        (residual + output).relu()
    }
}

pub struct IndelHeads {
    l5_3: nn::Linear,
    l5_4: nn::Linear,
    y_length_1: nn::Linear,
    y_length_2: nn::Linear,
}

pub struct ClassificationHeads {
    l4: nn::Linear,
    l5_1: nn::Linear,
    l5_2: nn::Linear,
    y_gt21: nn::Linear,
    y_genotype: nn::Linear,
    indel: Option<IndelHeads>,
}

impl ClassificationHeads {
    pub fn new(vs: &nn::Path, l4_input: i64, l4_units: i64, add_indel_length: bool) -> Self {
        let indel = add_indel_length.then(|| IndelHeads {
            l5_3: linear(&(vs / "L5_3"), l4_units, 128),
            l5_4: linear(&(vs / "L5_4"), l4_units, 128),
            y_length_1: linear(&(vs / "Y_indel_length_logits_1"), 128, 33),
            y_length_2: linear(&(vs / "Y_indel_length_logits_2"), 128, 33),
        });
        Self {
            l4: linear(&(vs / "L4"), l4_input, l4_units),
            l5_1: linear(&(vs / "L5_1"), l4_units, 128),
            l5_2: linear(&(vs / "L5_2"), l4_units, 128),
            y_gt21: linear(&(vs / "Y_gt21_logits"), 128, 21),
            y_genotype: linear(&(vs / "Y_genotype_logits"), 128, 3),
            indel,
        }
    }

    pub fn forward_t(&self, input: &Tensor, train: bool) -> Vec<Tensor> {
        let input = self.l4.forward(input).selu().dropout(0.5, train);
        let first = self.l5_1.forward(&input).selu().dropout(0.2, train);
        let second = self.l5_2.forward(&input).selu().dropout(0.2, train);
        let mut outputs = vec![
            self.y_gt21.forward(&first).selu().softmax(-1, Kind::Float),
            self.y_genotype
                .forward(&second)
                .selu()
                .softmax(-1, Kind::Float),
        ];
        if let Some(indel) = &self.indel {
            let third = indel.l5_3.forward(&input).selu().dropout(0.2, train);
            let fourth = indel.l5_4.forward(&input).selu().dropout(0.2, train);
            outputs.push(
                indel
                    .y_length_1
                    .forward(&third)
                    .selu()
                    .softmax(-1, Kind::Float),
            );
            outputs.push(
                indel
                    .y_length_2
                    .forward(&fourth)
                    .selu()
                    .softmax(-1, Kind::Float),
            );
        }
        outputs
    }
}

pub fn xavier_uniform(vs: &nn::Path, name: &str, dims: &[i64]) -> Tensor {
    let fan_in = dims[1];
    let fan_out = dims[0];
    let bound = (6.0 / (fan_in + fan_out) as f64).sqrt();
    vs.var(
        name,
        dims,
        Init::Uniform {
            lo: -bound,
            up: bound,
        },
    )
}

pub fn zeros(vs: &nn::Path, name: &str, dims: &[i64]) -> Tensor {
    vs.zeros(name, dims)
}
