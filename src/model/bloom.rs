use tch::{
    Kind, Tensor,
    nn::{self, Module},
};

use super::layers::{BasicConv2d, ResidualBlock, linear, xavier_uniform, zeros};

struct MultiheadAttention {
    in_proj_weight: Tensor,
    in_proj_bias: Tensor,
    out_proj_weight: Tensor,
    out_proj_bias: Tensor,
}

impl MultiheadAttention {
    fn new(vs: &nn::Path) -> Self {
        Self {
            in_proj_weight: xavier_uniform(vs, "in_proj_weight", &[288, 96]),
            in_proj_bias: zeros(vs, "in_proj_bias", &[288]),
            out_proj_weight: xavier_uniform(&(vs / "out_proj"), "weight", &[96, 96]),
            out_proj_bias: zeros(&(vs / "out_proj"), "bias", &[96]),
        }
    }

    fn forward_t(&self, input: &Tensor, train: bool) -> Tensor {
        let sizes = input.size();
        let (batch, sequence) = (sizes[0], sizes[1]);
        let projected = input.linear(&self.in_proj_weight, Some(&self.in_proj_bias));
        let chunks = projected.chunk(3, -1);
        let reshape = |tensor: &Tensor| {
            tensor
                .reshape([batch, sequence, 4, 24])
                .permute([0, 2, 1, 3])
        };
        let query = reshape(&chunks[0]);
        let key = reshape(&chunks[1]);
        let value = reshape(&chunks[2]);
        let scores = query.matmul(&key.transpose(-2, -1)) / 24.0_f64.sqrt();
        let probabilities = scores.softmax(-1, Kind::Float).dropout(0.1, train);
        let context = probabilities
            .matmul(&value)
            .permute([0, 2, 1, 3])
            .contiguous()
            .reshape([batch, sequence, 96]);
        context
            .linear(&self.out_proj_weight, Some(&self.out_proj_bias))
            .dropout(0.1, train)
    }
}

pub struct BloomModel {
    input_channels: i64,
    conv1: BasicConv2d,
    res1: ResidualBlock,
    conv2: BasicConv2d,
    res2: ResidualBlock,
    conv3: BasicConv2d,
    res3: ResidualBlock,
    attn_proj: nn::Conv2D,
    attn_mha: MultiheadAttention,
    attn_norm: nn::LayerNorm,
    dense1: nn::Linear,
    dense2: nn::Linear,
    logit: nn::Linear,
}

impl BloomModel {
    pub fn new(vs: &nn::Path, input_channels: i64) -> Self {
        use super::layers::conv2d;
        Self {
            input_channels,
            conv1: BasicConv2d::new(&(vs / "conv1"), input_channels, 48, true),
            res1: ResidualBlock::new(&(vs / "res1"), 48, 48, [1, 1], true),
            conv2: BasicConv2d::new(&(vs / "conv2"), 48, 96, true),
            res2: ResidualBlock::new(&(vs / "res2"), 96, 96, [1, 1], true),
            conv3: BasicConv2d::new(&(vs / "conv3"), 96, 160, true),
            res3: ResidualBlock::new(&(vs / "res3"), 160, 160, [1, 1], true),
            attn_proj: conv2d(
                &(vs / "attn_proj"),
                160,
                96,
                [1, 1],
                [1, 1],
                [0, 0],
                1,
                true,
            ),
            attn_mha: MultiheadAttention::new(&(vs / "attn_mha")),
            attn_norm: nn::layer_norm(
                &(vs / "attn_norm"),
                vec![96],
                nn::LayerNormConfig {
                    eps: 1e-5,
                    ..Default::default()
                },
            ),
            dense1: linear(&(vs / "dense1"), 96 * 2, 128),
            dense2: linear(&(vs / "dense2"), 128, 96),
            logit: linear(&(vs / "logit"), 96, 1),
        }
    }

    pub fn forward_t(&self, input: &Tensor, train: bool) -> Vec<Tensor> {
        let mut input = input.to_kind(Kind::Float) / 100.0;
        if input.dim() == 4 && input.size()[3] == self.input_channels {
            input = input.permute([0, 3, 1, 2]);
        }
        let mut output = self.conv1.forward_t(&input, train);
        output = self.res1.forward_t(&output, train);
        output = self.conv2.forward_t(&output, train);
        output = self.res2.forward_t(&output, train);
        output = self.conv3.forward_t(&output, train);
        output = self.res3.forward_t(&output, train);

        let attention_input = self.attn_proj.forward(&output).selu();
        let sizes = attention_input.size();
        let (batch, channels, height, width) = (sizes[0], sizes[1], sizes[2], sizes[3]);
        let sequence =
            attention_input
                .permute([0, 2, 3, 1])
                .reshape([batch, height * width, channels]);
        let attention = self.attn_mha.forward_t(&sequence, train);
        let sequence = self.attn_norm.forward(&(sequence + attention));
        let attention_map = sequence
            .reshape([batch, height, width, channels])
            .permute([0, 3, 1, 2]);
        let combined = Tensor::cat(&[&output, &attention_map], 1);
        let pooled = Tensor::cat(
            &[
                combined.mean_dim(&[2, 3][..], false, Kind::Float),
                combined.amax(&[2, 3][..], false),
            ],
            1,
        );
        let output = self.dense1.forward(&pooled.dropout(0.35, train)).selu();
        let output = self.dense2.forward(&output.dropout(0.25, train)).selu();
        vec![self.logit.forward(&output).sigmoid()]
    }
}
