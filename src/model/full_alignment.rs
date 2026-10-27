use tch::{Kind, Tensor, nn};

use super::layers::{BasicConv2d, ClassificationHeads, ResidualBlock};

pub struct FullAlignmentModel {
    input_channels: i64,
    conv1: BasicConv2d,
    res_block1: ResidualBlock,
    conv3: BasicConv2d,
    res_block2: ResidualBlock,
    conv5: BasicConv2d,
    res_block3: ResidualBlock,
    heads: ClassificationHeads,
}

impl FullAlignmentModel {
    pub fn new(vs: &nn::Path, add_indel_length: bool, input_channels: i64) -> Self {
        Self {
            input_channels,
            conv1: BasicConv2d::new(&(vs / "conv1"), input_channels, 64, false),
            res_block1: ResidualBlock::new(&(vs / "res_block1" / "0"), 64, 64, [1, 1], false),
            conv3: BasicConv2d::new(&(vs / "conv3"), 64, 128, false),
            res_block2: ResidualBlock::new(&(vs / "res_block2" / "0"), 128, 128, [1, 1], false),
            conv5: BasicConv2d::new(&(vs / "conv5"), 128, 256, false),
            res_block3: ResidualBlock::new(&(vs / "res_block3" / "0"), 256, 256, [1, 1], false),
            heads: ClassificationHeads::new(vs, 14 * 256, 256, add_indel_length),
        }
    }

    fn pyramid_pool(input: &Tensor) -> Tensor {
        let sizes = input.size();
        let height = sizes[2];
        let width = sizes[3];
        let mut pooled = Vec::with_capacity(3);
        for pool_size in [3, 2, 1] {
            let window_h = (height + pool_size - 1) / pool_size;
            let window_w = (width + pool_size - 1) / pool_size;
            let stride_h = window_h;
            let stride_w = window_w;
            let out_h = (height + stride_h - 1) / stride_h;
            let out_w = (width + stride_w - 1) / stride_w;
            let pad_h = ((out_h - 1) * stride_h + window_h - height).max(0);
            let pad_w = ((out_w - 1) * stride_w + window_w - width).max(0);
            let pad_top = pad_h / 2;
            let pad_bottom = pad_h - pad_top;
            let pad_left = pad_w / 2;
            let pad_right = pad_w - pad_left;
            let padded = if pad_h > 0 || pad_w > 0 {
                input.constant_pad_nd([pad_left, pad_right, pad_top, pad_bottom])
            } else {
                input.shallow_clone()
            };
            let pool = padded
                .max_pool2d(
                    [window_h, window_w],
                    [stride_h, stride_w],
                    [0, 0],
                    [1, 1],
                    false,
                )
                .permute([0, 2, 3, 1])
                .flatten(1, -1);
            pooled.push(pool);
        }
        Tensor::cat(&pooled, 1)
    }

    pub fn forward_t(&self, input: &Tensor, train: bool) -> Vec<Tensor> {
        let mut input = input.to_kind(Kind::Float) / 100.0;
        if input.dim() == 4 && input.size()[3] == self.input_channels {
            input = input.permute([0, 3, 1, 2]);
        }
        let mut output = self.conv1.forward_t(&input, train);
        output = self.res_block1.forward_t(&output, train);
        output = self.conv3.forward_t(&output, train);
        output = self.res_block2.forward_t(&output, train);
        output = self.conv5.forward_t(&output, train);
        output = self.res_block3.forward_t(&output, train);
        output = Self::pyramid_pool(&output).dropout(0.2, train);
        self.heads.forward_t(&output, train)
    }
}
