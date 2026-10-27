use tch::{
    Tensor,
    nn::{self, RNN},
};

use super::layers::{ClassificationHeads, uniform_init};

pub struct PileupModel {
    lstm1: nn::LSTM,
    lstm2: nn::LSTM,
    heads: ClassificationHeads,
}

impl PileupModel {
    pub fn new(vs: &nn::Path, add_indel_length: bool, input_channels: i64) -> Self {
        let rnn_config = |hidden| nn::RNNConfig {
            bidirectional: true,
            w_ih_init: uniform_init(hidden),
            w_hh_init: uniform_init(hidden),
            b_ih_init: Some(uniform_init(hidden)),
            b_hh_init: Some(uniform_init(hidden)),
            ..Default::default()
        };
        Self {
            lstm1: nn::lstm(&(vs / "LSTM1"), input_channels, 128, rnn_config(128)),
            lstm2: nn::lstm(&(vs / "LSTM2"), 256, 160, rnn_config(160)),
            heads: ClassificationHeads::new(vs, 160 * 2 * 33, 128, add_indel_length),
        }
    }

    pub fn forward_t(&self, input: &Tensor, train: bool) -> Vec<Tensor> {
        let input = input.to_kind(tch::Kind::Float);
        let (input, _) = self.lstm1.seq(&input);
        let (input, _) = self.lstm2.seq(&input);
        let input = input.dropout(0.2, train).flatten(1, -1);
        self.heads.forward_t(&input, train)
    }
}
