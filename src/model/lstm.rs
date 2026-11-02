use std::collections::HashMap;

use tch::{Device, Kind, Tensor};

pub struct NativeLstm {
    input_weights: Tensor,
    bias_ih: Tensor,
    bias_hh: Tensor,
    forward_weights: Tensor,
    reverse_weights: Tensor,
    features: i64,
    hidden: i64,
}

impl NativeLstm {
    pub fn new(
        weights: &HashMap<String, Tensor>,
        prefix: &str,
        features: i64,
        hidden: i64,
    ) -> Option<Self> {
        unsafe extern "C" {
            fn clair3_lstm_prepare();
        }
        unsafe { clair3_lstm_prepare() };
        let tensor = |name: &str| {
            weights
                .get(&format!("{prefix}{name}"))
                .map(Tensor::shallow_clone)
        };
        let forward_input = tensor("weight_ih_l0")?;
        let reverse_input = tensor("weight_ih_l0_reverse")?;
        let forward_bias_ih = tensor("bias_ih_l0")?;
        let forward_bias_hh = tensor("bias_hh_l0")?;
        let reverse_bias_ih = tensor("bias_ih_l0_reverse")?;
        let reverse_bias_hh = tensor("bias_hh_l0_reverse")?;
        let input_weights = Tensor::cat(&[forward_input, reverse_input], 0).contiguous();
        let bias_ih = Tensor::cat(&[forward_bias_ih, reverse_bias_ih], 0).contiguous();
        let bias_hh = Tensor::cat(&[forward_bias_hh, reverse_bias_hh], 0).contiguous();
        let forward_weights = tensor("weight_hh_l0")?.contiguous();
        let reverse_weights = tensor("weight_hh_l0_reverse")?.contiguous();
        let tensors = [
            &input_weights,
            &bias_ih,
            &bias_hh,
            &forward_weights,
            &reverse_weights,
        ];
        if tensors
            .iter()
            .any(|tensor| tensor.device() != Device::Cpu || tensor.kind() != Kind::Float)
        {
            return None;
        }
        if input_weights.size() != [8 * hidden, features]
            || bias_ih.size() != [8 * hidden]
            || bias_hh.size() != [8 * hidden]
            || forward_weights.size() != [4 * hidden, hidden]
            || reverse_weights.size() != [4 * hidden, hidden]
        {
            return None;
        }
        Some(Self {
            input_weights,
            bias_ih,
            bias_hh,
            forward_weights,
            reverse_weights,
            features,
            hidden,
        })
    }

    pub fn forward(&self, input: &Tensor) -> Option<Tensor> {
        let sizes = input.size();
        if sizes.len() != 3
            || sizes[2] != self.features
            || sizes[0] <= 0
            || sizes[1] <= 0
            || input.device() != Device::Cpu
            || input.kind() != Kind::Float
        {
            return None;
        }
        let input = input.contiguous();
        let (batch, time) = (sizes[0], sizes[1]);
        let output = Tensor::empty([time, batch, 2 * self.hidden], (Kind::Float, Device::Cpu));
        let mut input_gates = vec![0_f32; (batch * time * 8 * self.hidden) as usize];
        let mut hidden_gates = vec![0_f32; (batch * 4 * self.hidden) as usize];
        let mut states = vec![0_f32; (batch * 2 * self.hidden) as usize];
        unsafe extern "C" {
            fn clair3_lstm_forward(
                input: *const f32,
                input_weights: *const f32,
                bias_ih: *const f32,
                bias_hh: *const f32,
                forward_weights: *const f32,
                reverse_weights: *const f32,
                output: *mut f32,
                input_gates: *mut f32,
                hidden_gates: *mut f32,
                states: *mut f32,
                batch: i32,
                time: i32,
                features: i32,
                hidden: i32,
            ) -> i32;
        }
        let status = unsafe {
            clair3_lstm_forward(
                input.data_ptr().cast(),
                self.input_weights.data_ptr().cast(),
                self.bias_ih.data_ptr().cast(),
                self.bias_hh.data_ptr().cast(),
                self.forward_weights.data_ptr().cast(),
                self.reverse_weights.data_ptr().cast(),
                output.data_ptr().cast(),
                input_gates.as_mut_ptr(),
                hidden_gates.as_mut_ptr(),
                states.as_mut_ptr(),
                batch as i32,
                time as i32,
                self.features as i32,
                self.hidden as i32,
            )
        };
        if status == 0 {
            Some(output.transpose(0, 1))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tch::nn::{self, RNN};

    #[test]
    fn native_lstm_matches_libtorch() {
        let store = nn::VarStore::new(Device::Cpu);
        let layer = nn::lstm(
            store.root(),
            3,
            4,
            nn::RNNConfig {
                bidirectional: true,
                b_ih_init: Some(nn::Init::Const(0.125)),
                b_hh_init: Some(nn::Init::Const(-0.0625)),
                ..Default::default()
            },
        );
        let native = NativeLstm::new(&store.variables(), "", 3, 4).unwrap();
        for batch in [1, 2, 4] {
            for time in [1, 2, 5] {
                let input = Tensor::randn([batch, time, 3], (Kind::Float, Device::Cpu));
                let (expected, _) = layer.seq(&input);
                let Some(actual) = native.forward(&input) else {
                    return;
                };
                let difference = (actual - expected).abs().max().double_value(&[]);
                assert!(
                    difference < 0.0001,
                    "batch {batch}, time {time}: {difference}"
                );
            }
        }
    }

    #[test]
    fn deterministic_lstm_layout() {
        let store = nn::VarStore::new(Device::Cpu);
        let layer = nn::lstm(
            store.root(),
            3,
            4,
            nn::RNNConfig {
                bidirectional: true,
                ..Default::default()
            },
        );
        tch::no_grad(|| {
            for mut tensor in store.variables().into_values() {
                let _ = tensor.fill_(0.0);
            }
            let weights = store.variables();
            let _ = weights["weight_ih_l0"]
                .shallow_clone()
                .narrow(0, 8, 1)
                .narrow(1, 0, 1)
                .fill_(1.0);
            let _ = weights["weight_ih_l0_reverse"]
                .shallow_clone()
                .narrow(0, 8, 1)
                .narrow(1, 1, 1)
                .fill_(1.0);
        });
        let input = Tensor::from_slice(&[
            1.0_f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        ])
        .reshape([2, 2, 3]);
        let native = NativeLstm::new(&store.variables(), "", 3, 4).unwrap();
        let (expected, _) = layer.seq(&input);
        let Some(actual) = native.forward(&input) else {
            return;
        };
        let difference = (actual - expected).abs().max().double_value(&[]);
        assert!(difference < 0.0001, "maximum difference {difference}");
    }

    #[test]
    fn official_lstm_layers_match() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let checkpoint = root.join("target/test-data/models/ont/pileup.pt");
        if !checkpoint.is_file() {
            return;
        }
        let model = crate::model::Model::new(
            crate::model::ModelKind::Pileup,
            Device::Cpu,
            Some(18),
            false,
            Some(&checkpoint),
        )
        .unwrap();
        let weights = model.var_store.variables();
        let mut expected =
            Tensor::read_npy(root.join("tests/fixtures/models/npy/pileup_ont_input.npy"))
                .unwrap()
                .to_kind(Kind::Float);
        for (prefix, features, hidden) in [("LSTM1.", 18, 128), ("LSTM2.", 256, 160)] {
            let native = NativeLstm::new(&weights, prefix, features, hidden).unwrap();
            let Some(actual) = native.forward(&expected) else {
                return;
            };
            let names = [
                "weight_ih_l0",
                "weight_hh_l0",
                "bias_ih_l0",
                "bias_hh_l0",
                "weight_ih_l0_reverse",
                "weight_hh_l0_reverse",
                "bias_ih_l0_reverse",
                "bias_hh_l0_reverse",
            ];
            let parameters: Vec<_> = names
                .iter()
                .map(|name| weights[&format!("{prefix}{name}")].shallow_clone())
                .collect();
            let states = Tensor::zeros([2, expected.size()[0], hidden], (Kind::Float, Device::Cpu));
            let (reference, _, _) = expected.lstm(
                &[states.shallow_clone(), states],
                &parameters,
                true,
                1,
                0.0,
                true,
                true,
                true,
            );
            let difference = (actual.contiguous() - reference.contiguous())
                .abs()
                .max()
                .double_value(&[]);
            let first_difference = (actual.contiguous().narrow(1, 0, 1)
                - reference.contiguous().narrow(1, 0, 1))
            .abs()
            .max()
            .double_value(&[]);
            eprintln!(
                "{prefix} maximum difference {difference}, first timestep {first_difference}"
            );
            assert!(difference < 0.0001);
            expected = reference;
        }
    }
}
