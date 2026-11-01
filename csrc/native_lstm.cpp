#include <cstddef>
#include <cstdint>

#ifdef CLAIR3_NATIVE_LSTM
#include <Accelerate/Accelerate.h>
#include <ATen/cpu/vec/vec.h>
#include <cmath>
#include <cstdlib>
#include <vector>

extern "C" void clair3_lstm_prepare(void) {
    setenv("VECLIB_MAXIMUM_THREADS", "1", 0);
}

extern "C" int32_t clair3_lstm_forward(
    const float *input,
    const float *input_weights,
    const float *bias_ih,
    const float *bias_hh,
    const float *forward_weights,
    const float *reverse_weights,
    float *output,
    float *input_gates,
    float *hidden_gates,
    float *states,
    int32_t batch,
    int32_t time,
    int32_t features,
    int32_t hidden
) {
    if (!input || !input_weights || !bias_ih || !bias_hh || !forward_weights || !reverse_weights ||
        !output || !input_gates || !hidden_gates || !states || batch <= 0 || time <= 0 ||
        features <= 0 || hidden <= 0) {
        return 1;
    }
    using Vec = at::vec::Vectorized<float>;
    const auto sigmoid = [](const Vec& value) {
        return Vec(1.0f) / (Vec(1.0f) + (-value).exp());
    };
    const auto gate = [](const float *projected, const float *recurrent,
                        const float *input_bias, const float *hidden_bias, int32_t offset) {
        const auto input_term = Vec::loadu(projected + offset) + Vec::loadu(input_bias + offset);
        const auto hidden_term = Vec::loadu(recurrent + offset) + Vec::loadu(hidden_bias + offset);
        return input_term + hidden_term;
    };
    const int32_t directions = 2;
    const int32_t gate_count = 4 * hidden;
    const int32_t combined_gates = directions * gate_count;
    cblas_sgemm(
        CblasRowMajor, CblasNoTrans, CblasTrans, batch * time, combined_gates, features,
        1.0f, input, features, input_weights, features, 0.0f, input_gates, combined_gates
    );
    for (int32_t direction = 0; direction < directions; ++direction) {
        for (int32_t index = 0; index < directions * batch * hidden; ++index) {
            states[index] = 0.0f;
        }
        const float *weights = direction == 0 ? forward_weights : reverse_weights;
        float *h = states;
        float *c = states + static_cast<size_t>(batch) * hidden;
        for (int32_t step = 0; step < time; ++step) {
            const int32_t current = direction == 0 ? step : time - step - 1;
            cblas_sgemm(
                CblasRowMajor, CblasNoTrans, CblasTrans, batch, gate_count, hidden,
                1.0f, h, hidden, weights, hidden, 0.0f, hidden_gates, gate_count
            );
            for (int32_t row = 0; row < batch; ++row) {
                const float *projected = input_gates +
                    (static_cast<size_t>(row) * time + current) * combined_gates +
                    direction * gate_count;
                const float *recurrent = hidden_gates + static_cast<size_t>(row) * gate_count;
                const float *input_bias = bias_ih + direction * gate_count;
                const float *hidden_bias = bias_hh + direction * gate_count;
                float *row_h = h + static_cast<size_t>(row) * hidden;
                float *row_c = c + static_cast<size_t>(row) * hidden;
                float *row_output = output +
                    (static_cast<size_t>(current) * batch + row) * directions * hidden +
                    direction * hidden;
                int32_t column = 0;
                for (; column + Vec::size() <= hidden; column += Vec::size()) {
                    const auto input_gate = sigmoid(
                        gate(projected, recurrent, input_bias, hidden_bias, column)
                    );
                    const auto forget_gate = sigmoid(
                        gate(projected, recurrent, input_bias, hidden_bias, hidden + column)
                    );
                    const auto candidate = gate(
                        projected, recurrent, input_bias, hidden_bias, 2 * hidden + column
                    ).tanh();
                    const auto output_gate = sigmoid(
                        gate(projected, recurrent, input_bias, hidden_bias, 3 * hidden + column)
                    );
                    const auto previous = forget_gate * Vec::loadu(row_c + column);
                    const auto update = input_gate * candidate;
                    const auto cell = previous + update;
                    const auto next_h = output_gate * cell.tanh();
                    cell.store(row_c + column);
                    next_h.store(row_h + column);
                    next_h.store(row_output + column);
                }
                for (; column < hidden; ++column) {
                    const float input_gate = 1.0f / (1.0f + std::exp(-(
                        (projected[column] + input_bias[column]) +
                        (recurrent[column] + hidden_bias[column])
                    )));
                    const float forget_gate = 1.0f / (1.0f + std::exp(-(
                        (projected[hidden + column] + input_bias[hidden + column]) +
                        (recurrent[hidden + column] + hidden_bias[hidden + column])
                    )));
                    const float candidate = std::tanh(
                        (projected[2 * hidden + column] + input_bias[2 * hidden + column]) +
                        (recurrent[2 * hidden + column] + hidden_bias[2 * hidden + column])
                    );
                    const float output_gate = 1.0f / (1.0f + std::exp(-(
                        (projected[3 * hidden + column] + input_bias[3 * hidden + column]) +
                        (recurrent[3 * hidden + column] + hidden_bias[3 * hidden + column])
                    )));
                    const float previous = forget_gate * row_c[column];
                    const float update = input_gate * candidate;
                    const float cell = previous + update;
                    const float next_h = output_gate * std::tanh(cell);
                    row_c[column] = cell;
                    row_h[column] = next_h;
                    row_output[column] = next_h;
                }
            }
        }
    }
    return 0;
}

extern "C" int32_t clair3_gemm_conv(
    const float *input,
    const float *weights,
    const float *bias,
    float *output,
    float *columns,
    int32_t batch,
    int32_t channels,
    int32_t height,
    int32_t width,
    int32_t output_channels,
    int32_t stride,
    int32_t padding
) {
    if (!input || !weights || !bias || !output || !columns || batch <= 0 || channels <= 0 ||
        height <= 0 || width <= 0 || output_channels <= 0 || stride <= 0 || padding < 0) {
        return 1;
    }
    clair3_lstm_prepare();
    const int32_t output_height = (height + 2 * padding - 3) / stride + 1;
    const int32_t output_width = (width + 2 * padding - 3) / stride + 1;
    if (output_height <= 0 || output_width <= 0) {
        return 1;
    }
    const size_t rows = static_cast<size_t>(batch) * output_height * output_width;
    const size_t kernel = static_cast<size_t>(channels) * 9;
    for (int32_t sample = 0; sample < batch; ++sample) {
        for (int32_t oy = 0; oy < output_height; ++oy) {
            for (int32_t ox = 0; ox < output_width; ++ox) {
                const size_t row = (static_cast<size_t>(sample) * output_height + oy) * output_width + ox;
                for (int32_t channel = 0; channel < channels; ++channel) {
                    for (int32_t ky = 0; ky < 3; ++ky) {
                        const int32_t y = oy * stride + ky - padding;
                        for (int32_t kx = 0; kx < 3; ++kx) {
                            const int32_t x = ox * stride + kx - padding;
                            if (y >= 0 && y < height && x >= 0 && x < width) {
                                const size_t column = (static_cast<size_t>(channel) * 3 + ky) * 3 + kx;
                                columns[row * kernel + column] = input[
                                    (static_cast<size_t>(sample) * channels + channel) * height * width +
                                    static_cast<size_t>(y) * width + x
                                ];
                            }
                        }
                    }
                }
            }
        }
    }
    std::vector<float> product(rows * output_channels);
    cblas_sgemm(
        CblasRowMajor, CblasNoTrans, CblasTrans, rows, output_channels, kernel,
        1.0f, columns, kernel, weights, kernel, 0.0f, product.data(), output_channels
    );
    const size_t sample_rows = static_cast<size_t>(output_height) * output_width;
    for (int32_t sample = 0; sample < batch; ++sample) {
        for (int32_t channel = 0; channel < output_channels; ++channel) {
            for (size_t row = 0; row < sample_rows; ++row) {
                const size_t input_row = static_cast<size_t>(sample) * sample_rows + row;
                const size_t output_row =
                    (static_cast<size_t>(sample) * output_channels + channel) * sample_rows + row;
                output[output_row] = product[input_row * output_channels + channel] + bias[channel];
            }
        }
    }
    return 0;
}
#else
extern "C" void clair3_lstm_prepare(void) {}

extern "C" int32_t clair3_lstm_forward(
    const float *input, const float *input_weights, const float *bias_ih,
    const float *bias_hh, const float *forward_weights, const float *reverse_weights,
    float *output, float *input_gates, float *hidden_gates, float *states,
    int32_t batch, int32_t time, int32_t features, int32_t hidden
) {
    (void)input; (void)input_weights; (void)bias_ih; (void)bias_hh; (void)forward_weights;
    (void)reverse_weights; (void)output; (void)input_gates; (void)hidden_gates;
    (void)states; (void)batch; (void)time; (void)features; (void)hidden;
    return 1;
}

extern "C" int32_t clair3_gemm_conv(
    const float *input, const float *weights, const float *bias, float *output,
    float *columns, int32_t batch, int32_t channels, int32_t height, int32_t width,
    int32_t output_channels, int32_t stride, int32_t padding
) {
    (void)input; (void)weights; (void)bias; (void)output; (void)columns;
    (void)batch; (void)channels; (void)height; (void)width; (void)output_channels;
    (void)stride; (void)padding;
    return 1;
}
#endif
