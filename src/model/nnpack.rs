use std::{
    ffi::c_void,
    sync::{Mutex, OnceLock},
};

use tch::{Device, Kind, Tensor, nn};

#[derive(Clone, Copy)]
#[repr(C)]
struct Size {
    width: usize,
    height: usize,
}

#[derive(Clone, Copy)]
#[repr(C)]
struct Padding {
    top: usize,
    right: usize,
    bottom: usize,
    left: usize,
}

type Inference = unsafe extern "C" fn(
    i32,
    i32,
    usize,
    usize,
    Size,
    Padding,
    Size,
    Size,
    *const f32,
    *const f32,
    *const f32,
    *mut f32,
    *mut c_void,
    *mut usize,
    i32,
    *const c_void,
    *mut c_void,
    *mut c_void,
) -> i32;
type CreatePool = unsafe extern "C" fn(usize) -> *mut c_void;
type DestroyPool = unsafe extern "C" fn(*mut c_void);
type Initialize = unsafe extern "C" fn() -> i32;

#[derive(Clone, Copy)]
struct Api {
    _initialize: Initialize,
    inference: Inference,
    create_pool: Option<CreatePool>,
    destroy_pool: Option<DestroyPool>,
}

#[cfg(unix)]
fn api() -> Option<&'static Api> {
    static API: OnceLock<Option<Api>> = OnceLock::new();
    API.get_or_init(|| unsafe {
        unsafe fn symbol<T: Copy>(name: &[u8]) -> Option<T> {
            let pointer = unsafe { libc::dlsym(libc::RTLD_DEFAULT, name.as_ptr().cast()) };
            if pointer.is_null() || std::mem::size_of::<T>() != std::mem::size_of::<*mut c_void>() {
                return None;
            }
            Some(unsafe { std::mem::transmute_copy(&pointer) })
        }
        let initialize: Initialize = symbol(b"nnp_initialize\0")?;
        if initialize() != 0 {
            return None;
        }
        Some(Api {
            _initialize: initialize,
            inference: symbol(b"nnp_convolution_inference\0")?,
            create_pool: symbol(b"pthreadpool_create\0"),
            destroy_pool: symbol(b"pthreadpool_destroy\0"),
        })
    })
    .as_ref()
}

#[cfg(not(unix))]
fn api() -> Option<&'static Api> {
    None
}

struct Workspace {
    tensor: Tensor,
    bytes: usize,
}

impl Workspace {
    fn new(bytes: usize) -> Option<Self> {
        if bytes == 0 || !bytes.is_multiple_of(4) {
            return None;
        }
        Some(Self {
            tensor: Tensor::empty([bytes as i64 / 4], (Kind::Float, Device::Cpu)),
            bytes,
        })
    }

    fn pointer(&mut self) -> *mut c_void {
        self.tensor.data_ptr()
    }
}

struct Cache {
    packed: Workspace,
    workspace: Workspace,
    pool: usize,
    input_size: [i64; 2],
}

impl Drop for Cache {
    fn drop(&mut self) {
        if self.pool != 0
            && let Some(destroy) = api().and_then(|api| api.destroy_pool)
        {
            unsafe { destroy(self.pool as *mut c_void) };
        }
    }
}

impl Cache {
    fn new(
        api: &Api,
        stride: [i64; 2],
        padding: [i64; 2],
        input: &Tensor,
        weight: &Tensor,
        bias: &Tensor,
        output: &Tensor,
    ) -> Option<Self> {
        let input_size = input.size();
        let weight_size = weight.size();
        let size = Size {
            width: input_size[3] as usize,
            height: input_size[2] as usize,
        };
        let kernel_size = Size {
            width: weight_size[3] as usize,
            height: weight_size[2] as usize,
        };
        let subsampling = Size {
            width: stride[1] as usize,
            height: stride[0] as usize,
        };
        let padding = Padding {
            top: padding[0] as usize,
            right: padding[1] as usize,
            bottom: padding[0] as usize,
            left: padding[1] as usize,
        };
        let mut packed_bytes = 0;
        let status = unsafe {
            (api.inference)(
                0,
                2,
                input_size[1] as usize,
                weight_size[0] as usize,
                size,
                padding,
                kernel_size,
                subsampling,
                input.data_ptr().cast(),
                weight.data_ptr().cast(),
                bias.data_ptr().cast(),
                output.data_ptr().cast(),
                std::ptr::null_mut(),
                &mut packed_bytes,
                0,
                std::ptr::null(),
                std::ptr::null_mut(),
                std::ptr::null_mut(),
            )
        };
        if status != 0 {
            return None;
        }
        let mut packed = Workspace::new(packed_bytes)?;
        let status = unsafe {
            (api.inference)(
                0,
                2,
                input_size[1] as usize,
                weight_size[0] as usize,
                size,
                padding,
                kernel_size,
                subsampling,
                input.data_ptr().cast(),
                weight.data_ptr().cast(),
                bias.data_ptr().cast(),
                output.data_ptr().cast(),
                packed.pointer(),
                &mut packed.bytes,
                0,
                std::ptr::null(),
                std::ptr::null_mut(),
                std::ptr::null_mut(),
            )
        };
        if status != 0 {
            return None;
        }
        let mut workspace_bytes = 0;
        let status = unsafe {
            (api.inference)(
                0,
                3,
                input_size[1] as usize,
                weight_size[0] as usize,
                size,
                padding,
                kernel_size,
                subsampling,
                input.data_ptr().cast(),
                packed.tensor.data_ptr().cast(),
                bias.data_ptr().cast(),
                output.data_ptr().cast(),
                std::ptr::null_mut(),
                &mut workspace_bytes,
                0,
                std::ptr::null(),
                std::ptr::null_mut(),
                std::ptr::null_mut(),
            )
        };
        if status != 0 {
            return None;
        }
        let workspace = Workspace::new(workspace_bytes)?;
        let threads = std::env::var("CLAIR3_NNPACK_THREADS")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(1);
        let pool = if threads > 1 {
            (unsafe { api.create_pool?(threads) }) as usize
        } else {
            0
        };
        if threads > 1 && pool == 0 {
            return None;
        }
        if std::env::var_os("CLAIR3_PROFILE").is_some() {
            eprintln!(
                "NNPACK reusable {}x{} kernel, packed {} bytes, pool {}",
                weight_size[0], weight_size[1], packed_bytes, threads
            );
        }
        Some(Self {
            packed,
            workspace,
            pool,
            input_size: [input_size[2], input_size[3]],
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn forward(
        &mut self,
        api: &Api,
        layer: &nn::Conv2D,
        stride: [i64; 2],
        padding: [i64; 2],
        input: &Tensor,
        bias: &Tensor,
        output: &Tensor,
    ) -> bool {
        let input_size = input.size();
        let weight_size = layer.ws.size();
        if self.input_size != [input_size[2], input_size[3]] {
            return false;
        }
        let output_size = output.size();
        let input_stride = input_size[1] * input_size[2] * input_size[3];
        let output_stride = output_size[1] * output_size[2] * output_size[3];
        for batch in 0..input_size[0] {
            let status = unsafe {
                (api.inference)(
                    0,
                    3,
                    input_size[1] as usize,
                    weight_size[0] as usize,
                    Size {
                        width: input_size[3] as usize,
                        height: input_size[2] as usize,
                    },
                    Padding {
                        top: padding[0] as usize,
                        right: padding[1] as usize,
                        bottom: padding[0] as usize,
                        left: padding[1] as usize,
                    },
                    Size {
                        width: weight_size[3] as usize,
                        height: weight_size[2] as usize,
                    },
                    Size {
                        width: stride[1] as usize,
                        height: stride[0] as usize,
                    },
                    input
                        .data_ptr()
                        .cast::<f32>()
                        .add((batch * input_stride) as usize),
                    self.packed.tensor.data_ptr().cast(),
                    bias.data_ptr().cast(),
                    output
                        .data_ptr()
                        .cast::<f32>()
                        .add((batch * output_stride) as usize),
                    self.workspace.pointer(),
                    &mut self.workspace.bytes,
                    0,
                    std::ptr::null(),
                    self.pool as *mut c_void,
                    std::ptr::null_mut(),
                )
            };
            if status != 0 {
                return false;
            }
        }
        true
    }
}

pub struct NnpConv2D {
    pub(crate) layer: nn::Conv2D,
    stride: [i64; 2],
    padding: [i64; 2],
    cache: Mutex<Option<Cache>>,
}

impl NnpConv2D {
    pub fn new(layer: nn::Conv2D, stride: [i64; 2], padding: [i64; 2]) -> Self {
        Self {
            layer,
            stride,
            padding,
            cache: Mutex::new(None),
        }
    }

    pub fn invalidate(&self) {
        if let Ok(mut cache) = self.cache.lock() {
            *cache = None;
        }
    }

    fn forward_spatial(
        &self,
        api: &Api,
        input: &Tensor,
        weight: &Tensor,
        bias: &Tensor,
    ) -> Option<Tensor> {
        let input_size = input.size();
        let weight_size = weight.size();
        let height = (input_size[2] + 2 * self.padding[0] - weight_size[2]) / self.stride[0] + 1;
        let width = (input_size[3] + 2 * self.padding[1] - weight_size[3]) / self.stride[1] + 1;
        if height <= 0
            || width <= 0
            || self.stride[0] != self.stride[1]
            || self.padding[0] != self.padding[1]
        {
            return None;
        }
        let group_size = 8;
        let batch = input_size[0];
        let groups = batch / group_size;
        if groups == 0 {
            return None;
        }
        let border = self.stride[1] + 1;
        let vertical = self.stride == [1, 1] && weight_size[0] != 128;
        let input_pitch = input_size[3] + border;
        let input_pitch_y = input_size[2] + border;
        let output_pitch =
            (input_pitch + 2 * self.padding[1] - weight_size[3]) / self.stride[1] + 1;
        let output_pitch_y =
            (input_pitch_y + 2 * self.padding[0] - weight_size[2]) / self.stride[0] + 1;
        let result = Tensor::empty(
            [batch, weight_size[0], height, width],
            (Kind::Float, Device::Cpu),
        );
        let packed_input = if vertical {
            Tensor::zeros(
                [1, input_size[1], group_size * input_pitch_y, input_size[3]],
                (Kind::Float, Device::Cpu),
            )
        } else {
            Tensor::zeros(
                [1, input_size[1], input_size[2], group_size * input_pitch],
                (Kind::Float, Device::Cpu),
            )
        };
        let packed_output = if vertical {
            Tensor::empty(
                [1, weight_size[0], group_size * output_pitch_y, width],
                (Kind::Float, Device::Cpu),
            )
        } else {
            Tensor::empty(
                [1, weight_size[0], height, group_size * output_pitch],
                (Kind::Float, Device::Cpu),
            )
        };
        let mut guard = self.cache.lock().ok()?;
        if guard.is_none() {
            *guard = Cache::new(
                api,
                self.stride,
                self.padding,
                &packed_input,
                weight,
                bias,
                &packed_output,
            );
        }
        for first in (0..groups * group_size).step_by(group_size as usize) {
            for sample in 0..group_size {
                if vertical {
                    packed_input
                        .narrow(2, sample * input_pitch_y, input_size[2])
                        .copy_(&input.narrow(0, first + sample, 1));
                } else {
                    packed_input
                        .narrow(3, sample * input_pitch, input_size[3])
                        .copy_(&input.narrow(0, first + sample, 1));
                }
            }
            if !guard.as_mut()?.forward(
                api,
                &self.layer,
                self.stride,
                self.padding,
                &packed_input,
                bias,
                &packed_output,
            ) {
                return None;
            }
            for sample in 0..group_size {
                if vertical {
                    result
                        .narrow(0, first + sample, 1)
                        .copy_(&packed_output.narrow(2, sample * output_pitch_y, height));
                } else {
                    result
                        .narrow(0, first + sample, 1)
                        .copy_(&packed_output.narrow(3, sample * output_pitch, width));
                }
            }
        }
        if groups * group_size < batch {
            let tail = nn::Module::forward(
                &self.layer,
                &input.narrow(0, groups * group_size, batch - groups * group_size),
            );
            result
                .narrow(0, groups * group_size, batch - groups * group_size)
                .copy_(&tail);
        }
        Some(result)
    }

    pub fn forward(&self, input: &Tensor) -> Option<Tensor> {
        if input.device() != Device::Cpu
            || input.kind() != Kind::Float
            || self.layer.ws.device() != Device::Cpu
            || self.layer.ws.kind() != Kind::Float
        {
            return None;
        }
        let api = api()?;
        let input = input.contiguous();
        let weight = self.layer.ws.contiguous();
        let bias = match &self.layer.bs {
            Some(bias) => bias.contiguous(),
            None => Tensor::zeros([weight.size()[0]], (Kind::Float, Device::Cpu)),
        };
        let input_size = input.size();
        let weight_size = weight.size();
        if input_size.len() != 4
            || weight_size.len() != 4
            || weight_size[1] != input_size[1]
            || weight_size[2..] != [3, 3]
        {
            return None;
        }
        self.forward_spatial(api, &input, &weight, &bias)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tch::nn::Module;

    #[test]
    fn cached_convolution_matches_the_standard_backend() {
        let store = nn::VarStore::new(Device::Cpu);
        for stride in [[2, 2], [1, 1]] {
            let layer =
                super::super::layers::conv2d(&store.root(), 8, 64, [3, 3], stride, [1, 1], 1, true);
            let cached = NnpConv2D::new(layer, stride, [1, 1]);
            for batch in [16, 32, 35] {
                let input = Tensor::randn([batch, 8, 55, 33], (Kind::Float, Device::Cpu));
                let expected = cached.layer.forward(&input);
                let Some(actual) = cached.forward(&input) else {
                    return;
                };
                let expected = Vec::<f32>::try_from(&expected.flatten(0, -1)).unwrap();
                let actual = Vec::<f32>::try_from(&actual.flatten(0, -1)).unwrap();
                let maximum = expected
                    .iter()
                    .zip(&actual)
                    .map(|(left, right)| (left - right).abs())
                    .fold(0.0_f32, f32::max);
                let unequal = expected
                    .iter()
                    .zip(&actual)
                    .filter(|(left, right)| left != right)
                    .count();
                eprintln!("stride {stride:?}, batch {batch}: unequal {unequal}, maximum {maximum}");
                let tolerance = if stride == [1, 1] { 0.0002 } else { 0.0 };
                assert!(maximum <= tolerance, "maximum error {maximum}");
            }
        }
    }
}
