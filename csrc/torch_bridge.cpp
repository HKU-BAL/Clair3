#include <ATen/ATen.h>
#include <c10/core/GradMode.h>
#include <c10/util/intrusive_ptr.h>
#include <torch/csrc/distributed/c10d/Backend.hpp>
#include <torch/csrc/distributed/c10d/TCPStore.hpp>
#ifdef USE_C10D_GLOO
#include <torch/csrc/distributed/c10d/ProcessGroupGloo.hpp>
#endif
#ifdef USE_C10D_NCCL
#include <c10/cuda/CUDAFunctions.h>
#include <torch/csrc/distributed/c10d/ProcessGroupNCCL.hpp>
#endif

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {
using BackendHandle = c10::intrusive_ptr<c10d::Backend>;

void set_error(char* buffer, size_t size, const std::exception& error) {
    if (buffer != nullptr && size > 0) {
        std::snprintf(buffer, size, "%s", error.what());
    }
}

std::string env_or_default(const char* name, const char* fallback) {
    const char* value = std::getenv(name);
    return value == nullptr || *value == '\0' ? fallback : value;
}

template <typename F>
int guarded(char* error, size_t size, F&& action) {
    try {
        action();
        return 0;
    } catch (const std::exception& exception) {
        set_error(error, size, exception);
        return 1;
    }
}
} // namespace

extern "C" int clair3_pg_init(
    void** output,
    int rank,
    int world_size,
    int local_rank,
    int use_cuda,
    char* error,
    size_t error_size) {
    return guarded(error, error_size, [&] {
        c10d::TCPStoreOptions options;
        options.port = static_cast<std::uint16_t>(std::stoi(env_or_default("MASTER_PORT", "29500")));
        options.isServer = rank == 0;
        options.numWorkers = static_cast<std::size_t>(world_size);
        options.timeout = std::chrono::seconds(1800);
        auto store = c10::make_intrusive<c10d::TCPStore>(env_or_default("MASTER_ADDR", "127.0.0.1"), options);
#ifdef CLAIR3_WITH_CUDA
        if (use_cuda != 0) {
            c10::cuda::set_device(static_cast<c10::DeviceIndex>(local_rank));
            *output = new BackendHandle(c10::make_intrusive<c10d::ProcessGroupNCCL>(store, rank, world_size));
            return;
        }
#else
        (void)local_rank;
#endif
        if (use_cuda != 0) {
            throw std::runtime_error("this Clair3 build does not include CUDA distributed support");
        }
        auto gloo_options = c10d::ProcessGroupGloo::Options::create();
        gloo_options->devices.push_back(c10d::ProcessGroupGloo::createDefaultDevice());
        auto backend = c10::make_intrusive<c10d::ProcessGroupGloo>(store, rank, world_size, gloo_options);
        backend->setSequenceNumberForGroup();
        *output = new BackendHandle(std::move(backend));
    });
}

extern "C" int clair3_pg_allreduce(
    void* handle,
    void* tensor,
    int world_size,
    char* error,
    size_t error_size) {
    return guarded(error, error_size, [&] {
        auto* backend = static_cast<BackendHandle*>(handle);
        auto* value = static_cast<at::Tensor*>(tensor);
        c10::NoGradGuard no_grad;
        std::vector<at::Tensor> tensors{at::empty_like(*value).copy_(*value)};
        (*backend)->allreduce(tensors)->wait();
        *value = tensors.front() / world_size;
    });
}

extern "C" int clair3_pg_broadcast(
    void* handle,
    long long value,
    long long* output,
    int root,
    int local_rank,
    int use_cuda,
    char* error,
    size_t error_size) {
    return guarded(error, error_size, [&] {
        auto* backend = static_cast<BackendHandle*>(handle);
        at::Device device = at::Device(at::kCPU);
#ifdef CLAIR3_WITH_CUDA
        if (use_cuda != 0) {
            device = at::Device(at::kCUDA, static_cast<c10::DeviceIndex>(local_rank));
        }
#else
        (void)local_rank;
#endif
        at::Tensor tensor = at::zeros({1}, at::TensorOptions().dtype(at::kLong).device(device));
        tensor.fill_(static_cast<int64_t>(value));
        std::vector<at::Tensor> tensors{tensor};
        c10d::BroadcastOptions options;
        options.rootRank = root;
        (*backend)->broadcast(tensors, options)->wait();
        *output = tensors.front().item<int64_t>();
    });
}

extern "C" int clair3_pg_barrier(void* handle, char* error, size_t error_size) {
    return guarded(error, error_size, [&] {
        auto* backend = static_cast<BackendHandle*>(handle);
        (*backend)->barrier()->wait();
    });
}

extern "C" void clair3_pg_free(void* handle) {
    delete static_cast<BackendHandle*>(handle);
}
