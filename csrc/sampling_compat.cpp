#include <cstdint>
#include <cstdlib>
#include <new>

#if defined(__linux__)
struct GlibcState {
    random_data data{};
    char state[128]{};

    explicit GlibcState(unsigned int seed) {
        initstate_r(seed, state, sizeof(state), &data);
    }
};
#endif

extern "C" void *clair3_sampling_new(unsigned int seed) {
#if defined(__linux__)
    return new (std::nothrow) GlibcState(seed);
#else
    return new (std::nothrow) unsigned int(seed);
#endif
}

extern "C" int32_t clair3_sampling_next(void *state) {
#if defined(__linux__)
    auto *generator = static_cast<GlibcState *>(state);
    int32_t result = 0;
    random_r(&generator->data, &result);
    return result;
#else
    return rand_r(static_cast<unsigned int *>(state));
#endif
}

extern "C" void clair3_sampling_free(void *state) {
#if defined(__linux__)
    delete static_cast<GlibcState *>(state);
#else
    delete static_cast<unsigned int *>(state);
#endif
}
