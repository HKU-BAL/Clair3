use std::{env, fs, path::PathBuf, process::Command};

fn python_details(python: &str) -> (Vec<PathBuf>, Vec<PathBuf>, String) {
    let output = Command::new(python)
        .arg("-c")
        .arg(
            "import torch\nfrom torch.utils import cpp_extension\nprint('ABI:', torch._C._GLIBCXX_USE_CXX11_ABI)\nfor p in cpp_extension.include_paths(): print('INCLUDE:', p)\nfor p in cpp_extension.library_paths(): print('LIB:', p)",
        )
        .output()
        .expect("failed to locate LibTorch through Python");
    assert!(
        output.status.success(),
        "failed to import LibTorch through {python}"
    );
    let output = String::from_utf8(output.stdout).expect("non-UTF-8 LibTorch configuration");
    let mut includes = Vec::new();
    let mut libraries = Vec::new();
    let mut abi = "1".to_string();
    for line in output.lines() {
        if let Some(path) = line.strip_prefix("INCLUDE: ") {
            includes.push(PathBuf::from(path));
        } else if let Some(path) = line.strip_prefix("LIB: ") {
            libraries.push(PathBuf::from(path));
        } else if let Some(value) = line.strip_prefix("ABI: ") {
            abi = if matches!(value, "1" | "True" | "true") {
                "1".to_string()
            } else {
                "0".to_string()
            };
        }
    }
    (includes, libraries, abi)
}

fn execute(command: &mut Command) {
    let description = format!("{command:?}");
    let status = command
        .status()
        .unwrap_or_else(|error| panic!("failed to execute {description}: {error}"));
    assert!(status.success(), "command failed: {description}");
}

fn prepare_gloo() -> (PathBuf, PathBuf) {
    if let Ok(library) = env::var("GLOO_LIBRARY") {
        let library = PathBuf::from(library);
        let include = env::var("GLOO_INCLUDE_DIR")
            .map(PathBuf::from)
            .unwrap_or_else(|_| PathBuf::from("target/gloo-src"));
        assert!(
            library.is_file(),
            "Gloo library not found: {}",
            library.display()
        );
        assert!(
            include.join("gloo/algorithm.h").is_file(),
            "Gloo headers not found in {}",
            include.display()
        );
        return (include, library);
    }
    let out = PathBuf::from(env::var_os("OUT_DIR").expect("OUT_DIR"));
    let manifest = PathBuf::from(env::var_os("CARGO_MANIFEST_DIR").expect("manifest directory"));
    let vendored = manifest.join("third_party/gloo");
    let source = if vendored.join("gloo/algorithm.h").is_file() {
        vendored
    } else {
        out.join("gloo-src")
    };
    if !source.join("gloo/algorithm.h").is_file() {
        let archive = out.join("gloo-master.tar.gz");
        if !archive.is_file() {
            execute(
                Command::new("curl")
                    .arg("-fLsS")
                    .arg("--retry")
                    .arg("3")
                    .arg("-o")
                    .arg(&archive)
                    .arg("https://github.com/pytorch/gloo/archive/refs/heads/master.tar.gz"),
            );
        }
        execute(
            Command::new("tar")
                .arg("-xzf")
                .arg(&archive)
                .arg("-C")
                .arg(&out),
        );
        let extracted = out.join("gloo-master");
        if extracted.is_dir() {
            fs::rename(extracted, &source).expect("failed to prepare Gloo sources");
        }
    }
    let build = out.join("gloo-build");
    let library = if build.join("gloo/libgloo.a").is_file() {
        build.join("gloo/libgloo.a")
    } else {
        build.join("libgloo.a")
    };
    if !library.is_file() {
        execute(
            Command::new("cmake")
                .arg("-S")
                .arg(&source)
                .arg("-B")
                .arg(&build)
                .arg("-DBUILD_TEST=OFF")
                .arg("-DBUILD_BENCHMARK=OFF")
                .arg("-DGLOO_STATIC=ON")
                .arg("-DUSE_LIBUV=OFF")
                .arg("-DUSE_MPI=OFF")
                .arg("-DUSE_CUDA=OFF")
                .arg("-DCMAKE_POSITION_INDEPENDENT_CODE=ON")
                .arg("-DCMAKE_BUILD_TYPE=Release"),
        );
        execute(
            Command::new("cmake")
                .arg("--build")
                .arg(&build)
                .arg("--target")
                .arg("gloo")
                .arg("--parallel")
                .arg(env::var("NUM_JOBS").unwrap_or_else(|_| "2".to_string())),
        );
    }
    let library = if build.join("gloo/libgloo.a").is_file() {
        build.join("gloo/libgloo.a")
    } else {
        build.join("libgloo.a")
    };
    assert!(library.is_file(), "failed to build Gloo");
    (source, library)
}

fn downloaded_libtorch() -> Option<PathBuf> {
    let out_dir = PathBuf::from(env::var_os("OUT_DIR")?);
    let build_dir = out_dir.parent()?.parent()?;
    let mut candidates = Vec::new();
    for entry in fs::read_dir(build_dir).ok()? {
        let entry = entry.ok()?;
        if !entry
            .file_name()
            .to_string_lossy()
            .starts_with("torch-sys-")
        {
            continue;
        }
        let root = entry.path().join("out/libtorch/libtorch");
        if root.join("include/ATen/ATen.h").is_file() && root.join("lib").is_dir() {
            candidates.push(root);
        }
    }
    candidates.sort_by_key(|root| {
        fs::metadata(root)
            .and_then(|metadata| metadata.modified())
            .unwrap_or(std::time::SystemTime::UNIX_EPOCH)
    });
    candidates.pop()
}

fn prepare_phase_core() -> Option<PathBuf> {
    if env::var("CARGO_CFG_TARGET_OS").is_ok_and(|os| os == "windows") {
        return None;
    }
    let manifest = PathBuf::from(env::var_os("CARGO_MANIFEST_DIR")?);
    let root = manifest.join("third_party/whatshap");
    let python = env::var("PYTHON").unwrap_or_else(|_| "python3".to_string());
    let output = Command::new(python)
        .arg("-c")
        .arg("import sysconfig; print(sysconfig.get_paths()['include']); print(sysconfig.get_config_var('EXT_SUFFIX'))")
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let output = String::from_utf8(output.stdout).ok()?;
    let mut lines = output.lines();
    let include = PathBuf::from(lines.next()?);
    if !include.join("Python.h").is_file() {
        return None;
    }
    let suffix = lines.next().unwrap_or(".so");
    let sources = [
        "whatshap/core.cpp",
        "src/pedigree.cpp",
        "src/pedigreedptable.cpp",
        "src/pedigreecolumncostcomputer.cpp",
        "src/columnindexingiterator.cpp",
        "src/columnindexingscheme.cpp",
        "src/entry.cpp",
        "src/graycodes.cpp",
        "src/read.cpp",
        "src/readset.cpp",
        "src/columniterator.cpp",
        "src/indexset.cpp",
        "src/genotype.cpp",
        "src/binomial.cpp",
        "src/pedmecheuristic.cpp",
        "src/multinomial.cpp",
        "src/pedigreepartitions.cpp",
        "src/phredgenotypelikelihoods.cpp",
        "src/genotyper.cpp",
        "src/genotypedistribution.cpp",
        "src/genotypedptable.cpp",
        "src/genotypecolumncostcomputer.cpp",
        "src/backwardcolumniterator.cpp",
        "src/transitionprobabilitycomputer.cpp",
        "src/hapchat/basictypes.cpp",
        "src/hapchat/balancedcombinations.cpp",
        "src/hapchat/binomialcoefficient.cpp",
        "src/hapchat/hapchatcore.cpp",
        "src/hapchat/hapchatcolumniterator.cpp",
        "src/polyphase/allelematrix.cpp",
        "src/polyphase/clustereditingsolution.cpp",
        "src/polyphase/clustereditingsolver.cpp",
        "src/polyphase/edgeheap.cpp",
        "src/polyphase/inducedcostheuristic.cpp",
        "src/polyphase/staticsparsegraph.cpp",
        "src/polyphase/trianglesparsematrix.cpp",
        "src/polyphase/readscoring.cpp",
        "src/polyphase/haplothreader.cpp",
        "src/polyphase/tupleconverter.cpp",
        "src/caller.cpp",
    ];
    let compiler = cc::Build::new().cpp(true).get_compiler();
    let mut command = Command::new(compiler.path());
    command
        .args(compiler.args())
        .arg("-std=c++11")
        .arg("-O3")
        .arg("-fPIC")
        .arg("-shared")
        .arg("-I")
        .arg(root.join("src"))
        .arg("-I")
        .arg(include);
    if env::var("CARGO_CFG_TARGET_OS").is_ok_and(|os| os == "macos") {
        command.arg("-undefined").arg("dynamic_lookup");
    }
    command.args(sources.map(|source| root.join(source)));
    let output = PathBuf::from(env::var_os("OUT_DIR")?).join(format!("whatshap_core{suffix}"));
    command.arg("-o").arg(&output);
    execute(&mut command);
    Some(output)
}

fn prepare_readselector() -> Option<PathBuf> {
    if env::var("CARGO_CFG_TARGET_OS").is_ok_and(|os| os == "windows") {
        return None;
    }
    let manifest = PathBuf::from(env::var_os("CARGO_MANIFEST_DIR")?);
    let root = manifest.join("third_party/whatshap");
    let python = env::var("PYTHON").unwrap_or_else(|_| "python3".to_string());
    let output = Command::new(python)
        .arg("-c")
        .arg("import sysconfig; print(sysconfig.get_paths()['include']); print(sysconfig.get_config_var('EXT_SUFFIX'))")
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let output = String::from_utf8(output.stdout).ok()?;
    let mut lines = output.lines();
    let include = PathBuf::from(lines.next()?);
    if !include.join("Python.h").is_file() {
        return None;
    }
    let suffix = lines.next().unwrap_or(".so");
    let compiler = cc::Build::new().cpp(true).get_compiler();
    let mut command = Command::new(compiler.path());
    command
        .args(compiler.args())
        .arg("-std=c++11")
        .arg("-O3")
        .arg("-fPIC")
        .arg("-shared")
        .arg("-I")
        .arg(root.join("src"))
        .arg("-I")
        .arg(include);
    if env::var("CARGO_CFG_TARGET_OS").is_ok_and(|os| os == "macos") {
        command.arg("-undefined").arg("dynamic_lookup");
    }
    let output = PathBuf::from(env::var_os("OUT_DIR")?).join(format!("clair3_readselect{suffix}"));
    command
        .arg(root.join("fast_readselect.cpp"))
        .arg("-o")
        .arg(&output);
    execute(&mut command);
    Some(output)
}

fn main() {
    println!("cargo:rerun-if-changed=csrc/torch_bridge.cpp");
    println!("cargo:rerun-if-changed=third_party/whatshap");
    println!("cargo:rerun-if-changed=csrc/sampling_compat.cpp");
    for name in [
        "LIBTORCH",
        "LIBTORCH_USE_PYTORCH",
        "LIBTORCH_CXX11_ABI",
        "PYTHON",
        "DEP_TCH_LIBTORCH_LIB",
        "GLOO_LIBRARY",
        "GLOO_INCLUDE_DIR",
    ] {
        println!("cargo:rerun-if-env-changed={name}");
    }
    let phase_core = prepare_phase_core().unwrap_or_else(|| {
        let path = PathBuf::from(env::var_os("OUT_DIR").expect("OUT_DIR")).join("whatshap_core.so");
        fs::write(&path, []).expect("failed to prepare the fallback phase-core marker");
        path
    });
    println!("cargo:rustc-env=CLAIR3_PHASE_CORE={}", phase_core.display());
    let read_selector = prepare_readselector().unwrap_or_else(|| {
        let path =
            PathBuf::from(env::var_os("OUT_DIR").expect("OUT_DIR")).join("clair3_readselect.so");
        fs::write(&path, b"").expect("failed to prepare the fallback selector marker");
        path
    });
    println!(
        "cargo:rustc-env=CLAIR3_READSELECT={}",
        read_selector.display()
    );
    if env::var_os("CARGO_FEATURE_DISTRIBUTED").is_none() {
        return;
    }
    let mut includes = Vec::new();
    let mut libraries = Vec::new();
    let mut abi = "1".to_string();
    if let Ok(root) = env::var("LIBTORCH") {
        let root = PathBuf::from(root);
        includes.push(root.join("include"));
        includes.push(root.join("include/torch/csrc/api/include"));
        libraries.push(root.join("lib"));
    } else if env::var_os("LIBTORCH_USE_PYTORCH").is_some() {
        let python = env::var("PYTHON").unwrap_or_else(|_| "python3".to_string());
        let details = python_details(&python);
        includes = details.0;
        libraries = details.1;
        abi = details.2;
    }
    let abi = env::var("LIBTORCH_CXX11_ABI").unwrap_or(abi);
    assert!(
        matches!(abi.as_str(), "0" | "1"),
        "LIBTORCH_CXX11_ABI must be 0 or 1"
    );
    if let Some(path) = env::var_os("DEP_TCH_LIBTORCH_LIB") {
        let path = PathBuf::from(path);
        libraries.push(path.clone());
        if let Some(root) = path.parent() {
            includes.push(root.join("include"));
            includes.push(root.join("include/torch/csrc/api/include"));
        }
    }
    if includes.is_empty()
        && libraries.is_empty()
        && let Some(root) = downloaded_libtorch()
    {
        includes.push(root.join("include"));
        includes.push(root.join("include/torch/csrc/api/include"));
        libraries.push(root.join("lib"));
    }
    let (gloo_include, gloo_library) = prepare_gloo();
    includes.push(gloo_include);
    libraries.push(gloo_library.parent().unwrap().to_path_buf());
    includes.retain(|path| path.is_dir());
    includes.sort();
    includes.dedup();
    libraries.retain(|path| path.is_dir());
    libraries.sort();
    libraries.dedup();
    assert!(
        includes
            .iter()
            .any(|path| path.join("ATen/ATen.h").is_file()),
        "LibTorch ATen headers not found; set LIBTORCH to a complete distribution"
    );
    assert!(
        libraries.iter().any(|path| {
            path.join("libtorch.so").is_file()
                || path.join("libtorch.dylib").is_file()
                || path.join("torch.lib").is_file()
        }),
        "LibTorch libraries not found; set LIBTORCH to a complete distribution"
    );
    let cuda = env::var_os("CARGO_FEATURE_CUDA").is_some();
    let mut build = cc::Build::new();
    build
        .cpp(true)
        .pic(true)
        .std("c++20")
        .warnings(false)
        .includes(&includes)
        .define("USE_C10D_GLOO", "1")
        .define("GLOG_USE_GLOG_EXPORT", None)
        .define("_GLIBCXX_USE_CXX11_ABI", abi.as_str())
        .file("csrc/torch_bridge.cpp")
        .file("csrc/sampling_compat.cpp");
    if cuda {
        build.define("USE_C10D_NCCL", "1");
        build.define("CLAIR3_WITH_CUDA", "1");
    }
    build.compile("clair3_torch_bridge");
    for directory in libraries {
        println!("cargo:rustc-link-search=native={}", directory.display());
        println!("cargo:rustc-link-arg=-Wl,-rpath,{}", directory.display());
    }
    println!("cargo:rustc-link-lib=static=gloo");
    println!("cargo:rustc-link-lib=dylib=torch");
    println!("cargo:rustc-link-lib=dylib=torch_cpu");
    println!("cargo:rustc-link-lib=dylib=c10");
    if cuda {
        println!("cargo:rustc-link-lib=dylib=torch_cuda");
        println!("cargo:rustc-link-lib=dylib=c10_cuda");
    }
    if cfg!(target_os = "macos") {
        println!("cargo:rustc-link-arg=-Wl,-rpath,@loader_path/../lib/clair3");
    } else if cfg!(target_os = "linux") {
        println!("cargo:rustc-link-arg=-Wl,-rpath,$ORIGIN/../lib/clair3");
    }
}
