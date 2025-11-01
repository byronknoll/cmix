project('cmix', 'cpp',
  version: '21',
  default_options: [
    'cpp_std=c++14',
    'warning_level=3',
    'buildtype=release'
  ]
)

# compiler setup
cpp = meson.get_compiler('cpp')

# base compiler flags
base_args = []

# test and add optimization flags
release_args = []

# optimization flags:
# newer clang deprecates -Ofast in favor of '-O3 -ffast-math'
# Try in order of preference:: -Ofast, then -O3 + -ffast-math, then just -O3
if cpp.has_argument('-Ofast')
  release_args += ['-Ofast']
else
  # -Ofast not supported or deprecated, try the equivalent flags
  if cpp.has_argument('-O3')
    release_args += ['-O3']
  endif
  
  # -ffast-math enables non-conforming floating point optimizations
  # (same behavior as -Ofast over -O3)
  if cpp.has_argument('-ffast-math')
    release_args += ['-ffast-math']
  endif
endif

# -march=native optimizes for the host CPU architecture
# only makes sense for native builds, not cross-compilation
# for cross-compilation, specify target arch in cross-file
if not meson.is_cross_build()
  if cpp.has_argument('-march=native')
    release_args += ['-march=native']
  endif
else
  warning('Cross-compilation detected: skipping -march=native (use cross-file to set target architecture)')
endif

cmix_sources = files(
  'src/coder/decoder.cpp',
  'src/coder/encoder.cpp',
  'src/context-manager.cpp',
  'src/contexts/bit-context.cpp',
  'src/contexts/bracket-context.cpp',
  'src/contexts/combined-context.cpp',
  'src/contexts/context-hash.cpp',
  'src/contexts/indirect-hash.cpp',
  'src/contexts/interval-hash.cpp',
  'src/contexts/interval.cpp',
  'src/contexts/sparse.cpp',
  'src/mixer/byte-mixer.cpp',
  'src/mixer/lstm-layer.cpp',
  'src/mixer/lstm.cpp',
  'src/mixer/mixer-input.cpp',
  'src/mixer/mixer.cpp',
  'src/mixer/sigmoid.cpp',
  'src/mixer/sse.cpp',
  'src/models/bracket.cpp',
  'src/models/byte-model.cpp',
  'src/models/direct-hash.cpp',
  'src/models/direct.cpp',
  'src/models/indirect.cpp',
  'src/models/fxcmv1.cpp',
  'src/models/match.cpp',
  'src/models/paq8.cpp',
  'src/models/ppmd.cpp',
  'src/predictor.cpp',
  'src/preprocess/dictionary.cpp',
  'src/preprocess/preprocessor.cpp',
  'src/runner.cpp',
  'src/states/nonstationary.cpp',
  'src/states/run-map.cpp'
)

enwik9_preproc_sources = files(
  'src/enwik9-preproc/main.cpp'
)

remap_sources = files(
  'src/enwik9-preproc/article_remap.cpp'
)

if get_option('buildtype') == 'debug'
  compile_args = base_args
else
  compile_args = base_args + release_args
endif

# build executables
executable('cmix',
  cmix_sources,
  cpp_args: compile_args,
  install: true
)

executable('enwik9-preproc',
  enwik9_preproc_sources,
  cpp_args: compile_args,
  install: true
)

executable('remap',
  remap_sources,
  cpp_args: compile_args,
  install: true
)
