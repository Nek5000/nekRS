# GEMV

`GEMV` is a library for performing platform independent dense matrix-vector
products on the GPUs.

## Build Instructions

This project uses conda to manage dependencies (CMake, clang-format, clang-tidy
and other dependencies for documentation). Dependencies can be installed by
executing following commands after installing
[conda](https://docs.conda.io/en/latest/miniconda.html).
```sh
conda env create -f environment-dev.yml
conda activate gemv
```

Then simply run `gemv.sh` script to build and install the library.
```sh
./gemv.sh --install
```

You can format the source code with `clang-format` using the option `--format`.
```
./gemv.sh --format
```

Use `--help` to see all the options supported by `gemv.sh` script.
```
./gemv.sh --help
```
