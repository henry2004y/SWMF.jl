# SWMF.jl Documentation

## Overview

!!! note
    This package is still under development, so be careful for any future breaking changes!

SWMF data reader and converter in Julia.

This package inherits the ideas and code structures from its predecessor in IDL (developed by G.Toth) and Matlab (developed by H.Zhou), and was originally part of [VisAna](https://github.com/henry2004y/VisAnaJulia).
It can be combined with the VTK format converter [writeVTK](https://github.com/jipolanco/WriteVTK.jl) to generate files for Paraview and Tecplot.
By default the file size will be reduced with compression level 6, but the actual compression ratio depends on the original data.

This package provides the following functionalities:
  * simulation data reader
  * data format conversion
  * programming language interoperability

The ultimate goal is to build a convenient tool of reading and analyzing simulation outputs which is easy to install and easy to use.

!!! tip "Ready to use?"
    Feel free to contact the author for any help or collaboration!

## Installation
Install VisAna from the `julia REPL` prompt with
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/henry2004y/SWMF", rev="master"))
```

```@contents
Pages = [
    "man/guide.md",
    "man/examples.md",
    "man/functions.md",
    "man/types.md"
]
Depth = 1
```

## Benchmark

Data loading speed of a 2.4GB 3D binary file and 65KB 2D binary file on Macbook Pro with quad core 2.2 GHz Intel i7 and 16 GB 1600 MHz DDR3:

| Language |   tmax |  tmean |
|:-------|:------:|:------:|
| Julia  | 2.73s  |  1.32s |
| Python | 1.35s  |  1.34s |
| IDL    | 6.18s  |  6.08s |
| MATLAB | 16.02s | 10.60s |

| 65KB   | tmean [μs] |
|:-------|:---------:|
| Julia  | 163.36    |
| Python | 4390.95   |
| IDL    | 1970.29   |
| MATLAB | 19273.25  |

The Julia, IDL, and MATLAB version all shares the same kernel design. The timings are obtained for Julia v1.3.1, Python 3.7.6 + Numpy 1.18.1, IDL 8.5, and MATLAB R2018b.
For dynamic languages, the first time when function gets executed is usually also the slowest. Currently [spacepy](https://github.com/spacepy/spacepy) performs slightly better because of the well-optimized numpy library in C. For small data sizes, Julia is much faster than others.

## Developers

VisAna is developed by [Hongyang Zhou](https://github.com/henry2004y).
