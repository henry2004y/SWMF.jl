# SWMF.jl Documentation

## Overview

!!! note
    This package is still under development, so be careful for any future breaking changes!

SWMF data reader and converter in Julia.

This package is the inherited from its predecessor in IDL (developed by G.Toth) and Matlab (developed by H.Zhou), and was originally part of [VisAna](https://github.com/henry2004y/VisAnaJulia).
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

Data loading speed of a 2.4GB 3D binary file on Macbook Pro with quad core 2.2 GHz Intel i7 and 16 GB 1600 MHz DDR3:

| Language |  tmin |   tmax |  tmean |
|:-------|:-----:|:------:|:------:|
| Julia  | 1.44s | 3.77s  |  2.07s |
| IDL    | 6.00s | 6.18s  |  6.08s |
| MATLAB | 6.67s | 16.02s | 10.60s |

## Developers

VisAna is developed by [Hongyang Zhou](https://github.com/henry2004y).
