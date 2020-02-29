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

|        | tmin [s] | tmax [s] | tmean [s] |
|--------|--------|--------|-----------|
| Julia  | 1.44   | 3.77   | 2.07      |
| IDL    | 6.00   | 6.18   | 6.08      |
| MATLAB | 6.67   | 16.02  | 10.60     |

## Developers

VisAna is developed by [Hongyang Zhou](https://github.com/henry2004y).
