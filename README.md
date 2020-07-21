# SWMF
[![](https://travis-ci.com/henry2004y/SWMF.jl.svg?branch=master)][travis-url]
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![](https://img.shields.io/badge/docs-latest-blue.svg)][SWMF-doc]
[![][codecov-img]][codecov-url]

Fast SWMF data reader and converter using Julia.

This is originally part of the [VisAna.jl](https://github.com/henry2004y/VisAnaJulia) package.
Later this was moved out and became a stand-alone package.
Note that you must have Julia and this package installed to successfully use it in Python.

This package provides the following functionalities:
  * simulation data reader
  * 2D/3D region cut from the whole data
  * data format conversion from Tecplot to VTK
  * programming language interoperability

For more details, please check the [document][SWMF-doc].

## Prerequisites

Julia 1.0+

## Installation
```
using Pkg
Pkg.add(PackageSpec(url="https://github.com/henry2004y/SWMF", rev="master"))
```

## Usage
```
#using Pkg; Pkg.activate(".") # for dev only
using SWMF
```

See the [examples](https://henry2004y.github.io/SWMF/dev/man/examples/).

### Using from Python

In Python, you can easily take advantage of this package with the aid of [PyJulia](https://pyjulia.readthedocs.io/en/latest/).
After the installation, in the Python repl:
```python
from julia import SWMF
dir = 'test'
filename = '1d__raw_2_t25.60000_n00000258.out'
data = SWMF.readdata(filename, dir=dir)
```
There you have it! Enjoy!

## Benchmark

See the [benchmark](https://henry2004y.github.io/SWMF/dev/#Benchmark-1) in the document.

## Author

* **Hongyang Zhou** - *Initial work* - [henry2004y](https://github.com/henry2004y)

## Acknowledgments

* All the nice guys who share their codes!

[travis-url]: https://travis-ci.com/github/henry2004y/SWMF.jl/builds
[codecov-img]: https://codecov.io/gh/henry2004y/SWMF/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/henry2004y/SWMF
[SWMF-doc]: https://henry2004y.github.io/SWMF.jl/dev
