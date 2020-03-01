# Guide

Demos are provided for calling Matlab/Python directly from Julia for debugging and testing. This part will later be separated out for potential Python and Matlab users. Currently the plotting and interpolation needed during plotting are done in Python. For instance, the 3D scatterred interpolation is done via `Interpolate` in Scipy. Hopefully these additional dependencies will be cut down.

The VTK files does not have timestep information. To allow for further time series processing in Paraview, a script `create_pvd.jl` is provided for generating the pvd container.

In principle, I could also try some multi-block (VTM) type for data conversion.

## Tricks

- This is the first time I use Julia for reading general ascii/binary files. It was a pain at first due to the lack of examples and documents using any basic function like read/read!, but fortunately I figured them out myself. One trick in reading binary array data is the usage of view, or subarrays, in Julia. In order to achieve that, I have to implement my own `read!` function in addition to the base ones.
- Tecplot and VTK unstructured data formats have the same connectivity ordering for hexahedron, but different ordering for voxel (in VTK). A function `swaprows` is implemented to switch the orderings.
- Because of the embarrassing parallelism nature of postprocessing, it is quite easy to take advantage of parallel approaches to process the data.


## Support for more complicated grid structures

For the plotting, streamline tracing and particle tracing, a common problem is the grid and related interpolation process. I am envisioning a more general approach to deal with block-based and unstructured grid to provide fundamental support for all of these.

##

A real open-source project is a collaborated work not only from a bunch of people, but also a group of languages. In Julia, this can be achieved with the help of the [Package manager](https://julialang.github.io/Pkg.jl/dev/).

I want to have some C dependencies in my code instead of rewriting everything in Julia. This would serve as an attempt to quickly make things work.

Right now this seems to be a little bit difficult for me. I need to learn from experts. The tracing scheme in C is rewritten in Julia so I don't need to bother for now.
Checkout [BinaryBuilder](https://juliapackaging.github.io/BinaryBuilder.jl/latest/#Project-flow-1) for more information.


## Issues

I am enlightened by the way spacepy handles files. Instead of a separate header and data array, it may be better to build a more contained struct.
Also, the header could use `NamedTuple` instead of `Dict`.

And actually, by far there is no single use case where I need to read multiple files together. If you want to do so, just call the function twice.

Before v0.5.1, `readdata` function in Matlab for large data is 2 times faster than that in Julia. The reason is simply using `read` or `unsafe_read` in the lower level. The latter one is much faster. After the fix, Julia version performs 5 times faster than the Matlab version in reading binary data.

In the roadmap of PyCall 2.0, there will direct support for accessing Julia objects. I hesitate to do it myself, so let's just wait for it to come.

The support for a long string containing several filenames as inputs has been dropped. It should be substituted by an array of strings.

Right now the derived quantity plots are not supported. In order to achieve this, I may need:
- [x] A new function `get_var(data, filehead, string)` returning the derived variable
- [ ] A new plotting function that understands the derived data type

The first one is achieved by a trick I found on discourse, which basically identifies symbols as names to members in a struct.

Vector naming is messed up if you are using Tecplot VTK reader. For example, "B [nT]" --> "B [nT]_X", "B [nT]_Y", "B [nT]_Z". Not a big issue, but annoying.

There is a unit package in Julia [unitful](https://github.com/PainterQubits/Unitful.jl) for handling units. Take a look at that one if you really want to solve the unit problems.

I have encountered a very bad problem of corrupting binary *.vtu files. It turned out that the issue is the starting position of data is wrong because of the way I skip the header AUXDATA part. Sometimes the binary numbers may contain newline character that confuses the reader. It is now fixed. Later on the reading of the header part is completely rewritten to provide better support for a variety of Tecplot Ascii headers.

I have already made a lot of mistakes by mixing the row-major and column-major codes. Explicitly list all the parts that require extra care!

I have a new issue coming up with the interoperability with Python. I may need to split this package into pure IO and pure plotting to avoid the cross-dependency of Matplotlib. The idea is that PyPlot is only needed when I want to quickly scan through the data!

- [ ] Full coverage of tests
- [x] PyJulia support for manipulating data directly in Python
- [x] Derived variable support
- [x] General postprocessing script for concatenating and converting files.
- [x] Replace np.meshgrid with list comprehension
- [ ] Find a substitution of triangulation in Julia
- [ ] Allow dot syntax to get dictionary contents (Base.convert?)
- [ ] Binary library support
