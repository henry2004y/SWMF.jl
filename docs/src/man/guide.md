# Guide

The VTK files does not have timestep information. To allow for further time series processing in Paraview, a script `create_pvd.jl` is provided for generating the pvd container.

## Tricks

- This is the first time I use Julia for reading general ascii/binary files. It was a pain at first due to the lack of examples and documents using any basic function like read/read!, but fortunately I figured them out myself. One trick in reading binary array data is the usage of view, or subarrays, in Julia. In order to achieve that, I have to implement my own `read!` function in addition to the base ones.
- Tecplot and VTK unstructured data formats have the same connectivity ordering for hexahedron, but different ordering for voxel (in VTK). A function `swaprows` is implemented to switch the orderings.
- Because of the embarrassing parallelism nature of postprocessing, it is quite easy to take advantage of parallel approaches to process the data.

## Issues

At first I forgot to export the Data struct, so everytime when I modified the code and rerun plotdata, it will shout error at me, saying no type was found for the input type.

The current support of animation in Matplotlib is not good enough, especially for interactive plotting and scanning through multiple snapshots. The color range is constantly giving me headaches.

The current wrapper over Matplotlib makes it difficult to modify the plots afterwards, which especially causes problems when dealing with time series snapshots. The colorbar is so hard to fix.

In the roadmap of PyCall 2.0, there will direct support for accessing Julia objects. I hesitate to do it myself, so let's just wait for it to come.

The support for a long string containing several filenames as inputs has been dropped. It should be substituted by an array of strings.

Right now the derived quantity plots are not supported. In order to achieve this, I may need:
- [x] A new function `get_var(data, filehead, string)` returning the derived variable
- [ ] A new plotting function that understands the derived data type

The first one is achieved by a trick I found on discourse, which basically identifies symbols as names to members in a struct.

There is a user recipe in Plots. This is exactly what I am looking for, but more issues are coming up. I have created a new branch for this development.

I want to do scattered interpolation in Julia directly, but I have not found a simple solution to do this.

A direct wrapper over PyPlot function is possible, and would be more suitable for passing arguments. This may be a more plausible way to go than relying on recipes.

When doing processing in batch mode on a cluster, there's usually no need to render the plots on screen. There exists such a backend for this purpose:
```
using PyPlot
PyPlot.matplotlib.use("Agg")
```
However, notice that currently Agg backend does not support draw_artist. For example, you cannot add an anchored text to your figure.

Vector naming is messed up if you are using Tecplot VTK reader. For example, "B [nT]" --> "B [nT]_X", "B [nT]_Y", "B [nT]_Z". Not a big issue, but annoying.

There is a unit package in Julia [unitful](https://github.com/PainterQubits/Unitful.jl) for handling units. Take a look at that one if you really want to solve the unit problems.

I have encountered a very bad problem of corrupting binary *.vtu files. It turned out that the issue is the starting position of data is wrong because of the way I skip the header AUXDATA part. Sometimes the binary numbers may contain newline character that confuses the reader. It is now fixed.

- [x] Fixed colorbar control through Matplotlib
- [ ] Test suite for checking validity
- [x] Cuts from 3D data visualization besides contour
- [ ] Switch to Makie for 3D plotting and animation
- [ ] PyBase support for manipulating data directly in Python
- [x] Derived variable support
- [ ] General postprocessing script for concatenating and converting files.
- [x] Direct wrapper over matplotlib functions to get seamless API
- [x] Replace np.meshgrid with list comprehension
- [ ] Find a substitution of triangulation in Julia
- [ ] Allow dot syntax to get dictionary contents (Base.convert?)