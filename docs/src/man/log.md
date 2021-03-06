# Development Log

All the workflows here is not restricted to one type of model output. After being familiar with new ideas and new models, one can easily make use of existing samples and create reader of their own.
Because of the embarrassing parallelism nature of postprocessing, it is quite easy to take advantage of parallel approaches to process the data.

This is the first time I use Julia for reading general ascii/binary files. It was a pain at first due to the lack of examples and documents using any basic function like read/read!, but fortunately I figured them out myself. One trick in reading binary array data is the usage of view, or subarrays, in Julia. In order to achieve that, I have to implement my own `read!` function in addition to the base ones.
Before v0.5.1, `readdata` function in Matlab for large data is 2 times faster than that in Julia. The reason is simply using `read` or `unsafe_read` in the lower level. The latter one is much faster. After the fix, Julia version performs 5 times faster than the Matlab version in reading binary data.

I am enlightened by the way SpacePy handles files. Instead of a separate header and data array, it may be better to build a more contained struct.
Also, the header could use `NamedTuple` instead of `Dict`.

## Reading Multiple Files

And actually, by far there is no single use case where I need to read multiple files together. If you want to do so, just call the function twice.

## Interoperability

Demos are provided for calling Matlab/Python directly from Julia for debugging and testing. This part will later be separated out for potential Python and Matlab users. Currently the plotting and interpolation needed during plotting are done in Python. For instance, the 3D scatterred interpolation is done via `Interpolate` in Scipy. Hopefully these additional dependencies will be cut down.

In the current version of PyCall and PyJulia, there is already direct support for accessing Julia struct objects (noted as `jlwrap`).

I have a new issue coming up with the interoperability with Python. I may need to split this package into pure IO and pure plotting to avoid the cross-dependency of Matplotlib. The idea is that PyPlot is only needed when I want to quickly scan through the data!

## Units

There is a unit package in Julia [unitful](https://github.com/PainterQubits/Unitful.jl) for handling units.
Take a look at that one if you really want to solve the unit problem.
On top of that, [UnitfulRecipes.jl](https://github.com/jw3126/UnitfulRecipes.jl) provides integrity with [Plots.jl](http://docs.juliaplots.org/latest/).

The ideas of how to make abstractions is more important than the unit conversion itself.

## C Dependency

A real open-source project is a collaborated work not only from a bunch of people, but also a group of languages. In Julia, this can be achieved with the help of the [Package manager](https://julialang.github.io/Pkg.jl/dev/).

I want to have some C dependencies in my code instead of rewriting everything in Julia. This would serve as an attempt to quickly make things work.

Right now this seems to be a little bit difficult for me. I need to learn from experts. The tracing scheme in C is rewritten in Julia so I don't need to bother for now.
Checkout [BinaryBuilder](https://juliapackaging.github.io/BinaryBuilder.jl/latest/#Project-flow-1) for more information.
A nice example is given in [this C package](https://github.com/jakubwro/SineWaves.jl).

## VTK

The VTK files does not have timestep information. To allow for further time series processing in Paraview, a script `create_pvd.jl` is provided for generating the `pvd` container. After dicussing with the author of [WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl), this is finally achieved by adding a scalar that is not related to grid in VTK.

Currently, multi-block (VTM), rectilinear (VTR), and unstructured (VTU) conversion are supported.
By default the file size will be reduced with compression level 6, but the actual compression ratio depends on the original data.

For the plotting, streamline tracing and particle tracing, a common problem is the grid and related interpolation process. I am envisioning a more general approach to deal with block-based and unstructured grid to provide fundamental support for all of these.

Currently, BATSRUS output contains only cell center coordinates and cell center values (referring as `tcp` in `PARAM.in`), or node coordinates and interpolated nodal values (referring as `tec`). It is recommended to save in the `tcp` format to avoid interpolation. In principle VTK also supports the combination of node coordinates and cell center values, but it is not necessary here.

The simple multiblock VTK format conversion can only deal with data within a block, but it cannot figure out the connectivities between neighboring blocks. To fill the gap between blocks, we have to retrieve the tree data stored in `.tree` files. This requires a in-depth understanding of the grid tree structure in BATSRUS, and it took me a whole week to think, write and debug the code!

Information contained in `iTree_IA`:
* the status of the node (used, unused, to be refined, to be coarsened, etc.);
* the current, the maximum allowed and minimum allowed AMR levels for this node;
* the three integer coordinates with respect to the whole grid;
* the index of the parent node (if any);
* the indexes of the children nodes (if any);
* the processor index where the block is stored for active nodes;
* the local block index for active nodes.

Basically, it requires several steps:
1. Obtain the neighbor block indexes.
2. Obtain the relative AMR level between neighbor blocks.
3. Calculate the global cell indexes.
4. Fill in the connectivity list.

Several issues worth noticing:
* The blocks are load balanced in Morton curve ordering.
* Rules are needed to decide which block is in charge of writing the connectivities. Following the `ModWriteTecplot` implementation,
  * connectivities between blocks on the same AMR level are written by the left neighbors.
  * connectivities between blocks on different AMR levels are written by the finer blocks. This choice may simplify the logic in writing.
* 2D formatted outputs (with `dxPlot⩾0`) will generate duplicate points for 2D runs after postprocessing. I don't fully understand why.
* 2D unformatted outputs (with `dxPlot<0`) can generate correct point coordinates, but the output points in postprocessing are sorted based on `x+e*y+e^2*z` to allow removing duplicate points. This is mostly useful for 2D cuts of 3D simulations where duplicate points are possible to exist, and should be merged/averaged. This also means that it is literally impossible to figure out the original cell ordering, and impossible to convert to 2D VTU format (unless you construct some new connectivities based on triangulation)! However, it guarantees that the cut outputs from parallel runs are identical.
* 3D unformatted output as stored as it is.
* The simulation is typically done with multiple MPI processes. In the tree structure, each node only records the MPI process local block index. What makes it more complicated is that the local block indexes on a given MPI process may have skipped numbers because they are not in the physical simulation domain. This means that in order to get the true global block indexes, one needs to count how many blocks are used on all MPI processes with lower ranks, and also the order of the current node based on local block index on this MPI process.
* It would be pretty difficult to count the number of elements in the connectivity list. Appending to an array is a very expensive operation: counting the elements in the first loop, preallocating the array and then repeat the same computation to fill in the list results in over 1000 times speedup.
* The implementation in v0.2.3 inherits the bug from BATSRUS `ModWriteTecplot` that will generate duplicate voxel grid. For example, imagine a coarse block which is surrounded by refined blocks all around. On edge 8, the two refined face neighbor will both think they are in charge of writing the connectivity for the cells around the edge. This does not cause any serious issues as of now, so I don't fix it. Speaking of fixing that, there are two ways:
  * `fillCellNeighbors!` not only deals with coarsened neighbors, but also the refined neighbors. This means that literally we need to implement the `message_pass_cell` completely. We still do a two-round process, write everything without checking, and in the end before IO remove the duplicate rows. This may be even faster then implementing some complicated writing rules to avoid duplicating.
  *  Polish the writing rules even further. This requires me to consider all possible cases. As a first step, for the above example case, neither of the two face neighbor should write the connectivity, but instead the edge neighbor should write.
* Many function names are inherited from BATL. I may consider rename them in the future if more general task is required.

### Ordering of connectivity

Tecplot and VTK unstructured data formats have the same connectivity ordering for hexahedron, but different ordering for voxel (in VTK). A function `swaprows` is implemented to switch the orderings.

### Variable naming

Vector naming is messed up if you are using Tecplot VTK reader. For example, "B [nT]" --> "B [nT]_X", "B [nT]_Y", "B [nT]_Z". Not a big issue, but annoying.

### AUXDATA

I have encountered a very bad problem of corrupting binary *.vtu files. It turned out that the issue is the starting position of data is wrong because of the way I skip the header AUXDATA part. Sometimes the binary numbers may contain newline character that confuses the reader. It is now fixed.
Later on the reading of the header part is completely rewritten to provide better support for a variety of Tecplot Ascii headers.
All the AXUDATA information is now stored into global VTK data.

### Native VTK output

In the future versions of BATSRUS, we should be able to output VTK files directly with [VTKFortran](https://github.com/szaghi/VTKFortran). I won't do it now.

### Custom VTK reader

ParaView allows for custom Python reader. Examples can be found in Chapter 12.4 in the official manual, and an example of full Python plugins can be found at Kiware's [gitlab](https://gitlab.kitware.com/paraview/paraview/blob/master/Examples/Plugins/PythonAlgorithm/PythonAlgorithmExamples.py) page.

The XML package not only provide writing into XML files, but also reading XML structures. Therefore, if you want you can create a VTK reader.

### AMR Grid Structure

`vtkOverlappingAMR` implements a somewhat strict Berger-Collela AMR scheme:
1. All grids are Cartesian.
2. Grids at the same level do not overlap.
3. The refinement ratios, RL, between adjacent levels are integer (typically 2 or 4) and uniform within the same level.
4. Grid cells are never partially refined; i.e., each cell is refined to four quads in 2D or eight hexahedra in 3D.

Or in other words,
* Refinement ratio across levels is constant.
* Each block at levels > 0 need to be covered 100% by one parent block of
previous level.
* Some other restriction about what happens at the boundary.

You can directly use `vtkUniformGridAMR`, which does not impose any
restrictions. Most filters should work for this class - there just wouldn't
be any specialized filters such as the dual-grid contour / clip ones for
the `vtkOverlappingAMR`.

The `vtkAMRInformation` documentation consists only of
* Refinement ratio between AMR levels
* Grid spacing for each level
* The file block index for each block parent child information, if requested

![sample_2DAMR](../images/sample_2DAMR.jpg)
Sample 2D AMR Dataset with two levels and refinement ratio, RL=4. The root level (L0) consists of a single grid shown in black wireframe while the next level (L1) consists of two grids, depicted in green wireframe and red wireframe respectively. The two grids at L1 are projected from the root level to illustrate that the cells underneath are “hidden.”

In VTK, the collection of AMR grids is stored in a `vtkHierarchicalBoxDataSet` data-structure. Each grid, G(Li,k), is represented by a `vtkUniformGrid` data structure where the unique key pair (Li,k) denotes the corresponding level (Li) and the grid index within the level (k) with respect to the underlying hierarchical structure. An array historically known as `IBLANK`, stored as a cell attribute in `vtkUniformGrid`, denotes whether a cell is hidden or not. The blanking array is subsequently used by the mapper to hide lower resolution cells accordingly when visualizing the dataset.

To enable the execution of data queries without loading the entire dataset in memory, metadata information is employed. The metadata stores a minimal set of geometric information for each grid in the AMR hierarchy. Specifically, the AMR metadata, B(Li,k), corresponding to the grid G(Li,k), is represented using a `vtkAMRBox` object and it consists of the following information:

1. N={Nx, Ny, Nz} — the cell dimensions of the grid (since the data is cell-centered)
2. The grid spacing at level L, hL={hx,hy,hz}
3. The grid level Li and grid index k
4. The global dataset origin, X=(X0, Y0, Z0), i.e., the minimum origin from all grids in level L0
5. The LoCorner and HiCorner, which describe the low and high corners of the rectangular region covered by the corresponding grid in a virtual integer lattice with the same spacing (h) that covers the entire domain.

![sample_2DAMR](../images/sample_AMRmetadata.png)

Given the metadata information stored in the AMR box of each grid, the refinement ratio at each level can be easily computed using relationship (1) from Table 1. Further, the cartesian bounds the corresponding grid covers and the number of points and cells is also available (see relationships 2-4 in Table 1). Notably, geometric queries such as determining which cell contains a given point, or if a grid intersects a user-supplied slice plane, can be answered using just the metadata.

There is a `vtkAMRDualExtractionFilter`, which constructs a dual-mesh (i.e., the mesh constructed by connecting the cell-centers) over the computational domain.
If we can directly tell ParaView that the mesh we have is a dual-mesh, then the initial trial with multi-block data may work directly.

`AMRGaussianPulseSource`

See [Multi-Resolution Rendering with Overlapping AMR](https://www.paraview.org/ParaView/index.php/Multi-Resolution_Rendering_with_Overlapping_AMR) for the implementation of C++ reader in VTK.

## Support on Derived Variables

Right now the derived quantity plots are not supported. In order to achieve this, I may need:
- [x] A new function `get_var(data, filehead, string)` returning the derived variable
- [ ] A new plotting function that understands the derived data type

The first one is achieved by a trick I found on discourse, which basically identifies symbols as names to members in a struct.
This test feature is not ideal and will be dropped in later versions. 
This looks like the Python Calculator in ParaView.
I don't know how ParaView achieve this, but in Python it is pretty easy to parse a string and execute it as some new codes using `exec` function, as Yuxi did in his interface to the yt library.
Julia has equivalent metaprogramming capabilities, but there are many restrictions to apply it in practice so it is generally adviced to avoid evaluating expressions inside functions during runtime.
Another option is to create functions for derived variables.
This is actually pretty good both in terms of performance and conciseness.
The idea is to create a dictionary of derived variable names as keys and anonymous functions as values, and if the input string is found in the dictionary, call the corresponding function to obtain the value.
This has been successfully tested in my new scripts for processing Vlasiator outputs, and can be directly ported here for BATSRUS.

## Extra Notes

I have already made a lot of mistakes by mixing the row-major and column-major codes. Explicitly list all the parts that require extra care!

## Todo List

- [x] Full coverage of tests
- [x] PyJulia support for manipulating data directly in Python
- [x] Derived variable support
- [x] General postprocessing script for concatenating and converting files.
- [x] Replace np.meshgrid with list comprehension
- [x] Drop the support for a long string containing several filenames; substitute by an array of strings.
- [ ] Find a substitution of triangulation in Julia
- [ ] Allow dot syntax to get dictionary contents (Base.convert?)
- [ ] Binary library support
- [ ] VTK reader
