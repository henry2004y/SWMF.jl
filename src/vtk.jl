# Convert full BATSRUS AMR output to VTK.

using FortranFiles, WriteVTK

export convertTECtoVTU, convertIDLtoVTK, readhead, readtree, getConnectivity
export Batl

# Named indexes of iTree_IA
const status_   =  Int8(1) # used, unused, to be refined, to be coarsened, etc.
const level_    =  Int8(2) # grid AMR level for this node
const proc_     =  Int8(3) # processor index where the block is stored for active nodes
const block_    =  Int8(4) # local block index for this node
const minLevel_ =  Int8(5) # minimum AMR level allowed for this node
const maxLevel_ =  Int8(6) # maximum AMR level allowed for this node
const coord0_   =  Int8(6) # equal to coord1_-1
const coord1_   =  Int8(7) # coordinate of node in 1st dimension
const coord2_   =  Int8(8) # coordinate of node in 2nd dimension
const coord3_   =  Int8(9) # coordinate of node in 3rd dimension
const coordLast_=  Int8(9) # coord0_ + MaxDim (?)
const parent_   = Int8(10) # parent_ must be equal to child0_
const child0_   = Int8(10) #
const child1_   = Int8(child0_ + 1)
#const childLast_= child0_ + nChild # This can be replaced with end

# Possible values for the status variable
const unset_   = Int8(-100) # index for unset values (that are otherwise larger)
const unused_     = Int8(-1) # unused block (not a leaf)
const refine_     = Int8(-2) # parent block to be refined
const noCoarsen_  = Int8(-3) # block not to be coarsened
const coarsen_    = Int8(-4) # child block to be coarsened
const used_       = Int8( 1) # currently used block (leaf)
const refineNew_  = Int8( 2) # child block to be refined
const refined_    = Int8( 3) # refined child block
const coarsenNew_ = Int8( 4) # parent block to be coarsened
const coarsened_  = Int8( 5) # coarsened parent block

# Deepest AMR level relative to root nodes (limited by 32 bit integers)
const maxLevel = 30

# The maximum integer coordinate for a given level below root nodes
const maxCoord_I = Int32[2^i for i in 0:maxLevel]

"BATSRUS output standalone header information."
struct Head
   # Number of cells per block in each direction.
   # These values are set by the Config.pl script.
   # Set 1 for ignored directions!
   nI::Int32
   nJ::Int32
   nK::Int32
   # Maximum number of ghost cells set by Config.pl script.
   # Valid values are 0,1,2,3,4,5
   nG::Int32
   # Refinement ratios in the 3 dimensions. Either 1 or 2.
   # The values are set by the Config.pl script.
   iRatio::Int32
   jRatio::Int32
   kRatio::Int32
   # Number of dimensions in which grid adaptation is done
   nDimAmr::Int32
   # Number of children per node
   nChild::Int32
   # Number of root blocks
   nRoot_D::Vector{Int32}
   # Grid lower and upper edges
   CoordMin_D::Vector{Float64}
   CoordMax_D::Vector{Float64}
   # Indexes of AMR dimensions
   iDimAmr_D::Vector{Int32}
   isPeriodic_D::Vector{Bool}
   # Plotting resolution
   dxPlot_D::Vector{Float64}
end

struct Batl
   head::Head
   iTree_IA::Array{Int32,2}
   iRatio_D::Vector{Int32} # Array of refinement ratios
   nDim::Int8
end


"""
	convertTECtoVTU(head, data, connectivity, filename="out")

Convert unstructured Tecplot data to VTK. Note that if using voxel type data
in VTK, the connectivity sequence is different from Tecplot.
Note that the 3D connectivity sequence in Tecplot is the same with the
`hexahedron` type in VTK, but different with the `voxel` type.
The 2D connectivity sequence is the same as the `quad` type, but different with
the `pixel` type.
For example, in 3D the index conversion is:
```
# PLT to VTK voxel index_ = [1 2 4 3 5 6 8 7]
for i = 1:2
   connectivity = swaprows!(connectivity, 4*i-1, 4*i)
end
```
"""
function convertTECtoVTU(head, data, connectivity, filename="out")

   nVar = length(head.variables)
   points = @view data[1:head.nDim,:]
   cells = Vector{MeshCell{VTKCellType,Array{Int32,1}}}(undef,head.nCell)

   if head.nDim == 3
      @inbounds for i = 1:head.nCell
         cells[i] = MeshCell(VTKCellTypes.VTK_HEXAHEDRON, connectivity[:,i])
      end
   elseif head.nDim == 2
      @inbounds for i = 1:head.nCell
         cells[i] = MeshCell(VTKCellTypes.VTK_QUAD, connectivity[:,i])
      end
   end

   vtkfile = vtk_grid(filename, points, cells)

   for ivar = head.nDim+1:nVar
      if endswith(head.variables[ivar],"_x") # vector
         if head.nDim == 3
            var1 = @view data[ivar,:]
            var2 = @view data[ivar+1,:]
            var3 = @view data[ivar+2,:]
            namevar = replace(head.variables[ivar], "_x"=>"")
            vtk_point_data(vtkfile, (var1, var2, var3), namevar)
         elseif head.nDim == 2
            var1 = @view data[ivar,:]
            var2 = @view data[ivar+1,:]
            namevar = replace(head.variables[ivar], "_x"=>"")
            vtk_point_data(vtkfile, (var1, var2), namevar)
         end
      elseif endswith(head.variables[ivar],r"_y|_z")
         continue
      else
         var = @view data[ivar,:]
         vtk_point_data(vtkfile, var, head.variables[ivar])
      end
   end

   # Add meta data from Tecplot AUXDATA
   for i in 1:length(head.auxdata)
      vtkfile[head.auxdataname[i],VTKFieldData()] = head.auxdata[i]
   end

   outfiles = vtk_save(vtkfile)
end

"""
	convertIDLtoVTK(filename; dir=".", gridType=1, verbose=false)

Convert 3D BATSRUS *.out to VTK. If `gridType==1`, it converts to the 
rectilinear grid; if `gridType==2`, it converts to the structured grid.
If `filename` does not end with "out", it tries to find the ".info" and ".tree"
file with the same name tag and generates 3D unstructured VTU file.
"""
function convertIDLtoVTK(filename::AbstractString; dir=".", gridType=1,
   verbose=false)

   if endswith(filename, ".out")
      data = readdata(filename, dir=dir)

      nVar = length(data.head.wnames)

      outname = filename[1:end-4]
   
      if gridType == 1 # rectilinear grid
         x = @view data.x[:,1,1,1]
         y = @view data.x[1,:,1,2]
         z = @view data.x[1,1,:,3]
   
         outfiles = vtk_grid(outname, x,y,z) do vtk
            for ivar = 1:nVar
               if data.head.wnames[ivar][end] == 'x' # vector
                  var1 = @view data.w[:,:,:,ivar]
                  var2 = @view data.w[:,:,:,ivar+1]
                  var3 = @view data.w[:,:,:,ivar+2]
                  namevar = data.head.wnames[ivar][1:end-1]
                  vtk_point_data(vtk, (var1, var2, var3), namevar)
               elseif data.head.wnames[ivar][end] in ('y','z')
                  continue
               else
                  var = @view data.w[:,:,:,ivar]
                  vtk_point_data(vtk, var, data.head.wnames[ivar])
               end
            end
         end
      elseif gridType == 2 # structured grid
         xyz = permutedims(data.x, [4,1,2,3])
   
         outfiles = vtk_grid(outname, xyz) do vtk
            for ivar = 1:nVar
               if data.head.wnames[ivar][end] == 'x' # vector
                  var1 = @view data.w[:,:,:,ivar]
                  var2 = @view data.w[:,:,:,ivar+1]
                  var3 = @view data.w[:,:,:,ivar+2]
                  namevar = data.head.wnames[ivar][1:end-1]
                  vtk_point_data(vtk, (var1, var2, var3), namevar)
               elseif data.head.wnames[ivar][end] in ('y','z')
                  continue
               else
                  var = @view data.w[:,:,:,ivar]
                  vtk_point_data(vtk, var, data.head.wnames[ivar])
               end
            end
         end
      elseif gridType == 3 # unstructured grid
         @error "No tree information for conversion!"
      end

   else
      # info, tree, and out files
      data = readdata(filename*".out", dir=dir)
      batl = Batl(readhead(filename*".info"), readtree(filename)...)
      connectivity = getConnectivity(batl)

      outname = filename

      nDim = batl.nDim
      nVar = length(data.head.wnames)
      nCell = size(connectivity, 2)
      if nDim == 3
         if batl.head.dxPlot_D[1] ≥ 0.0
            @error "Why are there duplicate points?"
         else
            points = data.x[:,1,1,:]'
         end
      elseif nDim == 2
         if batl.head.dxPlot_D[1] ≥ 0.0 # points are not sorted in postproc.f90
            @error "Why are there duplicate points? Ask!"
            points = data.x[:,1,:]'
         else # points are sorted in postproc.f90!
            @error "point original order cannot be retrieved!"
         end
      end
      cells = Vector{MeshCell{VTKCellType,Array{Int32,1}}}(undef,nCell)
   
      if nDim == 3
         @inbounds for i = 1:nCell
            cells[i] = MeshCell(VTKCellTypes.VTK_HEXAHEDRON, connectivity[:,i])
         end
      elseif nDim == 2
         @inbounds for i = 1:nCell
            cells[i] = MeshCell(VTKCellTypes.VTK_QUAD, connectivity[:,i])
         end
      end
   
      vtkfile = vtk_grid(filename, points, cells)
   
      for ivar = 1:nVar
         if endswith(data.head.wnames[ivar], "x") # vector
            if nDim == 3
               var1 = @view data.w[:,1,1,ivar]
               var2 = @view data.w[:,1,1,ivar+1]
               var3 = @view data.w[:,1,1,ivar+2]
               var = (var1, var2, var3)
               vtkfile[data.head.wnames[ivar][1:end-1], VTKPointData()] = var
            elseif nDim == 2 # not sure how VTK handles 2D vector!
               var1 = @view data.w[:,1,ivar]
               var2 = @view data.w[:,1,ivar+1]
               vtkfile[data.head.wnames[ivar][1:end-1], VTKPointData()] = (var1, var2)
            end
         elseif endswith(data.head.wnames[ivar], r"y|z")
            continue
         else
            if nDim == 3
               var = @view data.w[:,1,1,ivar]
            elseif nDim == 2
               var = @view data.w[:,1,ivar]
            end
            vtkfile[data.head.wnames[ivar], VTKPointData()] = var
         end
      end

      outfiles = vtk_save(vtkfile)
   end

   verbose && @info "$(filename) finished conversion."
   return
end

"Return matrix X with swapped rows i and j."
function swaprows!(X, i, j)
   m, n = size(X)
   if (1 ≤ i ≤ n) && (1 ≤ j ≤ n)
      @inbounds @simd for k = 1:n
         X[i,k],X[j,k] = X[j,k],X[i,k]
      end
      return X
   else
      throw(BoundsError())
   end
end

"Return header from info file. Currently only designed for 2D and 3D."
function readhead(filehead)

   nDim = 3
   nI, nJ, nK = 1, 1, 1
   nRoot_D = Int32[1,1,1]
   # Make sure that the thickness is unity in the ignored dimensions.
   CoordMin_D = fill(-0.5, 3)
   CoordMax_D = fill(0.5, 3)
   isPeriodic_D = fill(false, 3)
   dxPlot_D = fill(0.0,3)

   open(filehead) do f
      while !eof(f)
         ln = readline(f)
         if startswith(ln, "#NDIM")
            nDim = parse(Int32, split(readline(f))[1])
         elseif startswith(ln, "#GRIDBLOCKSIZE")
            nI = parse(Int32, split(readline(f))[1])
            nJ = parse(Int32, split(readline(f))[1])
            if nDim == 2
               nK = Int32(1)
            elseif nDim == 3
               nK = parse(Int32, split(readline(f))[1])
            end
         elseif startswith(ln, "#ROOTBLOCK")
            nRoot_D[1] = parse(Int32, split(readline(f))[1])
            nRoot_D[2] = parse(Int32, split(readline(f))[1])
            if nDim == 3
               nRoot_D[3] = parse(Int32, split(readline(f))[1])
            end
         elseif startswith(ln, "#GRIDGEOMETRYLIMIT")
            typeGeometry = split(readline(f))[1]
            CoordMin_D[1] = parse(Float64, split(readline(f))[1])
            CoordMax_D[1] = parse(Float64, split(readline(f))[1])
            CoordMin_D[2] = parse(Float64, split(readline(f))[1])
            CoordMax_D[2] = parse(Float64, split(readline(f))[1])
            if nDim == 3
               CoordMin_D[3] = parse(Float64, split(readline(f))[1])
               CoordMax_D[3] = parse(Float64, split(readline(f))[1])
            end
         elseif startswith(ln, "#PERIODIC")
            str = split(readline(f))[1]
            isPeriodic_D[1] = str == "T" ? true : false
            str = split(readline(f))[1]
            isPeriodic_D[2] = str == "T" ? true : false
            if nDim == 3
               str = split(readline(f))[1]
               isPeriodic_D[3] = str == "T" ? true : false
            end
         elseif startswith(ln, "#PLOTRESOLUTION")
            dxPlot_D[1] = parse(Float64, split(readline(f))[1])
            dxPlot_D[2] = parse(Float64, split(readline(f))[1])
            if nDim == 3
               dxPlot_D[3] = parse(Float64, split(readline(f))[1])
            end
         end  
      end
   end

   nG = 2 # most outputs are 2nd order

   iRatio, jRatio, kRatio = min(2, nI), min(2, nJ), min(2, nK)

   nDimAmr = iRatio + jRatio + kRatio - 3
      
   nChild = 2^nDimAmr

   # Get size of domain (in generalized coordinates)
   DomainSize_D = CoordMax_D - CoordMin_D

   # Indexes of AMR dimensions.
   # The magic formulas should be correct from 1 to nDimAmr.
   iDimAmrTmp_D = Int32[ 1 + (2-iRatio)*(3-jRatio), 6-iRatio-jRatio, 3 ]
   iDimAmr_D = iDimAmrTmp_D[1:nDimAmr]

   Head(nI, nJ, nK, nG, iRatio, jRatio, kRatio, nDimAmr, nChild, nRoot_D,
      CoordMin_D, CoordMax_D, iDimAmr_D, isPeriodic_D, dxPlot_D)
end

"Return BATL tree structure."
function readtree(filetag)

   iRatio_D = Int32[1,1,1]
   nRoot_D = Int32[1,1,1]

   # Loading AMR tree
   f = FortranFile(filetag*".tree")

   nDim, nInfo, nNode = read(f, Int32, Int32, Int32)
   iRatio_D[1:nDim] = read(f, (Int32,nDim)) # Array of refinement ratios
   nRoot_D[1:nDim] = read(f, (Int32,nDim)) # The number of root nodes in all dimension
   iTree_IA = read(f, (Int32,(nInfo,nNode)))

   close(f)

   return iTree_IA, iRatio_D, nDim
end


"""
   find_grid_block(batl, xyz_D)

Return block index that contains a point. Input location should be given in 
Cartesian coordinates.
"""
function find_grid_block(batl::Batl, xyz_D)

   CoordMin_D = batl.head.CoordMin_D
   CoordMax_D = batl.head.CoordMax_D

   Coord_D = xyz_D

   # Calculate normalized coordinates for tree search
   CoordTree_D = (Coord_D - CoordMin_D) ./ (CoordMax_D - CoordMin_D)

   if any(CoordTree_D .< 0.0) || any(CoordTree_D .> 1.0)
      iBlock = unset_
      return iBlock
   end

   # Find node containing the point
   iNode = find_tree_node(batl, CoordTree_D)

   # Check if point was found
   if iNode > 0
      # Convert to block and processor indexes
      iBlock = batl.iTree_IA[block_,iNode]
   else
      iBlock = unset_
   end

   return iBlock
end


"""
   find_tree_node(batl, Coord_D)

Find the node that contains a point. The point coordinates should be given in 
generalized coordinates normalized by the domain size.
"""
function find_tree_node(batl::Batl, Coord_D)

   nRoot_D = batl.head.nRoot_D
   nDimAmr = batl.head.nDimAmr
   iDimAmr_D = batl.head.iDimAmr_D
   iTree_IA = batl.iTree_IA

   # Scale coordinates so that 1 ≤ Coord_D ≤ nRoot_D+1
   Coord_D = @. 1.0 + nRoot_D*max(0.0, min(1.0, Coord_D))

   # Get root node index
   iRoot_D = min.(floor.(Int32, Coord_D), nRoot_D)

   # Root node indexes are ordered
   iNode = iRoot_D[1] + nRoot_D[1]*((iRoot_D[2]-1) + nRoot_D[2]*(iRoot_D[3]-1))

   if iTree_IA[status_,iNode] == used_
      return iNode
   end

   nLevelMax, iNodeMorton_I = order_tree(batl)

   # Get normalized coordinates within root node and scale it up
   # to the largest resolution: 0 <= iCoord_D <= maxCoord_I(nLevelMax)-1
   iCoord_D = min.(maxCoord_I[nLevelMax+1] - 1, 
      floor.(Int32, (Coord_D[iDimAmr_D] - iRoot_D[iDimAmr_D])*maxCoord_I[nLevelMax+1]))

   # Go down the tree using bit information
   for iLevel = nLevelMax-1:-1:0
      # Get the binary bits based on the coordinates
      iBit_D = ibits.(iCoord_D, iLevel, 1)

      # Construct child index as iChild = Sum Bit_i*2^i
      iChild = sum(iBit_D.*maxCoord_I[1:nDimAmr]) + child1_
      iNode  = iTree_IA[iChild,iNode]

      if iTree_IA[status_,iNode] == used_
         return iNode
      end
   end

   # Did not find the point so set iNode as unset
   iNode = unset_

   return iNode
end

"Logical shifts as the Fortran instrinsic function."
function ibits(i, pos, len)
   # Treat it as 32 bits integer as in BATL.
   ds = digits(i, base=2, pad=32)[pos+1:pos+len]
   s = zero(eltype(i))
   for val in ds
      s = s * 2 + val
   end
   return s
end


"""
   order_tree(batl)

Return maximum AMR level in the used block and the Morton curve order.
Set iNodeMorton_I indirect index arrays according to
1. root node order
2. Morton ordering for each root node
"""
function order_tree(batl::Batl)

   nRoot_D = batl.head.nRoot_D
   iTree_IA = batl.iTree_IA

   nNode = prod(nRoot_D)
   iNode = 0
   iMorton = 0

   # Get used node counts
   nNodeUsed = count(==(used_), iTree_IA[status_,:])

   iNodeMorton_I = fill(Int32(unset_), nNodeUsed)

   for kRoot = 1:nRoot_D[3]
      for jRoot = 1:nRoot_D[2]
         for iRoot = 1:nRoot_D[1]
            # Root nodes are the first ones
            iNode += 1
            # All root nodes are handled as if they were first child
            iMorton = order_children!(batl, iNode, iMorton, iNodeMorton_I)
         end
      end
   end

   nNodeUsed = iMorton

   # Set min and max refinement levels
   nLevelMin = maxLevel
   nLevelMax = 0
   for iMorton = 1:nNodeUsed
      iNode = iNodeMorton_I[iMorton]
      iLevel = iTree_IA[level_,iNode]
      nLevelMin = min(iLevel, nLevelMin)
      nLevelMax = max(iLevel, nLevelMax)
   end

   return nLevelMax, iNodeMorton_I
end


"""
   order_children!(iNode)

Recursively apply Morton ordering for nodes below a root block.
Store result into iNodeMorton_I and iMortonNode_A using the iMorton index.
"""
function order_children!(batl::Batl, iNode, iMorton, iNodeMorton_I)

   iTree_IA = batl.iTree_IA

   if iTree_IA[status_, iNode] ≥ used_
      iMorton += 1
      iNodeMorton_I[iMorton] = iNode
   else
      for iChild = child1_:size(iTree_IA,1)
         iMorton = order_children!(batl, iTree_IA[iChild, iNode], iMorton, iNodeMorton_I)
      end
   end

   return iMorton
end

"Find neighbors for any node in the tree. Only for Cartesian now."
function find_neighbor_for_anynode(batl::Batl, iNode)

   nDim = batl.nDim
   nRoot_D = batl.head.nRoot_D
   iTree_IA = batl.iTree_IA
   iRatio_D = batl.iRatio_D

   IsLatitudeAxis = false
   IsSphericalAxis = false
   IsCylindricalAxis = false

   # Periodicity check is necessary for recovering the boundary conditions,
   # but it makes postprocessing harder. An easy solution is to turn it off. 
   #isPeriodic_D = batl.head.isPeriodic_D
   isPeriodic_D = [false, false, false]

   DiLevelNei_III = fill(Int8(0), (3,3,3))
   iNodeNei_III = fill(Int32(-1), (4,4,4))

   # Get AMR level of the node
   iLevel = iTree_IA[level_,iNode]

   # Calculate scaling factor from integer index to 0<x,y,z<1 real coordinates
   Scale_D = 1.0 ./ nRoot_D
   for i = 1:3
      if iRatio_D[i] == 2
         Scale_D[i] /= maxCoord_I[iLevel+1]
      end
   end

   # Fill in self-referring info
   iNodeNei_III[2:3,2:3,2:3] .= Int32(iNode)
   DiLevelNei_III[2,2,2]     = 0

   # Loop through neighbors
   for k = 0:3
      Dk = round(Int8, (k - 1.5)/1.5)
      if nDim < 3
         if k != 1 continue end
         z = 0.3
      else
         z = (iTree_IA[coord3_, iNode] + 0.4*k - 1.1)*Scale_D[3]
         if z > 1.0 || z < 0.0
            if isPeriodic_D[3]
               z = mod(z, 1.0)
            elseif !IsLatitudeAxis
               iNodeNei_III[:,:,k+1] .= unset_
               DiLevelNei_III[:,:,Dk+2] .= unset_
               continue
            end
         end
      end
      # store z for spherical axis
      z0 = z
      for j = 0:3
         z = z0
         Dj = round(Int8, (j - 1.5)/1.5)
         if nDim < 2
            if j!=1 continue end
            y = 0.3
         else
            y = (iTree_IA[coord2_, iNode] + 0.4*j - 1.1)*Scale_D[2]
            if y > 1.0 || y < 0.0
               if isPeriodic_D[2]
                  y = mod(y, 1.0)
               elseif IsSphericalAxis
                  # Push back theta and go around half way in phi
                  y = max(0.0, min(1.0, y))
                  z = mod(z0+0.5, 1.0)
               else
                  iNodeNei_III[:,j+1,k+1] .= unset_
                  DiLevelNei_III[:,Dj+2,Dk+2] .= unset_
                  continue
               end
            end
            if z0 > 1.0 || z0 < 0.0
               # Push back latitude and go around half way in longitude
               z = max(0.0, min(1.0, z0))
               y = mod(y+0.5, 1.0)
            end
         end
         # store y for cylindrical axis case
         y0 = y
         for i = 0:3
            # Exclude inner points
            if 0<i<3 && 0<j<3 && 0<k<3 continue end

            Di = round(Int8, (i - 1.5)/1.5)

            # If neighbor is not finer, fill in the i=2 or j=2 or k=2 elements
            if DiLevelNei_III[Di+2,Dj+2,Dk+2] ≥ 0
               if i == 2
                  iNodeNei_III[i+1,j+1,k+1] = iNodeNei_III[2,j+1,k+1]
                  continue
               end
               if j == 2
                  iNodeNei_III[i+1,j+1,k+1] = iNodeNei_III[i+1,2,k+1]
                  continue
               end
               if k == 2
                  iNodeNei_III[i+1,j+1,k+1] = iNodeNei_III[i+1,j+1,2]
                  continue
               end
            end

            x = (iTree_IA[coord1_, iNode] + 0.4*i - 1.1)*Scale_D[1]
            y = y0
            if x > 1.0 || x < 0.0
               if isPeriodic_D[1]
                  x = mod(x, 1.0)
               elseif IsCylindricalAxis && x < 0.0
                  # Push back radius and go around half way in phi direction
                  x = 0.0
                  y = mod(y0+0.5, 1.0)
               else
                  iNodeNei_III[i+1,j+1,k+1] = unset_
                  DiLevelNei_III[Di+2,Dj+2,Dk+2] = unset_
                  continue
               end
            end

            jNode = find_tree_node(batl, [x, y, z])

            iNodeNei_III[i+1,j+1,k+1] = jNode
            DiLevelNei_III[Di+2,Dj+2,Dk+2] = iLevel - iTree_IA[level_, jNode]
         end
      end
   end

   return iNodeNei_III, DiLevelNei_III
end

"""
Return the mapping from local block index to global node index.
The block index stored in the tree is the processor local index, which should be
converted into global block index.
Note that there are gaps between local used block indexes! 
"""
function block_to_node(batl::Batl)

   nLevelMax, iNodeMorton_I = order_tree(batl)

   nNodeUsed = length(iNodeMorton_I)

   iNode_B = fill(Int32(0), nNodeUsed)
   iTree_IA = batl.iTree_IA

   for iMorton = 1:nNodeUsed
      iNode = iNodeMorton_I[iMorton]
      iBlock = iTree_IA[block_,iNode]
      iNode_B[iBlock] = iNode
   end

   return iNode_B
end

"Return global block index for the node."
function nodeToGlobalBlock(batl::Batl, iNode::Int32, nBlock_P)

   iTree_IA = batl.iTree_IA

   iProc = iTree_IA[proc_, iNode]

   localNodes_B = findall(x->x==iProc, iTree_IA[proc_,:])
   localBlocks_B = iTree_IA[block_,localNodes_B]
   myRank = findfirst(x->x==iNode, localNodes_B)

   localBlocksSorted_B = sort(localBlocks_B)
   localOrder = findfirst(x->x==localBlocks_B[myRank], localBlocksSorted_B)

   globalBlock = localOrder + nBlock_P[iProc+1]

   return globalBlock
end


function nodeToGlobalBlock(batl::Batl, iNodes::Array, nBlock_P)

   iTree_IA = batl.iTree_IA

   globalBlocks = similar(iNodes, Int32)

   for iNode = 1:length(iNodes)   
      iProc = iTree_IA[proc_, iNodes[iNode]]
      localNodes_B = findall(x->x==iProc, iTree_IA[proc_,:])
      localBlocks_B = iTree_IA[block_,localNodes_B]
      myRank = findfirst(x->x==iNodes[iNode], localNodes_B)
      localBlocksSorted_B = sort(localBlocks_B)
      localOrder = findfirst(x->x==localBlocks_B[myRank], localBlocksSorted_B)
      globalBlocks[iNode] = localOrder + nBlock_P[iProc+1]
   end

   return globalBlocks
end


"Get cell connectivity list."
function getConnectivity(batl::Batl)

   iTree_IA = batl.iTree_IA
   nI, nJ, nK = batl.head.nI, batl.head.nJ, batl.head.nK
   nDim = batl.nDim
   if nDim == 3
      nConn = 8
   elseif nDim == 2
      nConn = 4
      @error "2D not working currently!"
   elseif nDim == 1
      @error "Detected 1D data, connectivity not implemented!"
   end

   iCell_G = Array{Int32,3}(undef, nI+2, nJ+2, nK+2)

   # It would be difficult to calculate this beforehand.
   # The current implementation is to compute this in two rounds, where the
   # first round only do the counting.
   nElem = 0
   nBlockBefore = 0

   # Count accumulated number of blocks for MPI ranks in ascending order
   # Local block indexes may have gaps in between!
   nProc = maximum(iTree_IA[proc_,:]) + 1
   nBlock_P = fill(Int32(0), nProc)
   if nProc > 1
      for iProc = 1:nProc-1
         nBlock_P[iProc+1] = nBlock_P[iProc] + count(==(iProc-1), iTree_IA[proc_,:]) 
      end
   end

   # Pre-allocate just to let Julia know this variable.
   connectivity = Array{Int32,2}(undef, nConn, nElem)

   for iRound = 1:2
      if iRound == 2
         connectivity = Array{Int32,2}(undef, nConn, nElem)
         iElem = 0
         nBlockBefore = 0 # Reset
      end
      for iProc = 0:nProc-1
         localNodes_B = findall(x->x==iProc, iTree_IA[proc_,:])
         # It needs to be sorted for filling the correct iCell_G indexes!
         localBlocks_B = iTree_IA[block_,localNodes_B]
         seq = sortperm(localBlocks_B)
         localNodes_B = localNodes_B[seq] # sorted according to local block index

         nBlock = length(localNodes_B) # number of blocks on this processor

         for iBlock = 1:nBlock
            # initial cell index
            iCell = nI*nJ*nK*(nBlockBefore + iBlock - 1)

            @inbounds for k = 1:nK, j = 1:nJ, i = 1:nI
               iCell += 1
               iCell_G[i+1,j+1,k+1] = iCell
            end

            iNodeNei_III, DiLevelNei_III = find_neighbor_for_anynode(batl, localNodes_B[iBlock]) # rename!

            fillCellNeighbors!(batl, iCell_G, DiLevelNei_III, iNodeNei_III, nBlock_P)

            # In ModWriteTecplot, the first three conditions are needed only for certain 2D cases.
            if DiLevelNei_III[3,3,2] < 0
               iCell_G[end,end,:] .= 0
            end
            if DiLevelNei_III[3,2,3] < 0
               iCell_G[end,:,end] .= 0
            end
            if DiLevelNei_III[2,3,3] < 0
               iCell_G[:,end,end] .= 0
            end
            if nDim == 3 && DiLevelNei_III[3,3,3] < 0
               iCell_G[end,end,end] = 0
            end

            iMin = DiLevelNei_III[1,2,2] == 1 ? 0 : 1
            jMin = DiLevelNei_III[2,1,2] == 1 ? 0 : 1
            kMin = DiLevelNei_III[2,2,1] == 1 ? 0 : 1

            iMax = DiLevelNei_III[3,2,2] < 0 ? nI-1 : nI
            jMax = DiLevelNei_III[2,3,2] < 0 ? nJ-1 : nJ
            kMax = DiLevelNei_III[2,2,3] < 0 ? nK-1 : nK

            if nDim == 3
               for k = kMin:kMax, j = jMin:jMax, i = iMin:iMax
                  # Skip bricks that are not fully inside/usable
                  if any(iCell_G[i+1:i+2,j+1:j+2,k+1:k+2] .== 0)
                     continue
                  end
                  if iRound == 1
                     nElem += 1
                  else
                     iElem += 1
                  end
                  if iRound == 2
                     connectivity[:,iElem] = [
                        iCell_G[i+1,j+1,k+1],
                        iCell_G[i+2,j+1,k+1],
                        iCell_G[i+2,j+2,k+1],
                        iCell_G[i+1,j+2,k+1],
                        iCell_G[i+1,j+1,k+2],
                        iCell_G[i+2,j+1,k+2],
                        iCell_G[i+2,j+2,k+2],
                        iCell_G[i+1,j+2,k+2]]
                  end
               end
            elseif nDim == 2
               for j = jMin:jMax, i = iMin:iMax
                  connectivity = hcat(connectivity,
                     [
                     iCell_G[i+1,j+1,2],
                     iCell_G[i+2,j+1,2],
                     iCell_G[i+2,j+2,2],
                     iCell_G[i+1,j+2,2]]
                  )
               end
            end
         end
         nBlockBefore += nBlock # for next processor in ascending order
      end
   end

   return connectivity
end

"Return sibling index (1-8) for the given block node matrix."
function getSibling(iNodeNei_III, iTree_IA)
   iMyNode = iNodeNei_III[2,2,2]
   iMotherNode = iTree_IA[parent_,iMyNode]
   motherNode = iTree_IA[:,iMotherNode]
   iSibling = findfirst(x->x==iMyNode, motherNode[child1_:end])
end


"Fill neighbor cell indexes for the given block. Only tested for 3D. The faces,
edges, and vertices are ordered from left to right in x-y-z sequentially."
function fillCellNeighbors!(batl, iCell_G, DiLevelNei_III, iNodeNei_III, nBlock_P)

   iTree_IA = batl.iTree_IA
   nI, nJ, nK = batl.head.nI, batl.head.nJ, batl.head.nK
   nDim = batl.nDim

   nIJ = nI*nJ
   nIJK = nI*nJ*nK

   ## Faces

   # -x face
   if DiLevelNei_III[1,2,2] == 1
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,2,2], nBlock_P)

      @inbounds for k = 1:nK, j = 1:nJ
         iSibling = getSibling(iNodeNei_III, iTree_IA)

         if iSibling == 1
            iCell_G[1,j+1,k+1] = nIJK*(neiBlock-1) +
               nI*(1 + floor(Int, (j-1)/2)) + nIJ*floor(Int, (k-1)/2)
         elseif iSibling == 3
            iCell_G[1,j+1,k+1] = nIJK*(neiBlock-1) +
               nIJ/2 + nI*(1 + floor(Int, (j-1)/2)) + nIJ*floor(Int, (k-1)/2)
         elseif iSibling == 5
            iCell_G[1,j+1,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nI*(1 + floor(Int, (j-1)/2)) + nIJ*floor(Int, (k-1)/2)
         elseif iSibling == 7
            iCell_G[1,j+1,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nIJ/2 + nI*(1 + floor(Int, (j-1)/2)) + nIJ*floor(Int, (k-1)/2)
         end
      end
   end

   # +x face
   if DiLevelNei_III[3,2,2] == 0
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,3,3], nBlock_P)

      @inbounds for k = 1:nK, j = 1:nJ
         iCell_G[end,j+1,k+1] = nIJK*(neiBlock-1) + 1 + nI*(j-1) + nIJ*(k-1)
      end
   elseif DiLevelNei_III[3,2,2] == 1
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,3,3], nBlock_P)

      @inbounds for k = 1:nK, j = 1:nJ
         iSibling = getSibling(iNodeNei_III, iTree_IA)

         if iSibling == 2
            iCell_G[end,j+1,k+1] = nIJK*(neiBlock-1) +
               1 + nI*floor(Int, (j-1)/2) + nIJ*floor(Int, (k-1)/2)
         elseif iSibling == 4
            iCell_G[end,j+1,k+1] = nIJK*(neiBlock-1) +
               nIJ/2 + 1 + nI*floor(Int, (j-1)/2) + nIJ*floor(Int, (k-1)/2)
         elseif iSibling == 6
            iCell_G[end,j+1,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + 1 + nI*floor(Int, (j-1)/2) + nIJ*floor(Int, (k-1)/2)
         elseif iSibling == 8
            iCell_G[end,j+1,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nIJ/2 + 1 + nI*floor(Int, (j-1)/2) + nIJ*floor(Int, (k-1)/2)
         end
      end
   end

   # -y face
   if DiLevelNei_III[2,1,2] == 1
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[2,1,2], nBlock_P)

      @inbounds for k = 1:nK, i = 1:nI
         iSibling = getSibling(iNodeNei_III, iTree_IA)

         if iSibling == 1
            iCell_G[i+1,1,k+1] = nIJK*(neiBlock-1) +
               nI*(nJ-1) + 1 + floor(Int, (i-1)/2) + nIJ*floor(Int, (k-1)/2)
         elseif iSibling == 2
            iCell_G[i+1,1,k+1] = nIJK*(neiBlock-1) +
               nI*(nJ-1) + nI/2 + 1 + floor(Int, (i-1)/2) + nIJ*floor(Int, (k-1)/2)
         elseif iSibling == 5
            iCell_G[i+1,1,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nI*(nJ-1) + 1 + floor(Int, (i-1)/2) + nIJ*floor(Int, (k-1)/2)
         elseif iSibling == 6
            iCell_G[i+1,1,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nI*(nJ-1) + nI/2 + 1 + floor(Int, (i-1)/2) + nIJ*floor(Int, (k-1)/2)
         end
      end
   end

   # +y face
   if DiLevelNei_III[2,3,2] == 0
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[3,4,3], nBlock_P)

      @inbounds for k = 1:nK, i = 1:nI
         iCell_G[i+1,end,k+1] = nIJK*(neiBlock-1) + i + nIJ*(k-1)
      end
   elseif DiLevelNei_III[2,3,2] == 1
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[3,4,3], nBlock_P)

      @inbounds for k = 1:nK, i = 1:nI
         iSibling = getSibling(iNodeNei_III, iTree_IA)

         if iSibling == 3
            iCell_G[i+1,end,k+1] = nIJK*(neiBlock-1) +
               1 + floor(Int, (i-1)/2) + nIJ*floor(Int, (k-1)/2)
         elseif iSibling == 4
            iCell_G[i+1,end,k+1] = nIJK*(neiBlock-1) +
               nI/2 + 1 + floor(Int, (i-1)/2) + nIJ*floor(Int, (k-1)/2)
         elseif iSibling == 7
            iCell_G[i+1,end,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + 1 + floor(Int, (i-1)/2) + nIJ*floor(Int, (k-1)/2)
         elseif iSibling == 8
            iCell_G[i+1,end,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nI/2 + 1 + floor(Int, (i-1)/2) + nIJ*floor(Int, (k-1)/2)
         end
      end
   end

   # -z face
   if DiLevelNei_III[2,2,1] == 1
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[2,2,1], nBlock_P)

      @inbounds for j = 1:nJ, i = 1:nI
         iSibling = getSibling(iNodeNei_III, iTree_IA)

         if iSibling == 1
            iCell_G[i+1,j+1,1] = nIJK*neiBlock -
               nIJ + 1 + floor(Int, (i-1)/2) + nI*floor(Int, (j-1)/2)
         elseif iSibling == 2
            iCell_G[i+1,j+1,1] = nIJK*neiBlock -
               nIJ + nI/2 + 1 + floor(Int, (i-1)/2) + nI*floor(Int, (j-1)/2)
         elseif iSibling == 3
            iCell_G[i+1,j+1,1] = nIJK*neiBlock -
               nIJ/2 + 1 + floor(Int, (i-1)/2) + nI*floor(Int, (j-1)/2)
         elseif iSibling == 4
            iCell_G[i+1,j+1,1] = nIJK*neiBlock -
               nIJ/2 + nI/2 + 1 + floor(Int, (i-1)/2) + nI*floor(Int, (j-1)/2)
         end
      end
   end

   # +z face
   if DiLevelNei_III[2,2,3] == 0
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[3,3,4], nBlock_P)

      @inbounds for j = 1:nJ, i = 1:nI
         iCell_G[i+1,j+1,end] = nIJK*(neiBlock-1) + i + nI*(j-1)
      end
   elseif DiLevelNei_III[2,2,3] == 1
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[3,3,4], nBlock_P)

      @inbounds for j = 1:nJ, i = 1:nI
         iSibling = getSibling(iNodeNei_III, iTree_IA)

         if iSibling == 5
            iCell_G[i+1,j+1,end] = nIJK*(neiBlock-1) +
               1 + floor(Int, (i-1)/2) + nI*floor(Int, (j-1)/2)
         elseif iSibling == 6
            iCell_G[i+1,j+1,end] = nIJK*(neiBlock-1) +
               nI/2 + 1 + floor(Int, (i-1)/2) + nI*floor(Int, (j-1)/2)
         elseif iSibling == 7
            iCell_G[i+1,j+1,end] = nIJK*(neiBlock-1) +
               nIJ/2 + 1 + floor(Int, (i-1)/2) + nI*floor(Int, (j-1)/2)
         elseif iSibling == 8
            iCell_G[i+1,j+1,end] = nIJK*(neiBlock-1) +
               nIJ/2 + nI/2 + 1 + floor(Int, (i-1)/2) + nI*floor(Int, (j-1)/2)
         end
      end   
   end

   ## Edges, in total 12

   # edge 1
   if DiLevelNei_III[2,1,1] == 0
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[2,1,1], nBlock_P)

      @inbounds for i = 1:nI
         iCell_G[i+1,1,1] = nIJK*neiBlock - nI + i
      end
   elseif DiLevelNei_III[2,1,1] in (1,2)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[2,1,1], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)
   
      iAMR = 2*DiLevelNei_III[2,1,1]

      if iSibling == 1
         @inbounds for i = 1:nI
            iCell_G[i+1,1,1] = nIJK*neiBlock -
               nI + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 2
         @inbounds for i = 1:nI
            iCell_G[i+1,1,1] = nIJK*neiBlock -
               nI/2 + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 3
         @inbounds for i = 1:nI
            iCell_G[i+1,1,1] = nIJK*neiBlock -
               nIJ/2 - nI + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 4
         @inbounds for i = 1:nI
            iCell_G[i+1,1,1] = nIJK*neiBlock -
               nIJ/2 - nI/2 + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 5
         @inbounds for i = 1:nI
            iCell_G[i+1,1,1] = nIJK*(neiBlock-1) +
               nIJ*(nK/2-1) + nI*(nJ-1) + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 6
         @inbounds for i = 1:nI
            iCell_G[i+1,1,1] = nIJK*(neiBlock-1) +
               nIJ*(nK/2-1) + nI*(nJ-1) + nI/2 + 1 + fld(i-1, iAMR)
         end
      end
   end

   # edge 2
   if DiLevelNei_III[2,3,1] == 0
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[2,4,1], nBlock_P)

      @inbounds for i = 1:nI
         iCell_G[i+1,end,1] = nIJK*neiBlock - nIJ + i
      end
   elseif DiLevelNei_III[2,3,1] in (1,2)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[2,4,1], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)
   
      iAMR = 2*DiLevelNei_III[2,3,1]

      if iSibling == 1
         @inbounds for i = 1:nI
            iCell_G[i+1,end,1] = nIJK*neiBlock -
               nIJ/2 + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 2
         @inbounds for i = 1:nI
            iCell_G[i+1,end,1] = nIJK*neiBlock -
               nIJ/2 + nI/2 + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 3
         @inbounds for i = 1:nI
            iCell_G[i+1,end,1] = nIJK*neiBlock -
               nIJ + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 4
         @inbounds for i = 1:nI
            iCell_G[i+1,end,1] = nIJK*neiBlock -
               nIJ + nI/2 + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 7
         @inbounds for i = 1:nI
            iCell_G[i+1,end,1] = nIJK*(neiBlock-1) +
               nIJ*(nK/2-1) + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 8
         @inbounds for i = 1:nI
            iCell_G[i+1,end,1] = nIJK*(neiBlock-1) +
               nIJ*(nK/2-1) + nI/2 + 1 + fld(i-1, iAMR)
         end
      end
   end

   # edge 3
   if DiLevelNei_III[2,1,3] == 0
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[2,1,4], nBlock_P)
   
      @inbounds for i = 1:nI
         iCell_G[i+1,1,end] = nIJK*(neiBlock-1) + nI*(nJ-1) + i
      end      
   elseif DiLevelNei_III[2,1,3] in (1,2)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[2,1,4], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)
   
      iAMR = 2*DiLevelNei_III[2,1,3]

      if iSibling == 1
         @inbounds for i = 1:nI
            iCell_G[i+1,1,end] = nIJK*(neiBlock-1) +
               nIJK/2 + nI*(nJ-1) + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 2
         @inbounds for i = 1:nI
            iCell_G[i+1,1,end] = nIJK*(neiBlock-1) +
               nIJK/2 + nI*(nJ-1) + nI/2 + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 5
         @inbounds for i = 1:nI
            iCell_G[i+1,1,end] = nIJK*(neiBlock-1) +
               nI*(nJ-1) + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 6
         @inbounds for i = 1:nI
            iCell_G[i+1,1,end] = nIJK*(neiBlock-1) +
               nI*(nJ-1) + nI/2 + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 7
         @inbounds for i = 1:nI
            iCell_G[i+1,1,end] = nIJK*(neiBlock-1) +
               nI*(nJ/2-1) + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 8
         @inbounds for i = 1:nI
            iCell_G[i+1,1,end] = nIJK*(neiBlock-1) +
               nI*(nJ/2-1) + nI/2 + 1 + fld(i-1, iAMR)
         end
      end
   end

   # edge 4
   if DiLevelNei_III[2,3,3] == 0
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[2,4,4], nBlock_P)
      @inbounds for i = 1:nI
         iCell_G[i+1,end,end] = nIJK*(neiBlock-1) + i
      end
   elseif DiLevelNei_III[2,3,3] in (1,2)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[2,4,4], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      iAMR = 2^DiLevelNei_III[2,3,3]
   
      if iSibling == 3
         @inbounds for i = 1:nI
            iCell_G[i+1,end,end] = nIJK*(neiBlock-1) +
               nIJK/2 + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 4
         @inbounds for i = 1:nI
            iCell_G[i+1,end,end] = nIJK*(neiBlock-1) +
               nIJK/2 + nI/2 + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 5
         @inbounds for i = 1:nI
            iCell_G[i+1,end,end] = nIJK*(neiBlock-1) +
               nIJ/2 + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 6
         @inbounds for i = 1:nI
            iCell_G[i+1,end,end] = nIJK*(neiBlock-1) +
               nIJ/2 + nI/2 + 1 + fld(i-1, iAMR)
         end
      elseif iSibling == 7
         @inbounds for i = 1:nI
            iCell_G[i+1,end,end] = nIJK*(neiBlock-1) +
               1 + fld(i-1, iAMR)
         end
      elseif iSibling == 8
         @inbounds for i = 1:nI
            iCell_G[i+1,end,end] = nIJK*(neiBlock-1) +
               nI/2 + 1 + fld(i-1, iAMR)
         end
      end
   end

   # edge 5
   if DiLevelNei_III[1,2,1] == 0
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,2,1], nBlock_P)

      @inbounds for j = 1:nJ
         iCell_G[1,j+1,1] = nIJK*neiBlock - nIJ + nI*j
      end
   elseif DiLevelNei_III[1,2,1] in (1,2)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,2,1], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      iAMR = 2^DiLevelNei_III[1,2,1]
   
      if iSibling == 1
         @inbounds for j = 1:nJ
            iCell_G[1,j+1,1] = nIJK*neiBlock -
               nIJ + nI*(1 + fld(j-1,iAMR))
         end
      elseif iSibling == 2
         @inbounds for j = 1:nJ
            iCell_G[1,j+1,1] = nIJK*neiBlock -
               nIJ + nI/2 + nI*fld(j-1,iAMR)
         end
      elseif iSibling == 3
         @inbounds for j = 1:nJ
            iCell_G[1,j+1,1] = nIJK*neiBlock -
               nIJ/2 + nI*(1 + fld(j-1,iAMR))
         end
      elseif iSibling == 4
         @inbounds for j = 1:nJ
            iCell_G[1,j+1,1] = nIJK*neiBlock -
               nIJ/2 + nI/2 + nI*fld(j-1,iAMR)
         end
      elseif iSibling == 5
         @inbounds for j = 1:nJ
            iCell_G[1,j+1,1] = nIJK*neiBlock -
               nIJ*(nK/2+1) + nI*(1 + fld(j-1,iAMR))
         end
      elseif iSibling == 7
         @inbounds for j = 1:nJ
            iCell_G[1,j+1,1] = nIJK*neiBlock -
               nIJ*(nK/2+1) + nIJ/2 + nI*(1 + fld(j-1,iAMR))
         end
      end
   end

   # edge 6
   if DiLevelNei_III[3,2,1] == 0
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,2,1], nBlock_P)

      @inbounds for j = 1:nJ
         iCell_G[end,j+1,1] = nIJK*neiBlock - nIJ + 1 + nI*(j-1)
      end
   elseif DiLevelNei_III[3,2,1] in (1,2)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,2,1], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      iAMR = 2^DiLevelNei_III[3,2,1]

      if iSibling == 1
         @inbounds for j = 1:nJ
            iCell_G[end,j+1,1] = nIJK*neiBlock -
               nIJ + nI/2 + 1 + nI*fld(j-1,iAMR)
         end
      elseif iSibling == 2
         @inbounds for j = 1:nJ
            iCell_G[end,j+1,1] = nIJK*neiBlock -
               nIJ + 1 + nI*fld(j-1,iAMR)
         end
      elseif iSibling == 3
         @inbounds for j = 1:nJ
            iCell_G[end,j+1,1] = nIJK*neiBlock -
               nIJ/2 + nI/2 + 1 + nI*fld(j-1,iAMR)
         end
      elseif iSibling == 4
         @inbounds for j = 1:nJ
            iCell_G[end,j+1,1] = nIJK*neiBlock -
               nIJ/2 + 1 + nI*fld(j-1,iAMR)
         end
      elseif iSibling == 6
         @inbounds for j = 1:nJ
            iCell_G[end,j+1,1] = nIJK*(neiBlock-1) +
               nIJ*(nK/2-1) + 1 + nI*fld(j-1,iAMR)
         end
      elseif iSibling == 8
         @inbounds for j = 1:nJ
            iCell_G[end,j+1,1] = nIJK*(neiBlock-1) +
               nIJ*(nK/2-1) + nIJ/2 + 1 + nI*fld(j-1,iAMR)
         end
      end
   end

   # edge 7
   if DiLevelNei_III[1,2,3] == 0
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,2,4], nBlock_P)
   
      @inbounds for j = 1:nJ
         iCell_G[1,j+1,end] = nIJK*(neiBlock-1) + nI*j
      end
   elseif DiLevelNei_III[1,2,3] in (1,2)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,2,4], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      iAMR = 2^DiLevelNei_III[1,2,3]

      if iSibling == 1
         @inbounds for j = 1:nJ
            iCell_G[1,j+1,end] = nIJK*(neiBlock-1) +
               nIJK/2 + nI*(1 + fld(j-1,iAMR))
         end
      elseif iSibling == 3
         @inbounds for j = 1:nJ
            iCell_G[1,j+1,end] = nIJK*(neiBlock-1) +
               nIJK/2 + nIJ/2 + nI*(1 + fld(j-1,iAMR))
         end
      elseif iSibling == 5
         @inbounds for j = 1:nJ
            iCell_G[1,j+1,end] = nIJK*(neiBlock-1) +
               nI*(1 + fld(j-1,iAMR))
         end
      elseif iSibling == 6
         @inbounds for j = 1:nJ
            iCell_G[1,j+1,end] = nIJK*(neiBlock-1) +
               nI*(0.5 + fld(j-1,iAMR))
         end
      elseif iSibling == 7
         @inbounds for j = 1:nJ
            iCell_G[1,j+1,end] = nIJK*(neiBlock-1) +
               nIJ/2 + nI*(1 + fld(j-1,iAMR))
         end
      elseif iSibling == 8
         @inbounds for j = 1:nJ
            iCell_G[1,j+1,end] = nIJK*(neiBlock-1) +
               nIJ/2 + nI*(0.5 + fld(j-1,iAMR))
         end
      end
   end

   # edge 8
   if DiLevelNei_III[3,2,3] == 0
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,2,4], nBlock_P)
      @inbounds for j = 1:nJ
         iCell_G[end,j+1,end] = nIJK*(neiBlock-1) + 1 + nI*(j-1)
      end
   elseif DiLevelNei_III[3,2,3] in (1,2)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,2,4], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      iAMR = 2^DiLevelNei_III[3,2,3]

      if iSibling == 2
         @inbounds for j = 1:nJ
            iCell_G[end,j+1,end] = nIJK*(neiBlock-1) +
               nIJK/2 + 1 + nI*fld(j-1,iAMR)
         end
      elseif iSibling == 4
         @inbounds for j = 1:nJ
            iCell_G[end,j+1,end] = nIJK*(neiBlock-1) +
               nIJK/2 + nIJ/2 + 1 + nI*fld(j-1,iAMR)
         end
      elseif iSibling == 5
         @inbounds for j = 1:nJ
            iCell_G[end,j+1,end] = nIJK*(neiBlock-1) +
               nI/2 + 1 + nI*fld(j-1,iAMR)
         end
      elseif iSibling == 6
         @inbounds for j = 1:nJ
            iCell_G[end,j+1,end] = nIJK*(neiBlock-1) +
               1 + nI*fld(j-1,iAMR)
         end
      elseif iSibling == 7
         @inbounds for j = 1:nJ
            iCell_G[end,j+1,end] = nIJK*(neiBlock-1) +
               nIJ/2 + nI/2 + 1 + nI*fld(j-1,iAMR)
         end
      elseif iSibling == 8
         @inbounds for j = 1:nJ
            iCell_G[end,j+1,end] = nIJK*(neiBlock-1) +
               nIJ/2 + 1 + nI*fld(j-1,iAMR)
         end
      end
   end

   # edge 9
   if DiLevelNei_III[1,1,2] == 0
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,1,2], nBlock_P)

      @inbounds for k = 1:nK
         iCell_G[1,1,k+1] = nIJK*(neiBlock-1) + nIJ*k
      end
   elseif DiLevelNei_III[1,1,2] in (1,2)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,1,2], nBlock_P)

      iAMR = 2^DiLevelNei_III[1,1,2]

      iSibling = getSibling(iNodeNei_III, iTree_IA)
      if iSibling == 1
         @inbounds for k = 1:nK
            iCell_G[1,1,k+1] = nIJK*(neiBlock-1) +
               nIJ*(1 + fld(k-1,iAMR))
         end
      elseif iSibling == 2
         @inbounds for k = 1:nK
            iCell_G[1,1,k+1] = nIJK*(neiBlock-1) +
               nI*(nJ-1) + nI/2 + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 3
         @inbounds for k = 1:nK
            iCell_G[1,1,k+1] = nIJK*(neiBlock-1) +
               nIJ/2 + nIJ*(1 + fld(k-1,iAMR))
         end
      elseif iSibling == 5
         @inbounds for k = 1:nK
            iCell_G[1,1,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nIJ*(1 + fld(k-1,iAMR))
         end
      elseif iSibling == 6
         @inbounds for k = 1:nK
            iCell_G[1,1,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nI*(nJ-1) + nI/2 + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 7
         @inbounds for k = 1:nK
            iCell_G[1,1,k+1] = nIJK*(neiBlock-1) +
               nIJ*nJ/2 + nIJ/2 + nIJ*fld(k-1,iAMR)
         end
      end
   end

   # edge 10
   if DiLevelNei_III[3,1,2] == 0
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,1,2], nBlock_P)

      @inbounds for k = 1:nK
         iCell_G[end,1,k+1] = nIJK*(neiBlock-1) + nI*(nJ-1) + 1 + nIJ*(k-1)
      end
   elseif DiLevelNei_III[3,1,2] in (1,2)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,1,2], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      iAMR = 2^DiLevelNei_III[3,1,2]

      if iSibling == 1
         @inbounds for k = 1:nK
            iCell_G[end,1,k+1] = nIJK*(neiBlock-1) +
               nI*(nJ-1) + nI/2 + 1 + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 2
         @inbounds for k = 1:nK
            iCell_G[end,1,k+1] = nIJK*(neiBlock-1) +
               nI*(nJ-1) + 1 + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 4
         @inbounds for k = 1:nK
            iCell_G[end,1,k+1] = nIJK*(neiBlock-1) +
               nI*(nJ/2-1) + 1 + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 5
         @inbounds for k = 1:nK
            iCell_G[end,1,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nI*(nJ-1) + nI/2 + 1 + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 6
         @inbounds for k = 1:nK
            iCell_G[end,1,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nI*(nJ-1) + 1 + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 8
         @inbounds for k = 1:nK
            iCell_G[end,1,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nI*(nJ/2-1) + 1 + nIJ*fld(k-1,iAMR)
         end
      end
   end

   # edge 11
   if DiLevelNei_III[1,3,2] == 0
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,4,2], nBlock_P)

      @inbounds for k = 1:nK
         iCell_G[1,end,k+1] = nIJK*(neiBlock-1) + nI + nIJ*(k-1)
      end
   elseif DiLevelNei_III[1,3,2] in (1,2)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,4,2], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      iAMR = 2^DiLevelNei_III[1,3,2]

      if iSibling == 1
         @inbounds for k = 1:nK
            iCell_G[1,end,k+1] = nIJK*(neiBlock-1) +
               nIJ/2 + nI + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 3
         @inbounds for k = 1:nK
            iCell_G[1,end,k+1] = nIJK*(neiBlock-1) +
               nI + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 4
         @inbounds for k = 1:nK
            iCell_G[1,end,k+1] = nIJK*(neiBlock-1) +
               nI/2 + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 5
         @inbounds for k = 1:nK
            iCell_G[1,end,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nIJ/2 + nI + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 7
         @inbounds for k = 1:nK
            iCell_G[1,end,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nI + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 8
         @inbounds for k = 1:nK
            iCell_G[1,end,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nI/2 + nIJ*fld(k-1,iAMR)
         end
      end
   end

   # edge 12
   if DiLevelNei_III[3,3,2] == 0
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,4,2], nBlock_P)
      @inbounds for k = 1:nK
         iCell_G[end,end,k+1] = nIJK*(neiBlock-1) + 1 + nIJ*(k-1)
      end
   elseif DiLevelNei_III[3,3,2] in (1,2)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,4,2], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      iAMR = 2^DiLevelNei_III[3,3,2]

      if iSibling == 2
         @inbounds for k = 1:nK
            iCell_G[end,end,k+1] = nIJK*(neiBlock-1) +
               nIJ/2 + 1 + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 3
         @inbounds for k = 1:nK
            iCell_G[end,end,k+1] = nIJK*(neiBlock-1) +
               nI/2 + 1 + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 4
         @inbounds for k = 1:nK
            iCell_G[end,end,k+1] = nIJK*(neiBlock-1) +
               1 + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 6
         @inbounds for k = 1:nK
            iCell_G[end,end,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nIJ/2 + 1 + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 7
         @inbounds for k = 1:nK
            iCell_G[end,end,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + nI/2 + 1 + nIJ*fld(k-1,iAMR)
         end
      elseif iSibling == 8
         @inbounds for k = 1:nK
            iCell_G[end,end,k+1] = nIJK*(neiBlock-1) +
               nIJK/2 + 1 + nIJ*fld(k-1,iAMR)
         end
      end
   end

   ## Node, in total 8

   # node 1
   if DiLevelNei_III[1,1,1] in (0,2,3)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,1,1], nBlock_P)

      iCell_G[1,1,1] = nIJK*neiBlock

   elseif DiLevelNei_III[1,1,1] == 1
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,1,1], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      if iSibling == 1
         iCell_G[1,1,1] = nIJK*neiBlock
      elseif iSibling == 2
         iCell_G[1,1,1] = nIJK*neiBlock - nI/2
      elseif iSibling == 3
         iCell_G[1,1,1] = nIJK*neiBlock - nIJ/2
      elseif iSibling == 4
         iCell_G[1,1,1] = nIJK*neiBlock - nIJ/2 - nI/2
      elseif iSibling == 5
         iCell_G[1,1,1] = nIJK*neiBlock - nIJK/2
      elseif iSibling == 6
         iCell_G[1,1,1] = nIJK*neiBlock - nIJK/2 - nI/2
      elseif iSibling == 7
         iCell_G[1,1,1] = nIJK*neiBlock - nIJK/2 - nIJ/2
      end
   end

   # node 2
   if DiLevelNei_III[3,1,1] in (0,2,3)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,1,1], nBlock_P)

      iCell_G[end,1,1] = nIJK*neiBlock - nI + 1

   elseif DiLevelNei_III[3,1,1] == 1
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,1,1], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      if iSibling == 1
         iCell_G[end,1,1] = nIJK*neiBlock - nI/2 + 1
      elseif iSibling == 2
         iCell_G[end,1,1] = nIJK*neiBlock - nI + 1
      elseif iSibling == 3
         iCell_G[end,1,1] = nIJK*neiBlock - nIJ/2 - nI/2 + 1
      elseif iSibling == 4
         iCell_G[end,1,1] = nIJK*neiBlock - nIJ/2 - nI + 1
      elseif iSibling == 5
         iCell_G[end,1,1] = nIJK*neiBlock - nIJK/2 - nI/2 + 1
      elseif iSibling == 6
         iCell_G[end,1,1] = nIJK*neiBlock - nIJK/2 - nI + 1
      elseif iSibling == 8
         iCell_G[end,1,1] = nIJK*neiBlock - nIJK/2 - nIJ/2 - nI + 1
      end
   end

   # node 3
   if DiLevelNei_III[1,3,1] in (0,2,3)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,4,1], nBlock_P)

      iCell_G[1,end,1] = nIJK*neiBlock - nIJ + nI

   elseif DiLevelNei_III[1,3,1] == 1
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,4,1], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      if iSibling == 1
         iCell_G[1,end,1] = nIJK*neiBlock - nIJ/2 + nI
      elseif iSibling == 2
         iCell_G[1,end,1] = nIJK*neiBlock - nIJ/2 + nI/2
      elseif iSibling == 3
         iCell_G[1,end,1] = nIJK*neiBlock - nIJ + nI
      elseif iSibling == 4
         iCell_G[1,end,1] = nIJK*neiBlock - nIJ + nI/2
      elseif iSibling == 5
         iCell_G[1,end,1] = nIJK*neiBlock - nIJK/2 - nIJ/2 + nI
      elseif iSibling == 7
         iCell_G[1,end,1] = nIJK*neiBlock - nIJK/2 - nIJ + nI
      elseif iSibling == 8
         iCell_G[1,end,1] = nIJK*neiBlock - nIJK/2 - nIJ + nI/2
      end
   end

   # node 4
   if DiLevelNei_III[3,3,1] in (0,2,3)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,4,1], nBlock_P)

      iCell_G[end,end,1] = nIJK*neiBlock - nIJ + 1

   elseif DiLevelNei_III[3,3,1] == 1
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,4,1], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      if iSibling == 1
         iCell_G[end,end,1] = nIJK*neiBlock - nIJ/2 + nI/2 + 1
      elseif iSibling == 2
         iCell_G[end,end,1] = nIJK*neiBlock - nIJ/2 + 1
      elseif iSibling == 3
         iCell_G[end,end,1] = nIJK*neiBlock - nIJ + nI/2 + 1
      elseif iSibling == 4
         iCell_G[end,end,1] = nIJK*neiBlock - nIJ + 1
      elseif iSibling == 6
         iCell_G[end,end,1] = nIJK*neiBlock - nIJK/2 - nIJ/2 + 1
      elseif iSibling == 7
         iCell_G[end,end,1] = nIJK*neiBlock - nIJK/2 - nIJ + nI/2 + 1
      elseif iSibling == 8
         iCell_G[end,end,1] = nIJK*neiBlock - nIJK/2 - nIJ + 1
      end
   end

   # node 5
   if DiLevelNei_III[1,1,3] in (0,2,3)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,1,4], nBlock_P)

      iCell_G[1,1,end] = nIJK*(neiBlock-1) + nIJ

   elseif DiLevelNei_III[1,1,3] == 1
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,1,4], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      if iSibling == 1
         iCell_G[1,1,end] = nIJK*(neiBlock-1) + nIJK/2 + nIJ
      elseif iSibling == 2
         iCell_G[1,1,end] = nIJK*(neiBlock-1) + nIJK/2 + nI*(nJ-1) + nI/2
      elseif iSibling == 3
         iCell_G[1,1,end] = nIJK*(neiBlock-1) + nIJK/2 + nIJ/2
      elseif iSibling == 5
         iCell_G[1,1,end] = nIJK*(neiBlock-1) + nIJ
      elseif iSibling == 6
         iCell_G[1,1,end] = nIJK*(neiBlock-1) + nI*(nJ-1) + nI/2
      elseif iSibling == 7
         iCell_G[1,1,end] = nIJK*(neiBlock-1) + nIJ/2
      elseif iSibling == 8
         iCell_G[1,1,end] = nIJK*(neiBlock-1) + nI*(nJ/2-1) + nI/2
      end
   end

   # node 6
   if DiLevelNei_III[3,1,3] in (0,2,3)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,1,4], nBlock_P)

      iCell_G[end,1,end] = nIJK*(neiBlock-1) + nI*(nJ-1) + 1

   elseif DiLevelNei_III[3,1,3] == 1
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,1,4], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      if iSibling == 1
         iCell_G[end,1,end] = nIJK*(neiBlock-1) + nIJK/2 + nI*(nJ-1) + nI/2 + 1
      elseif iSibling == 2
         iCell_G[end,1,end] = nIJK*(neiBlock-1) + nIJK/2 + nI*(nJ-1) + 1
      elseif iSibling == 4
         iCell_G[end,1,end] = nIJK*(neiBlock-1) + nIJK/2 + nI*(nJ/2-1) + 1
      elseif iSibling == 5
         iCell_G[end,1,end] = nIJK*(neiBlock-1) + nI*(nJ-1) + nI/2 + 1
      elseif iSibling == 6
         iCell_G[end,1,end] = nIJK*(neiBlock-1) + nI*(nJ-1) + 1
      elseif iSibling == 7
         iCell_G[end,1,end] = nIJK*(neiBlock-1) + nI*(nJ/2-1) + nI/2 + 1
      elseif iSibling == 8
         iCell_G[end,1,end] = nIJK*(neiBlock-1) + nI*(nJ/2-1) + 1
      end
   end

   # node 7
   if DiLevelNei_III[1,3,3] in (0,2,3)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,4,4], nBlock_P)

      iCell_G[1,end,end] = nIJK*(neiBlock-1) + nI

   elseif DiLevelNei_III[1,3,3] == 1
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[1,4,4], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      if iSibling == 1
         iCell_G[1,end,end] = nIJK*(neiBlock-1) + nIJK/2 + nIJ/2 + nI
      elseif iSibling == 3
         iCell_G[1,end,end] = nIJK*(neiBlock-1) + nIJK/2 + nI
      elseif iSibling == 4
         iCell_G[1,end,end] = nIJK*(neiBlock-1) + nIJK/2 + nI/2
      elseif iSibling == 5
         iCell_G[1,end,end] = nIJK*(neiBlock-1) + nIJ/2 + nI
      elseif iSibling == 6
         iCell_G[1,end,end] = nIJK*(neiBlock-1) + nIJ/2 + nI/2
      elseif iSibling == 7
         iCell_G[1,end,end] = nIJK*(neiBlock-1) + nI
      elseif iSibling == 8
         iCell_G[1,end,end] = nIJK*(neiBlock-1) + nI/2
      end
   end

   # node 8
   if DiLevelNei_III[3,3,3] in (0,2,3)
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,4,4], nBlock_P)

      iCell_G[end,end,end] = nIJK*(neiBlock-1) + 1

   elseif DiLevelNei_III[3,3,3] == 1
      neiBlock = nodeToGlobalBlock(batl, iNodeNei_III[4,4,4], nBlock_P)

      iSibling = getSibling(iNodeNei_III, iTree_IA)

      # Sibling 1 does not have to compute.
      if iSibling == 2
         iCell_G[end,end,end] = nIJK*(neiBlock-1) + nIJK/2 + nIJ/2 + 1
      elseif iSibling == 3
         iCell_G[end,end,end] = nIJK*(neiBlock-1) + nIJK/2 + nI/2 + 1
      elseif iSibling == 4
         iCell_G[end,end,end] = nIJK*(neiBlock-1) + nIJK/2 + 1
      elseif iSibling == 5
         iCell_G[end,end,end] = nIJK*(neiBlock-1) + nIJ/2 + nI/2 + 1
      elseif iSibling == 6
         iCell_G[end,end,end] = nIJK*(neiBlock-1) + nIJ/2 + 1
      elseif iSibling == 7
         iCell_G[end,end,end] = nIJK*(neiBlock-1) + nI/2 + 1
      elseif iSibling == 8
         iCell_G[end,end,end] = nIJK*(neiBlock-1) + 1
      end
   end

end