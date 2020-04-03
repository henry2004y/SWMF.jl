# Convert full BATSRUS AMR output to VTK.

using FortranFiles, SWMF, WriteVTK

##

# Number of cells per block in each direction.
# These values are set by the Config.pl script.
# Set 1 for ignored directions!
const nI, nJ, nK = 4, 4, 4

# Maximum number of ghost cells set by Config.pl script.
# Valid values are 0,1,2,3,4,5
const nG = 2

# Refinement ratios in the 3 dimensions. Either 1 or 2.
# The values are set by the Config.pl script.
const iRatio, jRatio, kRatio = min(2, nI), min(2, nJ), min(2, nK)

# Number of dimensions in which grid adaptation is done
const nDimAmr = iRatio + jRatio + kRatio - 3

# Number of children per node
const nChild = 2^nDimAmr

# Possible values for the status variable
const Unset_     = -100 # index for unset values (that are otherwise larger)
const Unused_      = -1 # unused block (not a leaf)
const Refine_      = -2 # parent block to be refined
const DontCoarsen  = -3 # block not to be coarsened
const Coarsen_     = -4 # child block to be coarsened
const Used_        =  1 # currently used block (leaf)
const RefineNew_   =  2 # child block to be refined
const Refined_     =  3 # refined child block
const CoarsenNew_  =  4 # parent block to be coarsened
const Coarsened_   =  5 # coarsened parent block

# Deepest AMR level relative to root nodes (limited by 32 bit integers)
const MaxLevel = 30

# Named indexes of iTree_IA
const Status_   =  1
const Level_    =  2 # grid level
const Proc_     =  3 # processor index
const Block_    =  4 # block index
const MinLevel_ =  5 # minimum level allowed
const MaxLevel_ =  6 # maximum level allowed
const Coord0_   =  6 # equal to Coord1_-1
const Coord1_   =  7 # coordinate of node in 1st dimension
const Coord2_   =  8 # coordinate of node in 2nd dimension
const Coord3_   =  9 # coordinate of node in 3rd dimension
const CoordLast_=  9 # Coord0_ + MaxDim (?)
const Parent_   = 10 # Parent_ must be equal to Child0_
const Child0_   = 10 #
const Child1_   = Child0_ + 1
const ChildLast_= Child0_ + nChild

#=
the status of the node (used, unused, to be refined, to be coarsened, etc.);
the current, the maximum allowed and minimum allowed AMR levels for this node;
the three integer coordinates with respect to the whole grid;
the index of the parent node (if any);
the indexes of the children nodes (if any);
the processor index where the block is stored for active nodes;
the local block index for active nodes.
=#

filetag = "3d__ful_2_t00000141_n00000142"

## Loading AMR tree
f = FortranFile(filetag*".tree")

nDim, nInfo, nNode = read(f, Int32, Int32, Int32)
iRatio_D = read(f, (Int32,nDim)) # Array of refinement ratios
nRoot_D = read(f, (Int32,nDim)) # The number of root nodes in all dimension
iTree_IA = read(f, (Int32,(nInfo,nNode)))

close(f)

# Get all the used blocks
blockused_ = [iNode for iNode in 1:nNode if iTree_IA[Status_,iNode] == 1]

## If you know how blocks are organized, you can get VTK!
data = readdata(filetag*".out")

# vtkMultiBlockDataSet
vtmfile = vtk_multiblock("vtm/my_vtm_file")

for (i,iblock) in enumerate(blockused_)
   x = data.x[64*(i-1)+1:64*(i-1)+4,:,:,1][:]
   y = data.x[64*(i-1)+1:4:64*(i-1)+13,:,:,2][:]
   z = data.x[64*(i-1)+1:16:64*(i-1)+49,:,:,3][:]
   p = data.w[64*(i-1)+1:64*(i-1)+64,:,:,8][:]

   vtkfile = vtk_grid(vtmfile, x, y, z)
   vtkfile["Pressure",VTKPointData()] = p
end

outfiles = vtk_save(vtmfile)

# vtkHierarchicalBoxDataSet
