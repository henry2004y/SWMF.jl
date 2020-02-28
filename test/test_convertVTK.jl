# Convert SWMF GM outputs to vtk files.
#
# Hongyang Zhou, hyzhou@umich.edu 07/23/2019

using SWMF, WriteVTK
try
   using MATLAB
catch
   println("MATLAB not found path...")
end

searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

#mypath = "."
#mykey  = "box_test.out"

#filenames = searchdir(mypath,mykey)

function convertBox2VTK_matlab(filenames::Array{String,1})

   for filename in filenames
      filehead, data = mxcall(:read_data,2,filename)

      data = data["file1"]

      x = data["x"][:,:,:,1]
      y = data["x"][:,:,:,2]
      z = data["x"][:,:,:,3]

      func = "Bx"
      func_ = findfirst(x->x==func, filehead["wnames"])
      if isnothing(func_) @error "Couldn't find variable $(func)!" end
      B = data["w"][:,:,:,func_:func_+2]
      B = permutedims(B, [4,1,2,3])

      func = "P"
      func_ = findfirst(x->x==func, filehead["wnames"])
      if isnothing(func_) @error "Couldn't find variable $(func)!" end
      P = data["w"][:,:,:,func_]

      outname = filename[1:end-4]

      # Rectilinear or structured grid
      outfiles = vtk_grid(outname, x, y, z) do vtk
         vtk_point_data(vtk, P, "P")
         vtk_point_data(vtk, B, "B")
      end
   end

end

function convertBox2VTK(filenames::Array{String,1}; dir=".", gridType=1)

   for filename in filenames
      head, data, list = readdata(filename, dir=dir, verbose=false)

      nVar = length(head[1][:wnames])

      outname = filename[1:end-4]

      if gridType == 1 # rectilinear grid
         x = @view data[1].x[:,1,1,1]
         y = @view data[1].x[1,:,1,2]
         z = @view data[1].x[1,1,:,3]

         outfiles = vtk_grid(dir*outname, x,y,z) do vtk
            for ivar = 1:nVar
               if head[1][:wnames][ivar][end] == 'x' # vector
                  var1 = @view data[1].w[:,:,:,ivar]
                  var2 = @view data[1].w[:,:,:,ivar+1]
                  var3 = @view data[1].w[:,:,:,ivar+2]
                  namevar = head[1][:wnames][ivar][1:end-1]
                  vtk_point_data(vtk, (var1, var2, var3), namevar)
               elseif head[1][:wnames][ivar][end] in ('y','z')
                  continue
               else
                  var = @view data[1].w[:,:,:,ivar]
                  vtk_point_data(vtk, var, head[1][:wnames][ivar])
               end
            end
         end
      elseif gridType == 2 # structured grid
         xyz = permutedims(data[1].x, [4,1,2,3])

         outfiles = vtk_grid(dir*outname, xyz) do vtk
            for ivar = 1:nVar
               if head[1][:wnames][ivar][end] == 'x' # vector
                  var1 = @view data[1].w[:,:,:,ivar]
                  var2 = @view data[1].w[:,:,:,ivar+1]
                  var3 = @view data[1].w[:,:,:,ivar+2]
                  namevar = head[1][:wnames][ivar][1:end-1]
                  vtk_point_data(vtk, (var1, var2, var3), namevar)
               elseif head[1][:wnames][ivar][end] in ('y','z')
                  continue
               else
                  var = @view data[1].w[:,:,:,ivar]
                  vtk_point_data(vtk, var, head[1][:wnames][ivar])
               end
            end
         end
      elseif gridType == 3 # unstructured grid, not finished
         vtkfile = vtk_grid(dir*outname, points, cells)
         vtk_cell_data(vtkfile, cdata, "my_cell_data")
         outfiles = vtk_save(vtkfile)
      end
      println(filename," finished conversion.")
   end

end

function convertTec2VTK()
   filename = "3d_ascii.dat"
   head, data, connectivity  = readtecdata(filename,false)

   points = @view data[1:3,:]
   cells = Vector{MeshCell{Array{Int32,1}}}(undef,head[:nCell])
   if head[:ndim] == 3
      @inbounds for i = 1:head[:nCell]
         cells[i] = MeshCell(VTKCellTypes.VTK_VOXEL, connectivity[:,i])
      end
   elseif head[:ndim] == 2
      @inbounds for i = 1:head[:nCell]
         cells[i] = MeshCell(VTKCellTypes.VTK_PIXEL, connectivity[:,i])
      end
   end

   # 777MB in Ascii to 147MB in VTK binary;
   # What about preplot? 234MB
   vtkfile = vtk_grid("test_unstructured1", points, cells)

   rho = @view data[4,:]
   p = @view data[14,:]

   vtk_point_data(vtkfile, rho, "Rho")
   vtk_point_data(vtkfile, p, "P")

   outfiles = vtk_save(vtkfile)
end

function test_bin()
   filename = "3d_bin.dat"

   head, data, connectivity  = readtecdata(filename,true)

   points = @view data[1:3,:]
   cells = Vector{MeshCell{Array{Int32,1}}}(undef,head[:nCell])
   if head[:ndim] == 3
      @inbounds for i = 1:head[:nCell]
         cells[i] = MeshCell(VTKCellTypes.VTK_VOXEL, connectivity[:,i])
      end
   elseif head[:ndim] == 2
      @inbounds for i = 1:head[:nCell]
         cells[i] = MeshCell(VTKCellTypes.VTK_PIXEL, connectivity[:,i])
      end
   end

   vtkfile = vtk_grid("test_unstructured2", points, cells)

   rho = @view data[4,:]
   p = @view data[14,:]

   vtk_point_data(vtkfile, rho, "Rho")
   vtk_point_data(vtkfile, p, "P")

   outfiles = vtk_save(vtkfile)

end


function test_cell()
   filename = "3d.tcp"

   head, data, connectivity  = readtecdata(filename,true)

   nVar = length(head[:variables])

   points = @view data[1:head[:ndim],:]
   cells = Vector{MeshCell{Array{Int32,1}}}(undef,head[:nCell])
   if head[:ndim] == 3
      # PLT to VTK index_ = [1 2 4 3 5 6 8 7]
      for i = 1:2
         connectivity = swaprows(connectivity, 4*i-1, 4*i)
      end
      @inbounds for i = 1:head[:nCell]
         cells[i] = MeshCell(VTKCellTypes.VTK_VOXEL, connectivity[:,i])
      end
   elseif head[:ndim] == 2
      @inbounds for i = 1:head[:nCell]
         cells[i] = MeshCell(VTKCellTypes.VTK_PIXEL, connectivity[:,i])
      end
   end

   vtkfile = vtk_grid("test_unstructured2", points, cells)

   for ivar = head[:ndim]+1:nVar
      if occursin("_x",head[:variables][ivar]) # vector
         var1 = @view data[ivar,:]
         var2 = @view data[ivar+1,:]
         var3 = @view data[ivar+2,:]
         namevar = replace(head[:variables][ivar], "_x"=>"")
         vtk_point_data(vtkfile, (var1, var2, var3), namevar)
      elseif occursin(r"(_y|_z)",head[:variables][ivar])
         continue
      else
         var = @view data[ivar,:]
         vtk_point_data(vtkfile, var, head[:variables][ivar])
      end
   end

   outfiles = vtk_save(vtkfile)

end

function swaprows(X, i, j)
   m, n = size(X)
   if (1 <= i <= n) && (1 <= j <= n)
      for k = 1:n
        @inbounds X[i,k],X[j,k] = X[j,k],X[i,k]
      end
      return X
   else
      throw(BoundsError())
   end
end


function test_ascii_convert()
   filename = "3d_ascii.dat"
   head, data, connectivity  = readtecdata(filename,false)
   convertVTK(head, data, connectivity)

end
