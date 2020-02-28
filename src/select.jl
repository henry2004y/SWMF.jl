export get_vars, cutdata, subvolume, subsurface


"""
	cutdata(data, head, var; plotrange=[-Inf,Inf,-Inf,Inf], cut=' ',
		cutPlaneIndex=1)

Get 2D plane cut data of 3D box data.
"""
function cutdata(data::Data, head::Dict, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], cut=' ', cutPlaneIndex=1)

   x,w = data.x, data.w
   VarIndex_ = findfirst(x->x==lowercase(var), lowercase.(head[:wnames]))
   isempty(VarIndex_) && error("$(var) not found in header variables!")

   X = @view x[:,:,:,1]
   Y = @view x[:,:,:,2]
   Z = @view x[:,:,:,3]

   W = w[:,:,:,VarIndex_]

   if cut ∈ ('x',' ')
      cut1 = @view X[cutPlaneIndex,:,:]
      cut2 = @view Y[cutPlaneIndex,:,:]
      W    = @view W[cutPlaneIndex,:,:]
   elseif cut ==  'y'
      cut1 = @view X[:,cutPlaneIndex,:]
      cut2 = @view Z[:,cutPlaneIndex,:]
      W    = @view W[:,cutPlaneIndex,:]
   elseif cut == 'z'
      cut1 = @view X[:,:,cutPlaneIndex]
      cut2 = @view Y[:,:,cutPlaneIndex]
      W    = @view W[:,:,cutPlaneIndex]
   end

   if !all(isinf.(plotrange))
      cut1, cut2, W = subsurface(cut1, cut2, W, plotrange)
   end

   return cut1, cut2, W
end

"""
	subsurface(x, y, data, limits)
	subsurface(x, y, u, v, limits)

Extract subset of 2D surface dataset.
This is a simplified version of subvolume.
"""
function subsurface(x, y, data, limits)

   if length(limits)!=4
      @error "Reduction must be [xmin xmax ymin ymax]"
   end

   if limits[1] > limits[2] || limits[3] > limits[4]
      @error "subsurface:InvalidReductionXRange"
   end

   sz = size(data)

   hx = x[:,1]
   hy = y[1,:]

   if isinf(limits[1]) limits[1] = minimum(hx) end
   if isinf(limits[3]) limits[3] = minimum(hy) end
   if isinf(limits[2]) limits[2] = maximum(hx) end
   if isinf(limits[4]) limits[4] = maximum(hy) end

   xind = findall(limits[1] .≤ hx .≤ limits[2])
   yind = findall(limits[3] .≤ hy .≤ limits[4])

   newdata = subdata(data, xind, yind, sz)

   newx = x[xind, yind]
   newy = y[xind, yind]

   return newx, newy, newdata
end

function subsurface(x, y, u, v, limits)

   if length(limits)!=4
      @error "Reduction must be [xmin xmax ymin ymax]"
   end

   if limits[1] > limits[2] || limits[3] > limits[4]
      @error "subsurface:InvalidReductionXRange"
   end

   sz = size(u)

   hx = x[:,1]
   hy = y[1,:]

   if isinf(limits[1]) limits[1] = minimum(hx) end
   if isinf(limits[3]) limits[3] = minimum(hy) end
   if isinf(limits[2]) limits[2] = maximum(hx) end
   if isinf(limits[4]) limits[4] = maximum(hy) end

   xind = findall(limits[1] .≤ hx .≤ limits[2])
   yind = findall(limits[3] .≤ hy .≤ limits[4])

   newu = subdata(u, xind, yind, sz)
   newv = subdata(v, xind, yind, sz)

   newx = x[xind, yind]
   newy = y[xind, yind]

   return newx, newy, newu, newv
end

"""
	subvolume(x, y, z, data, limits)
	subvolume(x, y, z, u, v, w, limits)

Extract subset of 3D dataset in ndgrid format.
"""
function subvolume(x, y, z, data, limits)
   if length(limits)!=6
      @error "Reduction must be [xmin xmax ymin ymax zmin zmax]"
   end

   if limits[1] > limits[2] || limits[3] > limits[4] || limits[5] > limits[6]
      @error "subvolume:InvalidReductionXRange"
   end

   sz = size(data)

   hx = x[:,1,1]
   hy = y[1,:,1]
   hz = z[1,1,:]

   if isinf(limits[1]) limits[1] = minimum(hx) end
   if isinf(limits[3]) limits[3] = minimum(hy) end
   if isinf(limits[5]) limits[5] = minimum(hz) end
   if isinf(limits[2]) limits[2] = maximum(hx) end
   if isinf(limits[4]) limits[4] = maximum(hy) end
   if isinf(limits[6]) limits[6] = maximum(hz) end

   xind = findall(limits[1] .≤ hx .≤ limits[2])
   yind = findall(limits[3] .≤ hy .≤ limits[4])
   zind = findall(limits[5] .≤ hy .≤ limits[6])

   newdata = subdata(data, xind, yind, zind, sz)

   newx = x[xind,yind,zind]
   newy = y[xind,yind,zind]
   newz = z[xind,yind,zind]

   return newx, newy, newz, newdata
end

function subvolume(x, y, z, u, v, w, limits)
   if length(limits)!=6
      @error "Reduction must be [xmin xmax ymin ymax zmin zmax]"
   end

   if limits[1] > limits[2] || limits[3] > limits[4] || limits[5] > limits[6]
      @error "subvolume:InvalidReductionXRange"
   end

   sz = size(u)

   hx = x[:,1,1]
   hy = y[1,:,1]
   hz = z[1,1,:]

   if isinf(limits[1]) limits[1] = minimum(hx) end
   if isinf(limits[3]) limits[3] = minimum(hy) end
   if isinf(limits[5]) limits[5] = minimum(hz) end
   if isinf(limits[2]) limits[2] = maximum(hx) end
   if isinf(limits[4]) limits[4] = maximum(hy) end
   if isinf(limits[6]) limits[6] = maximum(hz) end

   xind = findall(limits[1] .≤ hx .≤ limits[2])
   yind = findall(limits[3] .≤ hy .≤ limits[4])
   zind = findall(limits[5] .≤ hz .≤ limits[6])

   newu = subdata(u, xind, yind, zind, sz)
   newv = subdata(v, xind, yind, zind, sz)
   neww = subdata(w, xind, yind, zind, sz)

   newx = x[xind,yind,zind]
   newy = y[xind,yind,zind]
   newz = z[xind,yind,zind]

   return newx, newy, newz, newu, newv, neww
end

"""
	subdata(data, xind, yind, sz)
	subdata(data, xind, yind, zind, sz)

Return the sliced data based on indexes.
"""
function subdata(data, xind::Vector{Int}, yind::Vector{Int}, sz::Tuple{Int,Int})
   newdata = data[xind,yind]
   newsz = size(newdata)

   if length(sz) > 2
      newdata = reshape(newdata, (newsz[1:2]..., sz[3:end]))
   end

   return newdata
end

function subdata(data, xind::Vector{Int}, yind::Vector{Int}, zind::Vector{Int},
   sz::Tuple{Int,Int,Int})

   newdata = data[xind,yind,zind]
   newsz = size(newdata)

   if length(sz) > 3
      newdata = reshape(newdata, (newsz[1:3]..., sz[4:end]))
   end

   return newdata
end


function get_var(data::Data, head::Dict, var::AbstractString)
   VarIndex_ = findfirst(x->x==var,head[:wnames])

   ndim = head[:ndim]
   if ndim == 1
      w = data.w[:,VarIndex_]
   elseif ndim == 2
      w = data.w[:,:,VarIndex_]
   elseif ndim == 3
      w = data.w[:,:,:,VarIndex_]
   end
   w
end

function get_vars(data::Data, head::Dict, Names::Vector{AbstractString})

   dict = Dict()
   for name in Names
      dict[name] = get_var(data, head, name)
   end

   Vars(dict)
end

Base.getproperty(p::Vars, name::Symbol) = getfield(p, :data)[String(name)]
