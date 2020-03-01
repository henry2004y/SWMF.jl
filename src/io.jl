# All the IO related APIs.

export readdata, readlogdata, readtecdata, convertVTK

searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

"""
	readdata(filenameIn, (, dir=".", npict=1, verbose=false))

Read data from BATSRUS output files. Stores the npict-th snapshot from an ascii
or binary data file into the coordinates `x` and data `w` arrays.
Filenames can be provided with wildcards.

# Examples
```jldoctest
filename = "1d_raw*"
data = readdata(filename)
```
"""
function readdata(filenameIn::AbstractString; dir=".", npict=1, verbose=false)

   # Check the existence of files
	filenames = searchdir(dir, Regex(filenameIn)) # potential bugs
   if isempty(filenames)
      @error "readdata: no matching filename was found for $(filenameIn)"
   elseif length(filenames) > 1
      @error "Ambiguous filenames!"
   end

	filename = joinpath(dir, filenames[1])
   filelist, fileID, pictsize = getFileType(filename)

   if verbose
      @info "filename=$(filelist.name)\n"*"npict=$(filelist.npictinfiles)"
   end

   if any(filelist.npictinfiles - npict < 0)
      @error "file $(ifile): npict out of range!"
   end
   seekstart(fileID) # Rewind to start

   ## Read data from files

   # Skip npict-1 snapshots (because we only want npict snapshot)
   skip(fileID, pictsize*(npict-1))

   filehead = getfilehead(fileID, filelist.type)

   # Read data
   fileType = lowercase(filelist.type)
   if fileType == "ascii"
      x,w = getpictascii(fileID, filehead)
   elseif fileType == "binary"
      x,w = getpictreal(fileID, filehead, Float64)
   elseif fileType == "real4"
      x,w = getpictreal(fileID, filehead, Float32)
   else
      @error "get_pict: unknown filetype: $(filelist[ifile].type)"
   end

   #setunits(filehead,"")

   verbose && showhead(filelist, filehead)

	data = Data(filehead, x, w, filelist)

   verbose && @info "Finished reading $(filelist.name)"

   close(fileID)

   return data
end

"Read information from log file."
function readlogdata( filename::AbstractString )

   f = open(filename, "r")
   nLine = countlines(f) - 2
   seekstart(f)
   headline  = readline(f)
   variables = split(readline(f))
   ndim      = 1
   it        = 0
   t         = 0.0
   gencoord  = false
   nx        = 1
   nw        = length(variables)

   data = zeros(nw,nLine)
   for i = 1:nLine
      line = split(readline(f))
      data[:,i] = parse.(Float64,line)
   end

   close(f)

   head = (ndim=ndim, headline=headline, it=it, time=t, gencoord=gencoord,
      neqpar=neqpar, nw=nw, nx=nx, variables=variables)

   return head, data
end

"""
	readtecdata(filename, IsBinary=false, verbose=false)

Return header, data and connectivity from BATSRUS Tecplot outputs. Both binary
and ASCII formats are supported. The default is reading pure ASCII data.
# Examples
```jldoctest
filename = "3d_ascii.dat"
head, data, connectivity = readtecdata(filename)
```
"""
function readtecdata(filename::AbstractString; IsBinary=false, verbose=false)

   f = open(filename)

   nNode = Int32(0)
   nCell = Int32(0)
   ET = ""

   # Read Tecplot header
   ln = readline(f) |> strip
   if startswith(ln, "TITLE")
      title = match(r"\"(.*?)\"", split(ln,'=', keepempty=false)[2])[1]
   else
      @warn "No title provided."
   end
   ln = readline(f) |> strip
   if startswith(ln, "VARIABLES")
      # Read until another keyword appears
      varline = split(ln,'=')[2]

      ln = readline(f) |> strip
      while !startswith(ln, "ZONE")
         varline *= strip(ln)
         ln = readline(f) |> strip
      end
      VARS = split(varline, '\"')
      deleteat!(VARS,findall(x -> x in (" ",", ") || isempty(x), VARS))
   else
      @warn "No variable names provided."
   end

   if startswith(ln, "ZONE")
      ndim = 3 # default is 3D

      pt0 = position(f)

      while !isnothing(findfirst('=', ln))
         if startswith(ln, "AUXDATA") || startswith(ln,"DT")
            pt0 = position(f)
            # Are there better ways to identify end of header?
            if IsBinary && startswith(ln, "AUXDATA TIMESIMSHORT")
               break # last line of AUXDATA
            end
            ln = readline(f) |> strip # I should save the auxdata!
            continue
         end

         if !startswith(ln, "ZONE")
            zoneline = split(ln, ", ", keepempty=false)
         else # if the ZONE line has nothing, this won't work!
            zoneline = split(ln[6:end], ", ", keepempty=false)
            replace(zoneline[1], '"'=>"") # Remove the quotes in T
         end
         for zline in zoneline
            name, value = split(zline,'=', keepempty=false)
            name = uppercase(name)
            if name == "T" # ZONE title
               T = value
            elseif name in ("NODES","N")
               try
                  nNode = parse(Int32, value)
               catch
                  IsBinary = true
                  nNode = read(IOBuffer(value), Int32)
               end
            elseif name in ("ELEMENTS","E")
               try
                  nCell = parse(Int32, value)
               catch
                  IsBinary = true
                  nCell = read(IOBuffer(value), Int32)
               end
            elseif name in ("ET","ZONETYPE")
               if uppercase(value) in ("BRICK","FEBRICK")
                  ndim = 3
               elseif uppercase(value) == "FEQUADRILATERAL"
                  ndim = 2
               end
               ET = uppercase(value)
            end
         end
         ln = readline(f) |> strip
      end
   else
      @warn "No zone info provided."
   end

   seek(f, pt0)

   data = Array{Float32,2}(undef, length(VARS), nNode)

   if ndim == 3
	   connectivity = Array{Int32,2}(undef,8,nCell)
   elseif ndim == 2
	   connectivity = Array{Int32,2}(undef,4,nCell)
   end

   if IsBinary
	   @inbounds for i = 1:nNode
		   read!(f, @view data[:,i])
	   end
	   @inbounds for i = 1:nCell
         read!(f, @view connectivity[:,i])
      end
   else
   	@inbounds for i = 1:nNode
         x = readline(f)
         data[:,i] .= parse.(Float32, split(x))
      end
	   @inbounds for i = 1:nCell
		   x = readline(f)
		   connectivity[:,i] .= parse.(Int32, split(x))
	   end
   end

   close(f)

   head = (variables=VARS, nNode=nNode, nCell=nCell, ndim=ndim, ET=ET,
		title=title)

   return head, data, connectivity
end


"Obtain file type."
function getFileType(filename)

   fileID = open(filename, "r")
   bytes = filesize(filename)
   type  = ""

   # Check the appendix of file names
   # Gabor uses a trick: the first 4 bytes decides the file type
   if occursin(r"^.*\.(log)$", filename)
      type = "log"
      npictinfiles = 1
   elseif occursin(r"^.*\.(dat)$", filename)
      # Tecplot ascii format
      type = "dat"
      npictinfiles = 1
   else
      # Obtain filetype based on the length info in the first 4 bytes
      lenhead = read(fileID, Int32)

      if lenhead!=79 && lenhead!=500
         type = "ascii"
      else
         # The length of the 2nd line decides between real4 & real8
         # since it contains the time; which is real*8 | real*4
         skip(fileID, lenhead+4)
         len = read(fileID, Int32)
         if len == 20
            type = "real4"
         elseif len == 24
            type = "binary"
         else
            @error "Error in getFileTypes: strange unformatted file:
               $(filename)"
         end

         if lenhead == 500
            type = uppercase(type)
         end
      end
         # Obtain file size & number of snapshots
         seekstart(fileID)
         pictsize = getfilesize(fileID, type)
         npictinfiles = floor(Int, bytes / pictsize)
      end

      filelist = FileList(filename, type, bytes, npictinfiles)

   return filelist, fileID, pictsize
end

"""
	getfilehead(fileID, type, iargout=1)

Obtain the header information from BATSRUS output files.
# Input arguments
- `fileID::IOStream`: file identifier.
- `type::String`: file type in ["ascii", "real4", "binary", "log"].
- `iargout::Int`: 1 for output pictsize, 2 for output filehead.
# Output arguments
- `filehead::NamedTuple`: file header info.
"""
function getfilehead(fileID::IOStream, type::String)

   ftype = string(lowercase(type))

   if ftype == type lenstr = 79 else lenstr = 500 end

   # Read header
   pointer0 = position(fileID)

   if ftype == "ascii"
      headline = readline(fileID)
      line = readline(fileID)
      line = split(line)
      it = parse(Int,line[1])
      t = parse(Float64,line[2])
      ndim = parse(Int8,line[3])
      neqpar = parse(Int32,line[4])
      nw = parse(Int8,line[5])
      gencoord = ndim < 0
      ndim = abs(ndim)
      nx = parse.(Int64, split(readline(fileID)))
      if neqpar > 0
         eqpar = parse.(Float64, split(readline(fileID)))
      end
      varname = readline(fileID)
   elseif ftype ∈ ["real4","binary"]
      skip(fileID,4) # skip record start tag.
      headline = String(read(fileID, lenstr))
      skip(fileID,8) # skip record end/start tags.
      it = read(fileID, Int32)
      t = read(fileID, Float32)
      ndim = read(fileID, Int32)
      gencoord = (ndim < 0)
      ndim = abs(ndim)
      neqpar = read(fileID, Int32)
      nw = read(fileID, Int32)
      skip(fileID,8) # skip record end/start tags.
      nx = zeros(Int32, ndim)
      read!(fileID, nx)
      skip(fileID,8) # skip record end/start tags.
      if neqpar > 0
         eqpar = zeros(Float32,neqpar)
         read!(fileID, eqpar)
         skip(fileID,8) # skip record end/start tags.
      end
      varname = String(read(fileID, lenstr))
      skip(fileID,4) # skip record end tag.
   end

   # Header length
   pointer1 = position(fileID)
   headlen = pointer1 - pointer0

   # Calculate the snapshot size = header + data + recordmarks
   nxs = prod(nx)

   # Set variables array
   variables = split(varname) # returns a string array

	# Produce a wnames from the last file
   wnames = variables[ndim+1:ndim+nw]

	head = (ndim=ndim, headline=headline, it=it, time=t, gencoord=gencoord,
		neqpar=neqpar, nw=nw, nx=nx, eqpar=eqpar, variables=variables,
		wnames=wnames)
end

""" Return the size of file. """
function getfilesize(fileID::IOStream, type::String)

   pictsize = 0

   ndim = convert(Int32,1)
   tmp  = convert(Int32,1)
   nw   = convert(Int32,1)
   nxs  = 0

   ftype = string(lowercase(type))

   if ftype == type lenstr = 79 else lenstr = 500 end

   # Read header
   pointer0 = position(fileID)

   if ftype == "ascii"
      headline = readline(fileID)
      line = readline(fileID)
      line = split(line)
      it = parse(Int,line[1])
      time = parse(Float64,line[2])
      ndim = parse(Int8,line[3])
      neqpar = parse(Int32,line[4])
      nw = parse(Int8,line[5])
      gencoord = ndim < 0
      ndim = abs(ndim)
      nx = parse.(Int64, split(readline(fileID)))
      if neqpar > 0
         eqpar = parse.(Float64, split(readline(fileID)))
      end
      varname = readline(fileID)
   elseif ftype ∈ ["real4","binary"]
      skip(fileID,4)
      read(fileID,lenstr)
      skip(fileID,8)
      read(fileID,Int32)
      read(fileID,Float32)
      ndim = abs(read(fileID,Int32))
      tmp = read(fileID,Int32)
      nw = read(fileID,Int32)
      skip(fileID,8)
      nx = zeros(Int32,ndim)
      read!(fileID,nx)
      skip(fileID,8)
      if tmp > 0
         tmp2 = zeros(Float32,tmp)
         read!(fileID,tmp2)
         skip(fileID,8) # skip record end/start tags.
      end
      read(fileID, lenstr)
      skip(fileID,4)
   end

   # Header length
   pointer1 = position(fileID)
   headlen = pointer1 - pointer0

   # Calculate the snapshot size = header + data + recordmarks
   nxs = prod(nx)

   if ftype == "log"
      pictsize = 1
   elseif ftype == "ascii"
      pictsize = headlen + (18*(ndim+nw)+1)*nxs
   elseif ftype == "binary"
      pictsize = headlen + 8*(1+nw) + 8*(ndim+nw)*nxs
   elseif ftype == "real4"
      pictsize = headlen + 8*(1+nw) + 4*(ndim+nw)*nxs
   end

   return pictsize
end

# This will be supported in Julia v1.4. Remove it then.
import Base: read!

"""
	read!(s,a)

Read slices of arrays using subarrays, in addition to the built-in methods.
"""
function read!(s::IO, a::AbstractArray{T}) where T
   GC.@preserve a unsafe_read(s, pointer(a), sizeof(a))
   return a
end


"""
	getpictascii(fileID, filehead)

Read ascii format data.
"""
function getpictascii(fileID::IOStream, filehead::NamedTuple)

   ndim = filehead.ndim
   nw   = filehead.nw

   # Read coordinates & values row by row
   if ndim == 1 # 1D
      n1 = filehead.nx[1]
      x  = Array{Float64,2}(undef,n1,ndim)
      w  = Array{Float64,2}(undef,n1,nw)
      for ix = 1:n1
         temp = parse.(Float64, split(readline(fileID)))
         x[ix,:] .= temp[1]
         w[ix,:] .= temp[2:end]
      end
   elseif ndim == 2 # 2D
      n1, n2 = filehead.nx
      x  = Array{Float64,3}(undef,n1,n2,ndim)
      w  = Array{Float64,3}(undef,n1,n2,nw)
      for i = 1:n1, j = 1:n2
         temp = parse.(Float64, split(readline(fileID)))
         x[i,j,:] .= temp[1:2]
         w[i,j,:] .= temp[3:end]
      end
   elseif ndim == 3 # 3D
      n1, n2, n3 = filehead.nx
      x  = Array{Float64,4}(undef,n1,n2,n3,ndim)
      w  = Array{Float64,4}(undef,n1,n2,n3,nw)
      for i = 1:n1, j = 1:n2, k = 1:n3
         temp = parse.(Float64, split(readline(fileID)))
         x[i,j,k,:] .= temp[1:3]
         w[i,j,k,:] .= temp[4:end]
      end
   end

   return x, w
end


"""
	getpictreal(fileID, filehead)

Read real4/read8 format data.
"""
function getpictreal(fileID::IOStream, filehead::NamedTuple, T::DataType)

   ndim = filehead.ndim
   nw   = filehead.nw

   # Read coordinates & values
   if ndim == 1 # 1D
      n1 = filehead.nx[1]
      x  = Array{T,2}(undef,n1,ndim)
      w  = Array{T,2}(undef,n1,nw)
      skip(fileID,4) # skip record start tag.
      read!(fileID,x)
      skip(fileID,8) # skip record end/start tags.
      for iw = 1:nw
         read!(fileID, @view w[:,iw])
         skip(fileID,8) # skip record end/start tags.
      end
   elseif ndim == 2 # 2D
      n1, n2 = filehead.nx
      x  = Array{T,3}(undef,n1,n2,ndim)
      w  = Array{T,3}(undef,n1,n2,nw)
      skip(fileID,4) # skip record start tag.
      read!(fileID,x)
      skip(fileID,8) # skip record end/start tags.
      for iw = 1:nw
         read!(fileID, @view w[:,:,iw])
         skip(fileID,8) # skip record end/start tags.
      end
   elseif ndim == 3 # 3D
      n1, n2, n3 = filehead.nx
      x  = Array{T,4}(undef,n1,n2,n3,ndim)
      w  = Array{T,4}(undef,n1,n2,n3,nw)
      skip(fileID,4) # skip record start tag.
      read!(fileID,x)
      skip(fileID,8) # skip record end/start tags.
      for iw = 1:nw
         read!(fileID, @view w[:,:,:,iw])
         skip(fileID,8) # skip record end/start tags.
      end
   end

   return x,w
end

"""
	setunits(filehead, type, (distunit, Mion, Melectron))

Set the units for the output files。
If type is given as "SI", "CGS", "NORMALIZED", "PIC", "PLANETARY", "SOLAR", set
`typeunit = type`, otherwise try to guess from the fileheader.
Based on `typeunit` set units for distance [xSI], time [tSI], density [rhoSI],
pressure [pSI], magnetic field [bSI] and current density [jSI] in SI units.
Distance unit [rplanet | rstar], ion & electron mass in [amu] can be set with
optional distunit, Mion and Melectron.

Also calculate convenient constants ti0, cs0 ... for typical formulas.
This function needs to be improved!
"""
function setunits( filehead::NamedTuple, type::AbstractString; distunit=1.0,
   Mion=1.0, Melectron=1.0)

   # This is currently not used, so return here
   return

   ndim      = filehead.ndim
   headline  = filehead.headline
   neqpar    = filehead.neqpar
   nw        = filehead.nw
   eqpar     = filehead.eqpar
   variables = filehead.variables

   mu0SI = 4π*1e-7      # H/m
   cSI   = 2.9978e8       # speed of light, [m/s]
   mpSI  = 1.6726e-27     # kg
   eSI   = 1.602e-19      # elementary charge, [C]
   AuSI  = 149597870700   # m
   RsSI  = 6.957e8        # m
   Mi    = 1.0            # Ion mass, [amu]
   Me    = 1.0/1836.15    # Electron mass, [amu]
   gamma = 5/3            # Adiabatic index for first fluid
   gammae= 5/3            # Adiabatic index for electrons
   kbSI  = 1.38064852e-23 # Boltzmann constant, [m2 kg s-2 K-1]
   e0SI  = 8.8542e-12     # [F/m]

   # This part is used to guess the units.
   # To be honest; I don`t understand the logic here. For example
   # nPa & m/s may appear in the same headline?
   if type !== ""
      typeunit = uppercase(type)
   elseif occursin("PIC",filehead[:headline])
      typeunit = "PIC"
   elseif occursin(" AU ", filehead[:headline])
      typeunit = "OUTERHELIO"
   elseif occursin(r"(kg/m3)|(m/s)", filehead[:headline])
      typeunit = "SI"
   elseif occursin(r"(nPa)|( nT )", filehead[:headline])
      typeunit = "PLANETARY"
   elseif occursin(r"(dyne)|( G)", filehead[:headline])
      typeunit = "SOLAR"
   else
      typeunit = "NORMALIZED"
   end

   if typeunit == "SI"
      xSI   = 1.0             # m
      tSI   = 1.0             # s
      rhoSI = 1.0             # kg/m^3
      uSI   = 1.0             # m/s
      pSI   = 1.0             # Pa
      bSI   = 1.0             # T
      jSI   = 1.0             # A/m^2
   elseif typeunit == "CGS"
      xSI   = 0.01            # cm
      tSI   = 1.0             # s
      rhoSI = 1000.0          # g/cm^3
      uSI   = 0.01            # cm/s
      pSI   = 0.1             # dyne/cm^2
      bSI   = 1.0e-4          # G
      jSI   = 10*cSI          # Fr/s/cm^2
   elseif typeunit == "PIC"
      # Normalized PIC units
      xSI   = 1.0             # cm
      tSI   = 1.0             # s
      rhoSI = 1.0             # g/cm^3
      uSI   = 1.0             # cm/s
      pSI   = 1.0             # dyne/cm^2
      bSI   = 1.0             # G
      jSI   = 1.0             # Fr/s/cm^2
      c0    = 1.0             # speed of light always 1 for iPIC3D
   elseif typeunit == "NORMALIZED"
      xSI   = 1.0             # distance unit in SI
      tSI   = 1.0             # time unit in SI
      rhoSI = 1.0             # density unit in SI
      uSI   = 1.0             # velocity unit in SI
      pSI   = 1.0             # pressure unit in SI
      bSI   = sqrt(mu0SI)     # magnetic unit in SI
      jSI   = 1/sqrt(mu0SI)   # current unit in SI
      c0    = 1.0             # speed of light (for Boris correction)
   elseif typeunit == "PLANETARY"
      xSI   = 6378000         # Earth radius [default planet]
      tSI   = 1.0             # s
      rhoSI = mpSI*1e6        # mp/cm^3
      uSI   = 1e3             # km/s
      pSI   = 1e-9            # nPa
      bSI   = 1e-9            # nT
      jSI   = 1e-6            # muA/m^2
      c0    = cSI/uSI         # speed of light in velocity units
   elseif typeunit == "OUTERHELIO"
      xSI   = AuSI            # AU
      tSI   = 1.0             # s
      rhoSI = mpSI*1e6        # mp/cm^3
      uSI   = 1e3             # km/s
      pSI   = 1e-1            # dyne/cm^2
      bSI   = 1e-9            # nT
      jSI   = 1e-6            # muA/m^2
      c0    = cSI/uSI         # speed of light in velocity units
   elseif typeunit == "SOLAR"
      xSI   = RsSI            # radius of the Sun
      tSI   = 1.0             # s
      rhoSI = 1e3             # g/cm^3
      uSI   = 1e3             # km/s
      pSI   = 1e-1            # dyne/cm^2
      bSI   = 1e-4            # G
      jSI   = 1e-6            # muA/m^2
      c0    = cSI/uSI         # speed of light in velocity units
   else
      @error "invalid typeunit=$(typeunit)"
   end

   # Overwrite values if given by eqpar
   for ieqpar = 1:neqpar
      var = variables[ndim+nw+ieqpar]
      if var == "xSI"
         xSI   = eqpar[ieqpar]
      elseif var == "tSI"
         tSI   = eqpar[ieqpar]
      elseif var == "uSI"
         uSI   = eqpar[ieqpar]
      elseif var == "rhoSI"
         rhoSI = eqpar[ieqpar]
      elseif var == "mi"
         mi    = eqpar[ieqpar]
      elseif var == "m1"
         m1    = eqpar[ieqpar]
      elseif var == "me"
         me    = eqpar[ieqpar]
      elseif var == "qi"
         qi    = eqpar[ieqpar]
      elseif var == "q1"
         q1    = eqpar[ieqpar]
      elseif var == "qe"
         qe    = eqpar[ieqpar]
      elseif var == "g"
         gamma = eqpar[ieqpar]
      elseif var == "g1"
         gamma = eqpar[ieqpar]
      elseif var == "ge"
         ge    = eqpar[ieqpar]
      elseif var == "c"
         c     = eqpar[ieqpar]
      elseif var == "clight"
         clight= eqpar[ieqpar]
      elseif var == "r"
         r     = eqpar[ieqpar]
      elseif var == "rbody"
         rbody = eqpar[ieqpar]
      end
   end

   # Overwrite distance unit if given as an argument
   if !isempty(distunit) xSI = distunit end

   # Overwrite ion & electron masses if given as an argument
   if !isempty(Mion)      Mi = Mion end
   if !isempty(Melectron) Me = Melectron end

   # Calculate convenient conversion factors
   if typeunit == "NORMALIZED"
      ti0  = 1.0/Mi            # T      = p/rho*Mi           = ti0*p/rho
      cs0  = 1.0               # cs     = sqrt(gamma*p/rho)  = sqrt(gs*p/rho)
      mu0A = 1.0               # vA     = sqrt(b/rho)        = sqrt(bb/mu0A/rho)
      mu0  = 1.0               # beta   = p/(bb/2)           = p/(bb/(2*mu0))
      uH0  = Mi                # uH     = j/rho*Mi           = uH0*j/rho
      op0  = 1.0/Mi            # omegap = sqrt(rho)/Mi       = op0*sqrt(rho)
      oc0  = 1.0/Mi            # omegac = b/Mi               = oc0*b
      rg0  = √Mi               # rg = sqrt(p/rho)/b*sqrt(Mi) = rg0*sqrt(p/rho)/b
      di0  = c0*Mi             # di = c0/sqrt(rho)*Mi        = di0/sqrt(rho)
      ld0  = Mi                # ld = sqrt(p)/(rho*c0)*Mi    = ld0*sqrt(p)/rho
   elseif typeunit == "PIC"
      ti0  = 1.0/Mi            # T      = p/rho*Mi           = ti0*p/rho
      cs0  = 1.0               # cs     = sqrt(gamma*p/rho)  = sqrt(gs*p/rho)
      mu0A = 4*pi              # vA     = sqrt(b/(4*!pi*rho))= sqrt(bb/mu0A/rho)
      mu0  = 4*pi              # beta   = p/(bb/(8*!pi))     = p/(bb/(2*mu0))
      uH0  = Mi                # uH     = j/rho*Mi           = uH0*j/rho
      op0  = √(4π)/Mi          # omegap = sqrt(4*!pi*rho)/Mi = op0*sqrt(rho)
      oc0  = 1.0/Mi            # omegac = b/Mi               = oc0*b
      rg0  = √Mi               # rg = sqrt(p/rho)/b*sqrt(Mi) = rg0*sqrt(p/rho)/b
      di0  = 1.0/√(4π)         # di = 1/sqrt(4*!pi*rho)*Mi   = di0/sqrt(rho)
      ld0  = 1.0/√(4π)         # ld = sqrt(p/(4*!pi))/rho*Mi = ld0*sqrt(p)/rho
   else
      qom  = eSI/(Mi*mpSI); moq = 1/qom
      ti0  = mpSI/kbSI*pSI/rhoSI*Mi       # T[K]=p/(nk) = ti0*p/rho
      cs0  = pSI/rhoSI/uSI^2              # cs          = sqrt(gs*p/rho)
      mu0A = uSI^2*mu0SI*rhoSI*bSI^(-2)   # vA          = sqrt(bb/(mu0A*rho))
      mu0  = mu0SI*pSI*bSI^(-2)           # beta        = p/(bb/(2*mu0))
      uH0  = moq*jSI/rhoSI/uSI            # uH=j/(ne)   = uH0*j/rho
      op0  = qom*√(rhoSI/e0SI)*tSI        # omegap      = op0*sqrt(rho)
      oc0  = qom*bSI*tSI                  # omegac      = oc0*b
      rg0  = moq*√(pSI/rhoSI)/bSI/xSI/√(Mi)    # rg     = rg0*sqrt(p/rho)/b
      di0  = cSI/(op0/tSI)/xSI                 # di=c/omegap = di0/sqrt(rho)
      ld0  = moq*√(pSI)/rhoSI/xSI              # ld          = ld0*sqrt(p)/rho
   end

end

"""
	showhead(file, filehead)

Displaying file header information.
"""
function showhead(file::FileList, head::NamedTuple)
   @info "filename  = $(file.name)"
   @info "filetype  = $(file.type)"
   @info "headline  = $(head.headline)"
   @info "it        = $(head.it)"
   @info "time      = $(head.time)"
   @info "gencoord  = $(head.gencoord)"
   @info "ndim      = $(head.ndim)"
   @info "neqpar    = $(head.neqpar)"
   @info "nw        = $(head.nw)"
   @info "nx        = $(head.nx)"

   if head[:neqpar] > 0
      @info "parameters = $(head.eqpar)"
      @info "coord names= $(head.variables[1:head.ndim])"
      @info "var   names= $(head.variables[head.ndim+1:head.ndim+head.nw])"
      @info "param names= $(head.variables[head.ndim+head.nw+1:end])"
   end
end

"""
	convertVTK(head, data, connectivity, filename)

Convert 3D unstructured Tecplot data to VTK. Note that if using voxel type data
in VTK, the connectivity sequence is different from Tecplot.
"""
function convertVTK(head, data, connectivity, filename="3DBATSRUS")

   nVar = length(head.variables)

   points = @view data[1:head.ndim,:]
   cells = Vector{MeshCell{Array{Int32,1}}}(undef,head.nCell)
   if head.ndim == 3
      # PLT to VTK index_ = [1 2 4 3 5 6 8 7]
      for i = 1:2
         connectivity = swaprows!(connectivity, 4*i-1, 4*i)
      end
      @inbounds for i = 1:head.nCell
         cells[i] = MeshCell(VTKCellTypes.VTK_VOXEL, connectivity[:,i])
      end
   elseif head[:ndim] == 2
      @inbounds for i = 1:head.nCell
         cells[i] = MeshCell(VTKCellTypes.VTK_PIXEL, connectivity[:,i])
      end
   end

   vtkfile = vtk_grid(filename, points, cells)

   for ivar = head.ndim+1:nVar
      if occursin("_x",head.variables[ivar]) # vector
         var1 = @view data[ivar,:]
         var2 = @view data[ivar+1,:]
         var3 = @view data[ivar+2,:]
         namevar = replace(head.variables[ivar], "_x"=>"")
         vtk_point_data(vtkfile, (var1, var2, var3), namevar)
      elseif occursin(r"(_y|_z)",head.variables[ivar])
         continue
      else
         var = @view data[ivar,:]
         vtk_point_data(vtkfile, var, head.variables[ivar])
      end
   end

   outfiles = vtk_save(vtkfile)
end

function swaprows!(X, i, j)
   m, n = size(X)
   if (1 ≤ i ≤ n) && (1 ≤ j ≤ n)
      for k = 1:n
        @inbounds X[i,k],X[j,k] = X[j,k],X[i,k]
      end
      return X
   else
      throw(BoundsError())
   end
end
