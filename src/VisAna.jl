module VisAna
# BATSRUS data reader and analyzer.
#
# Hongyang Zhou, hyzhou@umich.edu

using Glob, Printf, PyCall, Dierckx, WriteVTK

export Data, FileList

struct Data{T}
   x::Array{T}
   w::Array{T}
end

struct FileList
   name::String
   type::String
   bytes::Int64
   npictinfiles::Int64
end

struct Vars
   data::Dict{String, Array{Float32}}
end

include("io.jl")
include("visual.jl")
include("trace.jl")
#include("visual_plot.jl")

end
