module SWMF
# BATSRUS data reader and analyzer.
#
# Hongyang Zhou, hyzhou@umich.edu

using Glob, Printf, WriteVTK

export Data, FileList, Vars

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

end
