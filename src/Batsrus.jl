module Batsrus
# BATSRUS data reader and analyzer.
#
# Hongyang Zhou, hyzhou@umich.edu

using Printf, WriteVTK

export Data, Vars

"""
	FileList

Type for the file information. Contains filename `name`, file type `type` of
ascii or binary, file size `bytes`, and number of snapshots `npictinfiles`.
"""
struct FileList
   name::String
   type::String
   bytes::Int64
   npictinfiles::Int64
end

"""
	Data

Primary data storage type, with fields `head` of header info, grid `x`, value
`w`, and file info `list`.
"""
struct Data{T<:AbstractFloat}
   head::NamedTuple
   x::Array{T}
   w::Array{T}
   list::FileList
end

struct Vars
   data::Dict{String, Array{Float32}}
end

include("io.jl")
include("select.jl")

end
