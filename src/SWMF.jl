module SWMF
# BATSRUS data reader and analyzer.
#
# Hongyang Zhou, hyzhou@umich.edu

using Printf, WriteVTK

export Data, Vars

struct FileList
   name::String
   type::String
   bytes::Int64
   npictinfiles::Int64
end

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
