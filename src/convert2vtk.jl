#!/usr/bin/env julia --startup-file=no
# Script for converting BATSRUS unstructured binary outputs to VTK formats.
# Example usage:
# julia -p 8 convert_parallel.jl
#
# Currently the relative location of the output files to be processed is "IO2".
#
# Hongyang Zhou, hyzhou@umich.edu

using Distributed
using Pkg
if "VisAna" ∈ keys(Pkg.installed())
   @everywhere using VisAna
else
   @warn "VisAna not installed. Consider installing the package by
   `Pkg.add(PackageSpec(url=\"https://github.com/henry2004y/VisAnaJulia\",
    rev=\"master\"))`"
end
@everywhere using Glob

dir = "IO2"

# Concatenate files
headers = glob("*.T", dir)
@sync @distributed for header in headers
   name = basename(header)[1:end-2]
   println("header=$(header)")
   nameIn = glob(name*"*.tec", dir)
   nameOut= joinpath(dir, name*".dat")
   run(pipeline(`cat $(nameIn)`, nameOut))
   rm(header)
   rm.(nameIn)
end

filenamesIn = "*.dat"
filenames = Vector{String}(undef,0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)
# Do not work on files that have already been converted
# This won't work if new files are generated with exactly the same name!
filenames = [fname for fname in filenames if ~isfile(fname[1:end-3]*"vtu")]

@sync @distributed for outname in filenames
   println("filename=$(outname)")
   head, data, connectivity = readtecdata(outname, false)
   convertVTK(head, data, connectivity, outname[1:end-4])
end

# Choose whether or not to delete the original *.dat files.
#datFiles = glob("*.dat", dir)
# Parallel:
#@sync @distributed for file in datFiles
#   rm(file)
#end
# Serial (might be faster if I/O bound):
#rm.(datfiles)