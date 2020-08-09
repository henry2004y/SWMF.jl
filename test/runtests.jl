# Test of BATSRUS data loader

ENV["MPLBACKEND"]="agg" # no GUI

using Batsrus, Test, SHA

function filecmp(path1::AbstractString, path2::AbstractString)
   stat1, stat2 = stat(path1), stat(path2)
   if !(isfile(stat1) && isfile(stat2)) || filesize(stat1) != filesize(stat2)
      return false # or should it throw if a file doesn't exist?
   end
   stat1 == stat2 && return true # same file
   open(path1, "r") do file1
      open(path2, "r") do file2
         buf1 = Vector{UInt8}(undef, 32768)
         buf2 = similar(buf1)
         while !eof(file1) && !eof(file2)
            n1 = readbytes!(file1, buf1)
            n2 = readbytes!(file2, buf2)
            n1 != n2 && return false
            0 != Base._memcmp(buf1, buf2, n1) && return false
         end
         return eof(file1) == eof(file2)
      end
   end
end

@testset "reading 1D ascii" begin
   filename = "1d__raw_2_t25.60000_n00000258.out"
   data = readdata(filename, dir="data", verbose=true)
   @test isa(data.head, NamedTuple)
   @test extrema(data.x) == (-127.5, 127.5)
   @test extrema(data.w) == (-0.79960780498, 1.9394335293)
end

@testset "reading 2D structured binary" begin
   filename = "z=0_raw_1_t25.60000_n00000258.out"
   data = readdata(filename, dir="data")
   @test data.head.time == 25.6f0
   @test extrema(data.x) == (-127.5f0, 127.5f0)
   @test extrema(data.w) == (-0.79985905f0, 1.9399388f0)
end

@testset "reading 2D unstructured binary" begin
   #filename = "z=0_raw_1_t25.60000_n00000258.out"
   #data = readdata(filename)
end

@testset "reading 3D structured binary" begin
   filename = "3d_raw.out"
   data = readdata(filename, dir="data")
   plotrange = [-50.0, 50.0, -0.5, 0.5]
   X, Z, p = cutdata(data, "p", cut='y', cutPlaneIndex=1, plotrange=plotrange)
   @test p[1] ≈ 0.560976f0
   @test p[2] ≈ 0.53704995f0
   vars = get_vars(data, ["p"]) 
   @test size(vars.p) == (8,8,8) 
end

@testset "log" begin
   logfilename = "data/log_n000001.log"
   head, data = readlogdata(logfilename)
   @test isa(head, NamedTuple)
   @test extrema(data) == (-0.105, 258.0)
end

@testset "vtk" begin
   @info("VTK conversion test.")
   filename = "data/3d_bin.dat"
   head, data, connectivity = readtecdata(filename)
   @test maximum(connectivity) ≤ head[:nNode] # check if it's read correctly
   convertVTK(head, data, connectivity)
   sha_str = bytes2hex(open(sha1, "out.vtu"))
   @test sha_str == "5b04747666542d802357dec183177f757754a254"
   rm("out.vtu")
end
