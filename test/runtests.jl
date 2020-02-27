# Test of SWMF data loader

# To do: add more test outputs and redesign the tests!

ENV["MPLBACKEND"]="agg" # no GUI

using SWMF, Test

@testset "reading 1D ascii" begin
   filename = "1d__raw_2_t25.60000_n00000258.out"
   head, data, list = readdata(filename, verbose=true)
   @test isa(head[1], Dict)
   @test isa(data[1], Data)
   @test isa(list[1], FileList)
end

@testset "reading 2D structured binary" begin
   filename = "z=0_raw_1_t25.60000_n00000258.out"
   head, data, list = readdata(filename)
   @test isa(head[1], Dict)
   @test isa(data[1], Data)
   @test isa(list[1], FileList)
end

@testset "reading 2D unstructured binary" begin
   #filename = "z=0_raw_1_t25.60000_n00000258.out"
   #head, data, list = readdata(filename)
end

@testset "reading 3D structured binary" begin
   #filename = "3d_structured.out";
   #head, data, list = readdata(filename,verbose=false);
end

@testset "log" begin
   logfilename = "log_n000001.log"
   head, data = readlogdata(logfilename)
   @test isa(head, Dict)
   @test isa(data, Array)
end

@testset "vtk" begin
   @info("VTK conversion test.")
   filename = "3d_bin.dat"
   head, data, connectivity  = readtecdata(filename, IsBinary=true)
   @test maximum(connectivity) ≤ head[:nNode] # check if it's read correctly
   convertVTK(head, data, connectivity)
   @test isfile("3DBATSRUS.vtu")
   rm("3DBATSRUS.vtu")
end
