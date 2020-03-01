# Test of SWMF data loader

# To do: add more test outputs and redesign the tests!

ENV["MPLBACKEND"]="agg" # no GUI

using SWMF, Test

@testset "reading 1D ascii" begin
   filename = "1d__raw_2_t25.60000_n00000258.out"
   data = readdata(filename, verbose=true)
   @test isa(data.head, NamedTuple)
   @test extrema(data.x) == (-127.5, 127.5)
   @test extrema(data.w) == (-0.79960780498, 1.9394335293)
end

@testset "reading 2D structured binary" begin
   filename = "z=0_raw_1_t25.60000_n00000258.out"
   data = readdata(filename)
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
   data = readdata(filename)
   plotrange = [-50.0, 50.0, -0.5, 0.5]
   X, Z, p = cutdata(data, "p", cut='y', cutPlaneIndex=1, plotrange=plotrange)
   @test p[1] ≈ 0.560976f0
   @test p[2] ≈ 0.53704995f0
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
