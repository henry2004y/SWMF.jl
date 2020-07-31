# Script for the benchmark between Julia, Python, IDL, and Matlab.
#
# Hongyang Zhou, hyzhou@umich.edu 03/28/2020

using PyPlot

tJulia = [163.36*1e-6, 180.8*1e-3, 1.32]
tPython= [4390.95*1e-6, 179.5*1e-3, 1.34]
tIDL   = [1970.29*1e-6, 453.5*1e-3, 6.08]
tMatlab= [19273.25*1e-6, 698.4*1e-3, 10.60]

tTable = ones(4,3)
tTable[2,:] = tPython ./ tJulia
tTable[3,:] = tIDL ./ tJulia
tTable[4,:] = tMatlab ./ tJulia

testset = 1:3

figure()
plt.rc("font", family = "serif", size = 14)
scatter(testset, tTable[1,:], s=200, marker="*", label="Julia")
scatter(testset, tTable[2,:], s=100, alpha=0.5, label="Python+Numpy")
scatter(testset, tTable[3,:], s=100, label="IDL")
scatter(testset, tTable[4,:], s=100, label="Matlab")
legend()
my_xticks = ["65KB", "317MB", "2.4GB"]
plt.xticks(testset, my_xticks)
yscale("log")
grid(b=true, which="major", color="k", linestyle="-", alpha=0.5)
grid(b=true, which="minor", color="k", linestyle="--", linewidth=0.2, alpha=0.5)

xlabel("File Size")
ylabel("Normalized Timing to Julia")

tight_layout()


## Timing for reading files in Julia
#=
using SWMF, BenchmarkTools
#filename = "3d_var_region0_0_t00001527_n00005528.out"
#filename = "z=0_raw_1_t25.60000_n00000258.out"
filename = "3d_var_region0_0_t00001205_n00037679.out"
dir = "../BATSRUS/VisAnaJulia"
#dir = "test/"
@btime data = readdata(filename, dir=dir);
=#
