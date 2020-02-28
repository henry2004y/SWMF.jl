# Calling Matlab functions in Julia.
#
# This requires that Matlab can be called from the command line, which means the
# MATLAB_PATH should be properly set.
# Hongyang Zhou, hyzhou@umich.edu 07/17/2019

using MATLAB

filename="/Users/hyzhou/SWMF/test/BATSRUS/run_test/RESULT/GM/"*
   "z=0_raw_1_t00.00000_n00000000.out";

filehead, data = mxcall(:read_data,2,filename);

data = data["file1"];

mat"plot_data($data,$filehead,'p','plotmode','contbar')"


x = -2:0.2:2
y = x'
z = x .* exp.(-x.^2 .- y.^2)
x = collect(x)
y = collect(y)
@mput x
@mput y
eval_string("z = x .* exp(-x.^2 - y.^2)")
px, py = mxcall(:gradient, 2, z)

quiver(x,y,px,py)
mat"""
contour($x,$y,$z)
hold on
quiver($x,$y,$px,$py)
hold off
"""

# Compare the Matlab plotting to the matplotlib one.
using PyPlot
x = -2:0.2:2
y = x'
z = x .* exp.(-x.^2 .- y.^2)
contour(x,y,z)
quiver(x,y,px,py)
