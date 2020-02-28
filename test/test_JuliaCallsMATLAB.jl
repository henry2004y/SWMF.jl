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


##
using PyCall, PyPlot

filename = "y*.outs"
filehead, data, filelist = readdata(filename,npict=20,verbose=false);

X = vec(data[1].x[:,:,1])
Y = vec(data[1].x[:,:,2])
W = vec(data[1].w[:,:,10])

# Perform linear interpolation of the data (x,y) on grid(xi,yi)
plotrange = zeros(4)
plotrange[1] = -3 # minimum(X)
plotrange[2] = 3  # maximum(X)
plotrange[3] = -3 # minimum(Y)
plotrange[4] = 3  # maximum(Y)

#plotinterval = X[2] - X[1]
plotinterval = 0.01

xi = range(plotrange[1], stop=plotrange[2], step=plotinterval)
yi = range(plotrange[3], stop=plotrange[4], step=plotinterval)

triang = matplotlib.tri.Triangulation(X,Y)
interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
np = pyimport("numpy")
Xi, Yi = np.meshgrid(xi, yi)
wi = interpolator(Xi, Yi)

# mask a circle in the middle:
interior = (Xi.^2 .+ Yi.^2 .< 1.0);
wi[interior] .= np.ma.masked

c = contourf(xi,yi,wi,50)
colorbar()
using MATLAB
px, py = mxcall(:gradient, 2, wi)
