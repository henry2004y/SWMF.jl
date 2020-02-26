# Python calls Julia
#
# You may encounter many issues with the version and dependency libraries between Julia and Python.
# An easy workaround is to use python-jl bundled with the PyJulia package.
# The arrays are automatically recognized as numpy ndarrays.

from julia import VisAna
filename = '1d__raw_2_t25.60000_n00000258.out'
head, data, filelist = VisAna.readdata(filename);

