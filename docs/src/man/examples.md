# Examples

## IDL format output loader

- Read data
```
filename = "1d_bin.out";
data = readdata(filename);
data = readdata(filename, verbose=true);
data = readdata(filename, npict=1);
data = readdata(filename, dir=".");
```

- 3D structured spherical coordinates
```
filename = "3d_structured.out";
data = readdata(filename, verbose=false);
```

- log file
```
logfilename = "shocktube.log";
head, data = readlogdata(logfilename)
```

## Derived variables
```
v = get_vars(data[1], head[1], ["Bx", "By", "Bz"])
B = @. sqrt(v.Bx^2 + v.By^2 + v.Bz^2)
```

## Output format conversion
ASCII Tecplot file:
```
filename = "3d_ascii.dat"
head, data, connectivity  = readtecdata(filename, IsBinary=false)
convertVTK(head, data, connectivity, outname)
```

Binary Tecplot file (`DOSAVETECBINARY=TRUE`):
```
filename = "3d_bin.dat"
head, data, connectivity  = readtecdata(filename,true)
convertVTK(head, data, connectivity, outname)
```

Multiple files:
```
using Glob
filenamesIn = "3d*.dat"
dir = "."
filenames = Vector{String}(undef,0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)
tec = readtecdata.(filenames, false) # head, data, connectivity
for (i, outname) in enumerate(filenames)
   convertVTK(tec[i][1], tec[i][2], tec[i][3], outname[1:end-4])
end
```

If each individual file size is large, consider using:
```
using Glob
filenamesIn = "3d*.dat"
dir = "."
filenames = Vector{String}(undef,0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)
for (i, outname) in enumerate(filenames)
   head, data, connectivity = readtecdata(filenames, false)
   convertVTK(head, data, connectivity, outname[1:end-4])
end
```

Multiple files in parallel:
```
using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("VisAnaJulia");
@everywhere using VisAna, Glob

filenamesIn = "cut*.dat"
dir = "."
filenames = Vector{String}(undef,0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)

@sync @distributed for outname in filenames
   println("filename=$(outname)")
   head, data, connectivity = readtecdata(outname, false)
   convertVTK(head, data, connectivity, outname[1:end-4])
end
```
