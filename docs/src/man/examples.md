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
v = get_vars(data, ["Bx", "By", "Bz"])
B = @. sqrt(v.Bx^2 + v.By^2 + v.Bz^2)
```

## Output format conversion
We can convert 2D/3D BATSRUS outputs `*.dat` to VTK formats. The default converted filename is `out.vtu`.

ASCII Tecplot file (supports both `tec` and `tcp`) and binary Tecplot file (set `DOSAVETECBINARY=TRUE` in BATSRUS `PARAM.in`):
```
filename = "x=0_mhd_1_n00000050.dat"
#filename = "3d_ascii.dat"
#filename = "3d_bin.dat"
head, data, connectivity = readtecdata(filename)
convertTECtoVTU(head, data, connectivity)
```

3D structured IDL file (`gridType=1` returns rectilinear `vtr` file, `gridType=2` returns structured `vts` file):
```
filename = "3d_structured.out"
convertIDLtoVTK(filename, gridType=1)
```

3D unstructured IDL file together with header and tree file:
```
filetag = "3d__var_1_n00002500"
convertIDLtoVTK(filetag)
```

!!! note
    The file suffix should not be provided for this to work correctly!

Multiple files:
```
using Glob
filenamesIn = "3d*.dat"
dir = "."
filenames = Vector{String}(undef,0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)
tec = readtecdata.(filenames) # head, data, connectivity
for (i, outname) in enumerate(filenames)
   convertTECtoVTU(tec[i][1], tec[i][2], tec[i][3], outname[1:end-4])
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
   head, data, connectivity = readtecdata(outname)
   convertTECtoVTU(head, data, connectivity, outname[1:end-4])
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
   head, data, connectivity = readtecdata(outname)
   convertTECtoVTU(head, data, connectivity, outname[1:end-4])
end
```
