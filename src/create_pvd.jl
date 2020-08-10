using Glob, LightXML

filenames = glob("*.vtu")

# create an empty XML document
xdoc = XMLDocument()

# create & attach a root node
xroot = create_root(xdoc, "VTKFile")

type = "Collection"
byte_order = "LittleEndian"
compressor = "vtkZLibDataCompressor"

if VERSION < v"1.4" 
   set_attributes(xroot; type=type, byte_order=byte_order, compressor=compressor)
else
   set_attributes(xroot; type, byte_order, compressor)
end

    
# create the first child
xs1 = new_child(xroot, "Collection")

for filename in filenames
   i_end  = findfirst("_n",filename)[1] - 1

   second = parse(Int32, filename[i_end-1:i_end])
   minute = parse(Int32, filename[i_end-3:i_end-2])
   timestep = 60*minute + second
   @show timestep

   e = new_child(xs1, "DataSet")

   set_attributes(e; timestep=timestep, group="", part="0", file=filename)
end

name_index = findfirst("_t",filenames[1])[1] - 1

save_file(xdoc, filenames[1][1:name_index]*".pvd")
