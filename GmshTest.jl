# This file reimplements gmsh/tutorial/t1.geo in Julia.
# See the corresponding Python tutorial for detailed comments.
using GridapGmsh: gmsh

gmsh.initialize()

gmsh.model.add("t1")

lc = 1e-2
gmsh.model.geo.addPoint(0, 0, 0, lc, 1) #point 1, position (0,0,0)
gmsh.model.geo.addPoint(.1, 0,  0, lc, 2) #point 2, position (0.1, 0, 0)
gmsh.model.geo.addPoint(.1, .3, 0, lc, 3) #point 3, position (0.1, 0.3, 0)

p4 = gmsh.model.geo.addPoint(0, .3, 0, lc) #point 4, position (0, 0.3, 0)

gmsh.model.geo.addLine(1, 2, 1) #line 1, from p1 to p2
gmsh.model.geo.addLine(3, 2, 2) #line 2, from p2 to p3
gmsh.model.geo.addLine(3, p4, 3) #line 3, from p3 to p4
gmsh.model.geo.addLine(4, 1, p4) #line 4, from p4 to p1

gmsh.model.geo.addCurveLoop([4, 1, -2, 3], 1) #Define closed loop, negative values defines reversed orientation
gmsh.model.geo.addPlaneSurface([1], 1) #Define surface inside loop 1

gmsh.model.geo.synchronize() # Sync CAD representation

gmsh.model.addPhysicalGroup(0, [1, 2], 1) # Physical point (dimension 0, point 1 and 2)
gmsh.model.addPhysicalGroup(1, [1, 2], 2) # Physical curve (dim 1, point 1, 2)
gmsh.model.addPhysicalGroup(2, [1], 6)

gmsh.model.setPhysicalName(2, 6, "My surface")

gmsh.model.mesh.generate(2)

path = "/Users/martinkristiansen/Desktop/Simula_2022/"
if @isdefined path
	gmsh.write(joinpath(path, "t1.msh"))
end

if "-gui" in ARGS
    gmsh.fltk.run()
end

gmsh.finalize()
