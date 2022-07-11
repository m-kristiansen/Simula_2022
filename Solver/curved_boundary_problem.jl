using GridapGmsh: gmsh

gmsh.initialize(["", "-clmax", string(0.5)])
gmsh.model.add("Curved_boundary")

# add box points to gmsh
bp1 = gmsh.model.geo.addPoint(0, 0, 0)
bp2 = gmsh.model.geo.addPoint(1, 0, 0)
bp3 = gmsh.model.geo.addPoint(1, 1, 0)
bp4 = gmsh.model.geo.addPoint(0, 1, 0)

# add box lines to gmsh
bl23 = gmsh.model.geo.addLine(bp2, bp3) #left box edge
bl34 = gmsh.model.geo.addLine(bp3, bp4) #top box edge
bl41 = gmsh.model.geo.addLine(bp4, bp1) #right box edge

# add curved bottom coundary
center = gmsh.model.geo.addPoint(0.5, -1, 0)
bl12 = gmsh.model.geo.addCircleArc(bp1, center, bp2) #bottom box curve

boxloop = gmsh.model.geo.addCurveLoop([bl12, bl23, bl34, bl41])
boxsurface = gmsh.model.geo.addPlaneSurface([boxloop])

gmsh.model.geo.synchronize() # Sync CAD representation

gmsh.model.addPhysicalGroup(2, [boxsurface], 1) #(dim, tag, label)
gmsh.model.setPhysicalName(2, 1, "my surface") #(dim, tag, label)

gmsh.model.addPhysicalGroup(1, [bl12], 2) #(dim, tag, label)
gmsh.model.addPhysicalGroup(1, [bl23], 3)
gmsh.model.addPhysicalGroup(1, [bl34], 4)
gmsh.model.addPhysicalGroup(1, [bl41], 5)
gmsh.model.setPhysicalName(1, 2, "bottom")
gmsh.model.setPhysicalName(1, 3, "left")
gmsh.model.setPhysicalName(1, 4, "top")
gmsh.model.setPhysicalName(1, 5, "right")

gmsh.model.mesh.generate(2)

gmsh.fltk.run()
gmsh.write("curved.msh")
gmsh.finalize()
