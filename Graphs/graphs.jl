using GridapGmsh: gmsh
using LazySets
using Random

function min_area_rect(points, hull)
    """ Get minimum area rectangle by rotation w.r.t the longest line of the convex hull.
    returns: Rotated input points and rotation matrix. """
    n = length(hull)
    # calculate edges
    edges = Vector{Float64}[] #empty array for storing edges [[length_x, lenth_y], ..]
    idx = [] #empty array to store edge indexes
    for i in 1:n
        for j in 1:n
            if (n-(i-1))>j
                append!(edges, [hull[n-(i-1)]-hull[j]])
                push!(idx, [n-(i-1), j]) #push allows to append vector [x,y]
            end
        end
    end
    edge_lengths = [sqrt(edges[i][1]^2+edges[i][2]^2) for i in 1:length(edges)]
    longest_edge_idx = findmax(edge_lengths)[2]


    x = edges[longest_edge_idx][1]
    y = edges[longest_edge_idx][2]

    angle = -atan(y, x) #Minus because we want clockwise rotation
    rot_matrix = [cos(angle) -sin(angle); sin(angle) cos(angle)] #Rotational matrix
    p = [[points[i][1] for i in 1:length(points)] [points[i][2] for i in 1:length(points)]] #reshape
    new_points = rot_matrix*p'

return new_points', rot_matrix
end


function get_bounding_box(x_values, y_values; padding=0.2, align::Bool=false, view::Bool=true)
""" Takes the x and y values corresponding to the nodes of the graph and
    Generates a bounding box around the graph and a mesh on the surface of the box
    padding prevents the graph from intersecting the bounding box."""

    if align
        #Rotate coordinates
        p = [[x_values[i], y_values[i]] for i in 1:length(x_values)]
        hull = convex_hull(p) #LazySets
        p_new, rot_matrix = min_area_rect(p, hull)
        x_values = p_new[:, 1]
        y_values = p_new[:, 2]
    end

    center_x = (findmax(x_values)[1]-findmin(x_values)[1])/2+findmin(x_values)[1]
    center_y = (findmax(y_values)[1]-findmin(y_values)[1])/2+findmin(y_values)[1]

    dx = findmax(center_x.-x_values)[1]+padding
    dy = findmax(center_y.-y_values)[1]+padding

    b1 = (center_x-dx, center_y-dy)
    b2 = (center_x-dx, center_y+dy)
    b3 = (center_x+dx, center_y+dy)
    b4 = (center_x+dx, center_y-dy)

    box = [b1, b2, b3, b4]
    if align
        #Rotate back
        b = [[box[i][1] for i in 1:length(box)] [box[i][2] for i in 1:length(box)]] #reshape
        b = (inv(rot_matrix)*b')
        box = [(b[1, 1], b[2, 1]), (b[1,2], b[2,2]), (b[1,3], b[2,3]), (b[1, 4], b[2,4])] #reshape
    end

    return box
end


function box_embed(name, graph_nodes, graph_lines; mesh_size=0.1, padding=0.2, align::Bool=false, view::Bool=true)
 """ Main function to compute bounding box and mesh for graph """
    gmsh.initialize(["", "-clmax", string(mesh_size)])

    gmsh.model.add("Graph")

    # Get bounding box from graph nodes
    b1, b2, b3, b4 = get_bounding_box(graph_nodes[:, 1], graph_nodes[:, 2], padding=padding, align=align, view=view)

    # add box points to gmsh
    bp1 = gmsh.model.geo.addPoint(b1[1], b1[2], 0)
    bp2 = gmsh.model.geo.addPoint(b2[1], b2[2], 0)
    bp3 = gmsh.model.geo.addPoint(b3[1], b3[2], 0)
    bp4 = gmsh.model.geo.addPoint(b4[1], b4[2], 0)

    # add box lines to gmsh
    bl12 = gmsh.model.geo.addLine(bp1, bp2) #bottom box edge
    bl23 = gmsh.model.geo.addLine(bp2, bp3) #left box edge
    bl34 = gmsh.model.geo.addLine(bp3, bp4) #top box edge
    bl41 = gmsh.model.geo.addLine(bp4, bp1) #right box edge

    boxloop = gmsh.model.geo.addCurveLoop([bl12, bl23, bl34, bl41])
    boxsurface = gmsh.model.geo.addPlaneSurface([boxloop])

    gmsh.model.geo.synchronize() # Sync CAD representation

    gmsh.model.addPhysicalGroup(2, [boxsurface], 1) #(dim, tag, label)
    gmsh.model.setPhysicalName(2, 1, "My surface") #(dim, tag, label)

    gmsh.model.addPhysicalGroup(1, [bl12], 2) #(dim, tag, label)
    gmsh.model.addPhysicalGroup(1, [bl23], 3)
    gmsh.model.addPhysicalGroup(1, [bl34], 4)
    gmsh.model.addPhysicalGroup(1, [bl41], 5)
    gmsh.model.setPhysicalName(1, 2, "Bottom")
    gmsh.model.setPhysicalName(1, 3, "Left")
    gmsh.model.setPhysicalName(1, 4, "Top")
    gmsh.model.setPhysicalName(1, 5, "Right")

    gmsh.model.geo.synchronize() # Sync CAD representation

    # At this point we are done with the tissue, now we need vasculature
    # So add the graph

    # add graph points to gmsh
    graph_points = [gmsh.model.geo.addPoint(x..., 0) for x in eachrow(graph_nodes)]

    # add graph lines to gmsh
    graph_lines = [gmsh.model.geo.addLine(graph_points[i], graph_points[j])
                   for (i, j) in eachrow(graph_lines)]

    gmsh.model.geo.synchronize() # Sync CAD representation

    gmsh.model.geo.addPhysicalGroup(1, graph_lines, 6)
    gmsh.model.setPhysicalName(1, 6, "Graph")
    gmsh.model.geo.addPhysicalGroup(0, [5], 7) #first graph node is tagged after box corners i.e. 5
    gmsh.model.setPhysicalName(0, 7, "starting_point")

    gmsh.model.geo.synchronize() # Sync CAD representation

    #embed the graph in the mesh
    gmsh.model.mesh.embed(1, graph_lines, 2, boxsurface) # (dim object_1, object_1, dim object_2, object_2)
    gmsh.model.mesh.generate(2)

    if view
        gmsh.fltk.run()
    end

    gmsh.write(name)
    gmsh.finalize()

    return name, b4, b2
end

#tester
if PROGRAM_FILE == basename(@__FILE__)
    include("RRT.jl")
    points  = RRT(30, 0.7, [0.0, 0.0], 1)
    (graph_nodes, graph_edges) = connect_RRT(points, 1.0)
    filename = "graph_1.msh"

    box_embed(filename, graph_nodes, graph_edges, mesh_size=1, padding=0.2, align = true, view = true)
end

"""
# Stress test
include("random_graph.jl")

(graph_nodes, graph_edges) = random_graph(10)
println(graph_edges)

box_embed(graph_nodes, graph_edges, padding=0.01, view = true)
"""
