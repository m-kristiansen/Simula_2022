using GridapGmsh: gmsh
using LazySets

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

    angle = atan(x, y)
    rot_matrix = [[cos(angle),sin(angle)] [-sin(angle), cos(angle)]] #Rotational matrix
    p = [[points[i][1] for i in 1:length(points)] [points[i][2] for i in 1:length(points)]] #reshape
    new_points = rot_matrix*p'


return new_points', rot_matrix
end


function get_bounding_box(x_values, y_values; padding=0.2, align::Bool=false, view::Bool=true)
""" Takes the x and y values corresponding to the nodes of the graph and
    Generates a bounding box around the graph and a mesh on the surface of the box
    padding prevents the graph from intersecting the bounding box."""

    if align
        p = [[x_values[i], y_values[i]] for i in 1:length(x_values)]
        hull = convex_hull(p)
        print(hull)
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
        b = [[box[i][1] for i in 1:length(box)] [box[i][2] for i in 1:length(box)]] #reshape
        b = (inv(rot_matrix)*b')
        box = [(b[1, 1], b[2, 1]), (b[1,2], b[2,2]), (b[1,3], b[2,3]), (b[1, 4], b[2,4])] #reshape
    end
    
    return box
end


function box_embed(graph_nodes, graph_edges; padding=0.2, align::Bool=false, view::Bool=true)

    gmsh.initialize()

    gmsh.model.add("Graph_1")

    # Get the box first
    b1, b2, b3, b4 = get_bounding_box(graph_nodes[:, 1], graph_nodes[:, 2], padding=padding, align=align, view=view)

    # Now start defining model
    bp1 = gmsh.model.geo.addPoint(b1[1], b1[2], 0)
    bp2 = gmsh.model.geo.addPoint(b2[1], b2[2], 0)
    bp3 = gmsh.model.geo.addPoint(b3[1], b3[2], 0)
    bp4 = gmsh.model.geo.addPoint(b4[1], b4[2], 0)

    bl12 = gmsh.model.geo.addLine(bp1, bp2)
    bl23 = gmsh.model.geo.addLine(bp2, bp3)
    bl34 = gmsh.model.geo.addLine(bp3, bp4)
    bl41 = gmsh.model.geo.addLine(bp4, bp1)

    boxloop = gmsh.model.geo.addCurveLoop([bl12, bl23, bl34, bl41])

    boxsurface = gmsh.model.geo.addPlaneSurface([boxloop])

    gmsh.model.geo.synchronize() # Sync CAD representation
    
    gmsh.model.addPhysicalGroup(2, [boxsurface], 2) #(dim, tag, label)
    gmsh.model.setPhysicalName(2, 2, "My surface") #(dim, tag, label)

    gmsh.model.addPhysicalGroup(1, [bl12, bl23, bl34, bl41], 3) #(dim, tag, label)
    gmsh.model.setPhysicalName(1, 3, "sides") #(dim, tag, label)

    gmsh.model.geo.synchronize() # Sync CAD representation
    
    # At this point we are done with the tissue, now we need vasculature
    # So add the graph

    graph_points = [gmsh.model.geo.addPoint(x..., 0) for x in eachrow(graph_nodes)]
    
    graph_lines = [gmsh.model.geo.addLine(graph_points[i], graph_points[j])
                   for (i, j) in eachrow(graph_edges)]
    print(graph_lines)
    gmsh.model.geo.synchronize() # Sync CAD representation
    
    gmsh.model.geo.addPhysicalGroup(1, graph_lines, 1)
    # gmsh.model.geo.setPhysicalName(1, 1, "Graph")

    gmsh.model.geo.synchronize() # Sync CAD representation

    gmsh.model.mesh.embed(1, graph_lines, 2, boxsurface)
    
    if view
        gmsh.fltk.run()
    end

    gmsh.model.mesh.generate(2)

    gmsh.finalize()
end


X = [0.5, 0.5, 0.25, 0.75]
Y = [0, 0.25, 0.5, 0.5]

graph_nodes = hcat(X, Y)
graph_edges = [1 2;
               2 3;
	           2 4]

include("random_graph.jl")

(graph_nodes, graph_edges) = random_graph(30)

box_embed(graph_nodes, graph_edges, padding=0.01, align = true, view = true)
	       
# #Define nodes
# https://de.mathworks.com/matlabcentral/communitycontests/contests/4/entries/5346