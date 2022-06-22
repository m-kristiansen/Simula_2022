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


function box_embed(x_values, y_values; padding=0.2, align::Bool=false, view::Bool=true)
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

    bp1 = gmsh.model.geo.addPoint(box[1][1], box[1][2], 0)
    bp2 = gmsh.model.geo.addPoint(box[2][1], box[2][2], 0)
    bp3 = gmsh.model.geo.addPoint(box[3][1], box[3][2], 0)
    bp4 = gmsh.model.geo.addPoint(box[4][1], box[4][2], 0)

    bl12 = gmsh.model.geo.addLine(bp1, bp2)
    bl23 = gmsh.model.geo.addLine(bp2, bp3)
    bl34 = gmsh.model.geo.addLine(bp3, bp4)
    bl41 = gmsh.model.geo.addLine(bp4, bp1)

    boxloop = gmsh.model.geo.addCurveLoop([bl12, bl23, bl34, bl41])

    boxsurface = gmsh.model.geo.addPlaneSurface([boxloop])

    gmsh.model.addPhysicalGroup(2, [boxsurface], 2) #(dim, tag, label)
    gmsh.model.setPhysicalName(2, 2, "My surface") #(dim, tag, label)

    gmsh.model.addPhysicalGroup(1, [bl12, bl23, bl34, bl41], 3) #(dim, tag, label)
    gmsh.model.setPhysicalName(1, 3, "sides") #(dim, tag, label)


    gmsh.model.geo.synchronize() # Sync CAD representation
    gmsh.model.mesh.generate(2)

    if view
        gmsh.fltk.run()
    end

    gmsh.finalize()
end


gmsh.initialize()

gmsh.model.add("Graph_1")

X = rand(4, 1)
Y = rand(4, 1)

#Define nodes
p1 = gmsh.model.geo.addPoint(X[1], Y[1], 0) #point 1, position
p2 = gmsh.model.geo.addPoint(X[2], Y[2], 0) #point 2, position
p3 = gmsh.model.geo.addPoint(X[3], Y[3], 0) #point 3, position
p4 = gmsh.model.geo.addPoint(X[4], Y[4], 0) #point 4, position

#Define lines connecting nodes
l12 = gmsh.model.geo.addLine(p1, p2) #line 1, from p1 to p4
l23 = gmsh.model.geo.addLine(p2, p3) #line 2, from p2 to p4
l34 = gmsh.model.geo.addLine(p3, p4) #line 3, from p4 to p3

gmsh.model.geo.addPhysicalGroup(1, [l12, l23, l34], 1)
#gmsh.model.geo.setPhysicalName(1, 1, "Graph")


box_embed(X, Y, padding=0.01, align = true, view = true)



"""
if distance from node_i < padding
    node_i_x = max_node_x
"""
