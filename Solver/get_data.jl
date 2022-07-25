include("../Graphs/RRT.jl")
include("../Graphs/graphs.jl")
include("solver.jl")
include("Binary_mask.jl")
include("Interpolate.jl")

using PyCall
np = pyimport("numpy") #need histogram2d and saver

num_points = 30             #Number of rapid random tree points to generate graph
connectivity = 1            #connectivity of graph
mesh_size = 1               #Resolution of grid
padding = 0.01              #distance graph edge to bounding box edges
partition = (128, 128)
align = false               #align bounding box?
view = false                #view msh file?
filename = "graph_1.msh"    #box embed stores file which is later used in the solver


#store bad seeds
intersecting_seeds = Int[] #store bad seeds
negativepoint_seeds = Int[]

for seed in 1:1000
    points = RRT(num_points, connectivity, seed)
    graph_nodes, graph_edges, POI = connect_RRT(points)

    if all(>=(0), points) == false
        push!(negativepoint_seeds, seed)
        continue #jump seed with negative point value, (found occurence at seed 43, 30 points 1.0 connectivity)
    end

    graph_nodes, graph_edges, POI = connect_RRT(points)

    if length(POI) > 0
        push!(intersecting_seeds, seed)
        continue #jumps seed if intersecting graph
    end

    name, xminmax, yminmax, dirichlet_tags = box_embed(filename, graph_nodes, graph_edges, mesh_size = mesh_size,padding=padding,align=align,view=view)
    #Define source terms
    f(x) = 0
    f̂(x) = 0
    g(x) = 0

    #Solve
    println("inside solver")
    uh, ûh, λh = solver(name, f,f̂,g, dirichlet_tags, write = false)
    domain = (xminmax[2], xminmax[1], yminmax[1], yminmax[2]) #normalization makes this always ≈ (0-padding, 1+padding)
    println("binary mask")
    input = binary_mask(graph_nodes, graph_edges, padding, partition)
    println("interpolating")
    target = interpolate(partition, domain, uh)

    if seed == 1
        global input_images = np.expand_dims(input, axis = 0)
        global target_images = np.expand_dims(target, axis = 0)
    else
        input_images = np.append(input_images, np.expand_dims(input, axis = 0), axis = 0)
        target_images = np.append(target_images, np.expand_dims(target, axis = 0), axis = 0)
    end
    println("done with seed: ", seed)
end

np.save("../Data/input_images_1", input_images)
np.save("../Data/target_images_1", target_images)
println("Intersections at seed: ", intersecting_seeds, "  Negative point(s) at seed: ", negativepoint_seeds[:])
