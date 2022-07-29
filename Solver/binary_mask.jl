using PyCall
np = pyimport("numpy") #need histogram2d

#Get input data
function getEquidistantPoints(p1, p2, parts)
    "Get points between points, parts is the number of new points"
    new_points = zip(LinRange(p1[1], p2[1], parts+2)[2: end-1], LinRange(p1[2], p2[2], parts+2)[2: end-1])
    new_points = hcat(collect.(new_points)...)' #1x2 matrix
    return new_points
end

function increase_resolution(nodes, lines, num_points, padding; static = true)
    "Append EquidistantPoints to nodes, since the histogram is computed for the graph_nodes, with top
    right corner value equal max node, keep_range adds fictive node s.t. the binary grid is static and the
    graph moves within it."
    x = nodes[:, 1]
    y = nodes[:, 2]
    for (i, j) in eachrow(lines)
        new_points = getEquidistantPoints(nodes[i, :], nodes[j, :], num_points)
        x = vcat(x, new_points[:, 1]) # append, placement doesnt matter
        y = vcat(y, new_points[:, 2])
    end
    if static
        push!(x, 1.0+padding)
        push!(y, 1.0+padding)
        push!(x, 0.0-padding)
        push!(y, 0.0-padding)
    end
    return x, y
end

function binary_mask(graph_nodes, graph_edges, padding, partition)
""" Main method, uses the above methods to create binary mask. The EquidistantPoints are increased until no change is
made to the mask, i.e. the binary mask is continous."""
    #first iteration
    x = graph_nodes[:, 1]
    y = graph_nodes[:, 2]
    push!(x, 1.0+padding) #ensures correct size
    push!(y, 1.0+padding)
    push!(x, 0.0-padding) #ensures correct size
    push!(y, 0.0-padding)


    M_new, xedges, yedges = np.histogram2d(x, y, np.array(partition)) #Make 2D image
    M_new[end, end] = 0 #remove fictive node
    M_new[1,1] = 0
    M_new = @. ifelse(M_new == 0, M_new, 1) #Make binary
    M = zeros(partition[1], partition[2]) #just to start

    i = 1
    while (all(M_new == M) == false)
        M_prev = M
        M = M_new
        x,y = increase_resolution(graph_nodes, graph_edges, i, padding)
        M_new = np.histogram2d(x, y, partition)[1]
        M_new[end, end] = 0 #remove fictive node
        M_new[1,1] = 0
        M_new = @. ifelse(M_new == 0, M_new, 1)
        #heatmap(M_prev, show = true)
        i+=1
        if i == 300
            # We break after 300 new points between each node pair. Reason: seed 3504
            break
        end
    end
    println("num_iter: ",i)
    return M_new
end

false && begin
# Visual stress test
    using Plots
    using Random
    include("../Graphs/RRT.jl")
    for seed in rand(1:1000, 100)
        points = RRT(30, 1, seed)
        graph_nodes, graph_edges = connect_RRT(points)
        img = binary_mask(graph_nodes, graph_edges, 0.01, (128,128))
        heatmap(img', show = true) #heatmap transposes image
    end
    function create_figure(seeds)
        points = RRT(30, 1, seeds[1])
        graph_nodes, graph_edges = connect_RRT(points)
        img = binary_mask(graph_nodes, graph_edges, 0.01, (128,128))
        a = heatmap(img', colorbar = false) #heatmap transposes image
        points = RRT(30, 1, seeds[2])
        graph_nodes, graph_edges = connect_RRT(points)
        img = binary_mask(graph_nodes, graph_edges, 0.01, (128,128))
        b = heatmap(img', colorbar = false) #heatmap transposes image
        points = RRT(30, 1, seeds[3])
        graph_nodes, graph_edges = connect_RRT(points)
        img = binary_mask(graph_nodes, graph_edges, 0.01, (128,128))
        c = heatmap(img', colorbar = false) #heatmap transposes image
        points = RRT(30, 1, seeds[4])
        graph_nodes, graph_edges = connect_RRT(points)
        img = binary_mask(graph_nodes, graph_edges, 0.01, (128,128))
        d = heatmap(img', colorbar = false) #heatmap transposes image
        e = plot(a, b, c, d, layout=(1,4), show = true, size =(1200, 310))
        savefig(e, "4_binary_masks")
    end
    #create_figure([1,2,3,4])
    
end
