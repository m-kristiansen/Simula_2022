using PyCall, Plots

"""
Graph based on Delaunay triangulation of a point cloud (0, 1)^2
"""
function random_graph(nvtx::Integer)
    X = rand(nvtx, 2)
    return random_graph(X)
end


"""
Graph based on Delaunay triangulation of a point cloud
"""
function random_graph(X)
    spatial = pyimport("scipy.spatial")
    tri = spatial.Delaunay(X)
    simplices = py"$(tri).simplices"
    # Shift by 1 because ...
    simplices = simplices .+ 1

    edge2cells = Dict{NTuple{2, Integer}, Vector{Integer}}()
    # The idea now is to get the edges of the triangulation and kick
    # out boundary one
    for (cid, cell) in enumerate(eachrow(simplices))
    	 cell_edges = [(cell[1], cell[2]), (cell[2], cell[3]), (cell[3], cell[1])]
	   
        for edge in cell_edges
	        # In the mapping we sort them
	        if edge[1] < edge[2]
	            key = edge
            else
	            key = (edge[2], edge[1])
            end
            
            if key ∈ keys(edge2cells)  
                push!(edge2cells[key], cid)
            else
                edge2cells[key] = [cid]
            end
         end
     end

     graph_edges = collect(filter(key -> length(edge2cells[key]) > 1, keys(edge2cells)))
     # Let's have this as matrix
     graph_edges = hcat(collect.(graph_edges)...)'
     # It may happen that in clipping boundary some vertex is left isolated
     used = collect(reduce(union, map(Set, graph_edges)))
     
     length(used) == size(X, 1) && return (X, graph_edges)

     # Otherwise redefine edges ...
     graph_nodes = X[used, :]
     # ... in the new numbering
     old2new = Dict(old => new for (new, old) in enumerate(used))

     for (row, edge) in enumerate(eachrow(graph_edges))
        graph_edges[row, 1] = old2new[edge[1]]
        graph_edges[row, 2] = old2new[edge[2]]
     end

     return (graph_nodes, graph_edges)
end


false && begin
    (graph_nodes, graph_edges) = random_graph(50)

    p = plot(graph_nodes[:, 1], graph_nodes[:, 2], seriestype=:scatter, legend=false)
    for edge in eachrow(graph_edges)
        local x = [graph_nodes[edge[1], 1], graph_nodes[edge[2], 1]]
        local y = [graph_nodes[edge[1], 2], graph_nodes[edge[2], 2]]
        plot!(x, y, legend=false)
    end
    display(p)
end