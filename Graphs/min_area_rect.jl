using Plots, LazySets
using Random
Random.seed!(2018)

include("random_graph.jl")

(graph_nodes, graph_edges) = random_graph(6)

points = [[graph_nodes[i, 1], graph_nodes[i, 2]] for i in 1:length(graph_nodes[:, 1])] #reshape

hull = convex_hull(points)

function minimum_bounding_rectangle(points, hull)
    n = length(hull)
    # calculate edges
    edges = Vector{Float64}[] #empty array for storing edges [[length_x, lenth_y], ..]
    idx = [] #empty array to store hull index
    for i in 1:n
        for j in 1:n
            if (n-(i-1))>j
                append!(edges, [hull[n-(i-1)]-hull[j]])
                println("---")
                println([hull[n-(i-1)]-hull[j]], [n-(i-1), j])
                push!(idx, [n-(i-1), j]) #push allows to append vector [x,y]
            end
        end
    end
    edge_lengths = [sqrt(edges[i][1]^2+edges[i][2]^2) for i in 1:length(edges)]
    longest_edge_idx = findmax(edge_lengths)[2]

    #Bounding box axis correpsponds to longest edge, need to extend the edges to enclose all points
    x = edges[longest_edge_idx][1]
    y = edges[longest_edge_idx][2]
    println("---")
    println(hull)
    #println("---")
    #println(edges[longest_edge_idx])
    #println("---")
    angle = atan(y, x) #Minus because we want clockwise rotation
    print(angle)
    rot_matrix = [cos(angle) sin(angle); -sin(angle) cos(angle)] #Rotational matrix
    p = [[points[i][1] for i in 1:length(points)] [points[i][2] for i in 1:length(points)]] #reshape
    new_points = rot_matrix*p'


return new_points', rot_matrix
end


new_points, rot_matrix = minimum_bounding_rectangle(points, hull)
npoints = [[new_points[i, 1], new_points[i, 2]] for i in 1:length(new_points[:, 1])]


#test convex_hull
p = plot([Singleton(bi) for bi in points])
plot!([Singleton(ni) for ni in npoints])
plot!(p, VPolygon(hull), alpha=0.2)
savefig(p, "Convex_hull.png")
