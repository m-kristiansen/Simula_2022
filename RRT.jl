using Distributions
using Plots
using LinearAlgebra

function RRT(num_points, connectivity, starting_point)
    N = num_points
    x0 = starting_point[1]
    y0 = starting_point[2]
    points = zeros(N,2)
    points[1, :] = [x0, y0]
    a = hcat([x0 for i in 1:N], [x0 for i in 1:N])
    for i in 2:N
        b = (connectivity*N)*rand(Uniform(0, 1), 1,2)
        w = sum((a.-b).^2, dims = 2)
        idx = map(argmin, eachcol(w))[1]
        n = a[idx, :]
        k = atan(b[2]-n[2], b[1]-n[1])
        c = [n[1]+cos(k), n[2]+sin(k)]
        a[i, :] = c
        points[i, :] = [(c[1]+n[1])/2, (c[2]+n[2])/2]
    end
return points
end


function connect_RRT(points, max_norm)
    connections = [] #store connections [(point1, point1), ..]
    N = size(points[:, 1])[1]
    for i in 2:N #skip starting point
        iter = 0
        for j in 1:N
            if norm(hcat(points[i, :]-points[j, :])) < max_norm
                iter += 1
                if iter < 2
                    push!(connections, (j, i))
                end
            end
        end
    end

return points, hcat(collect.(connections)...)'
end

"""
points  = RRT(10, 0.9, [1.0, 1.0])
nodes, lines = connect_RRT(points, 1.0)

p = scatter(nodes[:, 1], nodes[:, 2])
for (i, j) in lines
    plot!([nodes[i, 1], nodes[j, 1]], [nodes[i, 2], nodes[j, 2]], color = :blue)
end
savefig(p, "RRT.png")
"""
