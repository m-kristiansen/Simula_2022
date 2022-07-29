using Distributions
using Random
using Plots
using LinearAlgebra:norm

function RRT(num_points, connectivity, seed; starting_point=[0.0,0.0])
    Random.seed!(seed)
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


function line_intersection(p1, p2, p3, p4)
    """
    return true and point of intersection if lines intersect
    """

    x1, y1 = p1 # a point on the first line
    x2, y2 = p2 # another point on the first line
    x3, y3 = p3 # a point on the second line
    x4, y4 = p4 # another point on the second line

    #check paralell:
    if (x1-x2)*(y3-y4) == (y1-y2)*(x3-x4)
        return false
    end

    px = ((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
    py = ((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))

    #scalar value desciding point of intersection on line 1
    t = ((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)) #scalar value desciding point of intersection on line 1

    #scalar value desciding point of intersection on line 2
    u = ((x1-x3)*(y1-y2)-(y1-y3)*(x1-x2))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))

    s = [t, u]

    if all(>=(0), s) & all(<=(1), s)
        return true, [px,py]
    end

return false, []
end



function connect_RRT(points, ;max_norm::Float64 = 1.0)
    connections = [] #store connections [(point1, point1), ..]
    POI = [] #if intersect, append here
    nodes = points
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
    #check if intersect
    for i in 1:length(connections)
        for j in 1:length(connections)
            if i<j
                if (connections[i][1] ∉ connections[j]) & (connections[i][2] ∉ connections[j])
                    bool, p = line_intersection(points[connections[i][1], :], points[connections[i][2], :],
                        points[connections[j][1], :], points[connections[j][2], :])
                    if bool
                        push!(POI, p)
                        println("intersection at: ", POI)
                    end
                end
            end
        end
    end
    nodes = nodes/maximum(nodes) #normalize
return nodes, hcat(collect.(connections)...)', POI
end

false && begin
    #visual stress test
    for seed in rand(1:1000, 100)
        println("seed: ", seed)
        global points = RRT(30, 1, seed)
        nodes, lines, POI = connect_RRT(points)

        scatter(nodes[:, 1], nodes[:, 2])

        for (i,j) in eachrow(lines)
            plot!([nodes[i, 1], nodes[j, 1]],[[nodes[i, 2], nodes[j, 2]]], color = :blue, show = true, legend = false,
            xlim = (0,1.01), ylim = (0,1.01))
        end
        sleep(1.0) #time delay between plots
    end
    #Figure
    function create_figure(seeds)
        points = RRT(30, 1, seeds[1])
        nodes, lines, POI = connect_RRT(points)
        a = scatter(nodes[:, 1], nodes[:, 2])
        for (i,j) in eachrow(lines)
            plot!([nodes[i, 1], nodes[j, 1]],[[nodes[i, 2], nodes[j, 2]]], color = :blue, legend = false)
        end
        points = RRT(30, 1, seeds[2])
        nodes, lines, POI = connect_RRT(points)
        b = scatter(nodes[:, 1], nodes[:, 2])
        for (i,j) in eachrow(lines)
            plot!([nodes[i, 1], nodes[j, 1]],[[nodes[i, 2], nodes[j, 2]]], color = :blue, legend = false)
        end
        points = RRT(30, 1, seeds[3])
        nodes, lines, POI = connect_RRT(points)
        c = scatter(nodes[:, 1], nodes[:, 2])
        for (i,j) in eachrow(lines)
            plot!([nodes[i, 1], nodes[j, 1]],[[nodes[i, 2], nodes[j, 2]]], color = :blue, legend = false)
        end
        points = RRT(30, 1, seeds[4])
        nodes, lines, POI = connect_RRT(points)
        d = scatter(nodes[:, 1], nodes[:, 2])
        for (i,j) in eachrow(lines)
            plot!([nodes[i, 1], nodes[j, 1]],[[nodes[i, 2], nodes[j, 2]]], color = :blue, legend = false)
        end
        e = plot(a, b, c, d, layout=(1,4), show = true, size =(1200, 300))
        savefig(e, "4_graphs")
    end
    #create_figure([1,2,3,4])
end
