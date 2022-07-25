using Gridap
using Gridap.CellData

function interpolate(partition, domain, uh, write::Bool = false)
    """ Takes a singlefield FE function and interpolates on square grid, the returned value is a one channel image
    with size partition containing the values of uh. """

    𝒯₁ = CartesianDiscreteModel(domain, partition)

    for i in 1:100 # the length of this is highly ambigous
        try
            sm = KDTreeSearch(num_nearest_vertices=i) #this is default = 1 but we increase if AssertionError
            iuh = Interpolable(uh; tol=1E-12, searchmethod=sm)
            reffe_u = ReferenceFE(lagrangian, Float64, 0)
            V₁ = FESpace(𝒯₁, reffe_u; conformity=:L2)
            global fu = interpolate_everywhere(iuh,V₁)
            break #break loop if assertion error not thrown
        catch
            println("trying num_nearest_vertices= ", i+1)
            continue
        end
    end

    if write
        Gridap.writevtk(get_triangulation(fu), "interpolate", cellfields=["fu"=>fu])
    end

    imgvector = vcat(fu.cell_dof_values...) #collect in [len(cell_dof_values), 1] vector
    img = reshape(imgvector, partition) #reshape to square matrix with (x, y) orientation
return img
end


false && begin
# stress test
    include("../Graphs/RRT.jl")
    include("solver.jl")
    include("../Graphs/graphs.jl")

    partition = (128,128)
    num_points = 30             #Number of rapid random tree points to generate graph
    connectivity = 1            #connectivity of graph
    mesh_size = 1               #Resolution of grid
    padding = 0.01              #distance graph edge to bounding box edges
    align = false               #align bounding box?
    view = false                #view msh file?
    filename = "graph_1.msh"    #box embed stores file which is later used in the solver

    intersecting_seeds = Int[] #store bad seeds
    negativepoints_seeds = Int[]
    for seed in rand(1:1000, 20)
        points = RRT(num_points, connectivity, seed)
        if all(>=(0), points) == false
            push!(negativepoints_seeds, seed)
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
        uh, ûh, λh = solver(name, f,f̂,g, dirichlet_tags, write = true)
        domain = (xminmax[2], xminmax[1], yminmax[1], yminmax[2]) #normalization makes this always equal (0-padding, 1+padding)
        img = interpolate(partition, domain, uh)
        heatmap(img', show = true)
    end
    println("Intersections at seed: ", intersecting_seeds, "  Negative point(s) at seed: ", negativepoints_seeds[:])
end
