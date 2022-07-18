using Test
using Gridap
using Gridap.CellData
using Gridap.Visualization

path = "/Users/martinkristiansen/Desktop/Simula_2022"

include(path*"/Graphs/RRT.jl")
include(path*"/Graphs/graphs.jl")
include(path*"/Solver/solver.jl")

seed = 1                    #Random seed
num_points = 30             #Number of rapid random tree points to generate graph
connectivity = 0.6          #connectivity of graph
starting_point = [0.2, 0.2] #starting point for graph, include padding
max_norm = 1.0              #max distance to connect points
mesh_size = 1               #Resolution of grid
padding = 0.2               #distance graph edge to bounding box edges
align = false               #align bounding box?
view = false                #view msh file?
filename = "graph_1.msh"    #box embed stores file which is later used in the solver

#Generate random graph and embed in mesh
points = RRT(num_points, connectivity, starting_point, seed)
(graph_nodes, graph_edges) = connect_RRT(points, max_norm)
name, xminmax, yminmax = box_embed(filename, graph_nodes, graph_edges, mesh_size = mesh_size,padding=padding,align=align,view=view)

#Define source terms
f(x) = 0.01*(x[1]+x[2])
f̂(x) = 1
g(x) = exp(-(x[1]+x[2]))

#Solve
uh, ûh, λh = solver(name, f,f̂,g, write = false)


domain = (xminmax[2], xminmax[1], yminmax[1], yminmax[2])

partition = (20,20) #not bigger than xmax/mesh_size, ymax/mesh_size

# Get target data
𝒯₁ = CartesianDiscreteModel(domain, partition)
iuh = Interpolable(uh)
reffe_u = ReferenceFE(lagrangian, Float64, 0)
V₁ = FESpace(𝒯₁, reffe_u; conformity=:L2)

fu = interpolate(iuh,V₁)
#Gridap.writevtk(get_triangulation(fu), "interpolate", cellfields=["fu"=>fu])
imgvector = vcat(fu.cell_dof_values...) #collect in [len(cell_dof_values), 1] vector
img = reshape(imgvector, partition) #reshape to square matrix with (x, y) orientation

#test target data
xmin = xminmax[2]
xmax = xminmax[1]
ymin = yminmax[1]
ymax = yminmax[2]
# img = img' to get same oriantation as write vtk in heatmap (Plots) https://github.com/JuliaPlots/Makie.jl/issues/205
println("min= ", (xmin, ymin))
println("max= ", (xmax, ymax))
@test img[1,1] == fu(VectorValue(xmin, ymin))
@test img[1, end] == fu(VectorValue(xmin, ymax))
@test img[end, 1] == fu(VectorValue(xmax, ymin))
@test img[end, end] == fu(VectorValue(xmax, ymax))

#Get input data
function getEquidistantPoints(p1, p2, parts)
    new_points = zip(LinRange(p1[1], p2[1], parts+2)[2: end-1], LinRange(p1[2], p2[2], parts+2)[2: end-1])
    new_points = hcat(collect.(new_points)...)' #1x2 matrix
    return new_points
end

function increase_resolution(nodes, lines, num_points)
    x = nodes[:, 1]
    y = nodes[:, 2]
    for (i, j) in eachrow(lines)
        new_points = getEquidistantPoints(nodes[i, :], nodes[j, :], num_points)
        x = vcat(x, new_points[:, 1]) # append, placement doesnt matter
        y = vcat(y, new_points[:, 2])
    end
    return x,y
end

using PyCall
np = pyimport("numpy") #need histogram2d

partition = [30, 30]

x = graph_nodes[:, 1]
y = graph_nodes[:, 2]

M_new, xedges, yedges = np.histogram2d(x, y, partition)           #Make 2D image
M_new = @. ifelse(M_new == 0, M_new, 1)                           #Make binary


function cont_image(M_new)
    M = np.zeros(partition) #just to start
    i = 1
    while (all(M_new == M) == false)
        M = M_new
        x,y = increase_resolution(graph_nodes, graph_edges, i)
        M_new = np.histogram2d(x, y, partition)[1]
        M_new = @. ifelse(M_new == 0, M_new, 1)
        i+=1
    end
    return M_new
end

graph_image = cont_image(M_new)
