using Gridap
using GridapGmsh
include("graph_distance.jl")
"""
Solve the coupled 2D/1D problem using the Gridap functionality on gmsh generated file
-Δu = f                on Ω
-Δû - λ = f̂            on Γ       λ = -∇u⋅ν
u-û = g                on Γ

"""
function solver(filename,f, f̂, g, dirichlet_tags; κ=1, κ̂=1, write::Bool=true)
    model = GmshDiscreteModel(filename)
    #model = simplexify(model)

    order = 2 #Choose order of local interpolation, order 0 => piecewise linear, order 1 => hat function etc..

    #Define local basis functions
    reffe_u = ReferenceFE(lagrangian, Float64, order)
    reffe_û = ReferenceFE(lagrangian, Float64, order)
    reffe_μ = ReferenceFE(lagrangian,Float64, order-2)

    #where should our solution live
    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model, tags = ["Graph"])
    ∂Ω = BoundaryTriangulation(model, tags = ["Bottom", "Left", "Right", "Top"])

    V = TestFESpace(Ω,reffe_u, conformity=:H1)
    V̂ = TestFESpace(Γ,reffe_û,dirichlet_tags = dirichlet_tags, conformity=:H1)
    M = TestFESpace(Γ,reffe_μ, conformity=:L2) #create test space for lagrange multiplier

    Y = MultiFieldFESpace([V,V̂,M])

    dc(x) = 1/(1+x) #takes distance from root to tip as input and returns the dirchlet value to be used

    distances = zeros(length(dirichlet_tags))
    dirichlet_cond = [1.0]
    for i in 2:length(dirichlet_tags)
        distances[i] = graph_distance(model, "Graph", dirichlet_tags[1], dirichlet_tags[i])
        push!(dirichlet_cond, dc(distances[i]))
    end

    # Define global trial space
    U = TrialFESpace(V) #include u to enforce dirichlet tradtionally
    Û = TrialFESpace(V̂, dirichlet_cond) # 1 dirichlet on start of graph
    Λ = TrialFESpace(M) #trial space for lagrange multiplier
    X = MultiFieldFESpace([U,Û,Λ])

    #integration measure
    degree = order
    dΩ = Measure(Ω,degree)
    dΓ = Measure(Γ,degree)
    d∂Ω = Measure(∂Ω, degree)
    ν = get_normal_vector(∂Ω)

    h(x) = VectorValue(0,0) #Neumann on bounding box


    a((u,û,λ),(v,v̂,μ)) = ∫(κ*∇(u)⋅∇(v))dΩ + ∫(λ*(v-v̂))dΓ + ∫(κ̂*∇(û)⋅∇(v̂))dΓ + ∫((u-û)*μ)dΓ

    l((v,v̂,μ)) = ∫(f*v)dΩ + ∫(f̂*v̂)dΓ + ∫(g*μ)dΓ + ∫((h⋅ν)*v)d∂Ω

    # Build affine FE operator
    op = AffineFEOperator(a,l,X,Y)

    # Solve
    ls = LUSolver()
    solver = LinearFESolver(ls)
    wh = solve(solver, op)
    uh, ûh, λh = wh

    if write
        Gridap.writevtk(Ω, "../Data/foo", cellfields=["uh" => uh])
        Gridap.writevtk(Γ, "../Data/bar", cellfields=["uhath" => ûh, "ph"=>λh])
    end
    return uh, ûh, λh
end
