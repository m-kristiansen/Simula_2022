using Gridap
using GridapGmsh
"""
We want to solve the coupled 2D/1D problem
-Δu = f                on Ω
-Δû - λ = f̂            on Γ       λ = -∇u⋅ν
u-û = g                on Γ

"""
function solver(filename,f, f̂, g; write::Bool=true)
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
    V̂ = TestFESpace(Γ,reffe_û,dirichlet_tags = ["starting_point", "tips"], conformity=:H1)
    M = TestFESpace(Γ,reffe_μ, conformity=:L2) #create test space for lagrange multiplier

    Y = MultiFieldFESpace([V,V̂,M])

    dc(x) = 1/(1+0.01*sqrt((x[1]/12)^2+(x[2]/9)^2))
    # Define global trial space
    U = TrialFESpace(V) #include u to enforce dirichlet tradtionally
    Û = TrialFESpace(V̂, [1, dc]) # 1 dirichlet on start of graph
    Λ = TrialFESpace(M) #trial space for lagrange multiplier
    X = MultiFieldFESpace([U,Û,Λ])

    #integration measure
    degree = order
    dΩ = Measure(Ω,degree)
    dΓ = Measure(Γ,degree)
    d∂Ω = Measure(∂Ω, degree)
    ν = get_normal_vector(∂Ω)

    h(x) = VectorValue(0,0) #Neumann on bounding box


    a((u,û,λ),(v,v̂,μ)) = ∫(∇(u)⋅∇(v))dΩ + ∫(λ*(v-v̂))dΓ + ∫(∇(û)⋅∇(v̂))dΓ + ∫((u-û)*μ)dΓ

    l((v,v̂,μ)) = ∫(f*v)dΩ + ∫(f̂*v̂)dΓ + ∫(g*μ)dΓ + ∫((h⋅ν)*v)d∂Ω

    # Build affine FE operator
    op = AffineFEOperator(a,l,X,Y)

    # Solve
    ls = LUSolver()
    solver = LinearFESolver(ls)
    wh = solve(solver, op)
    uh, ûh, λh = wh

    if write
        Gridap.writevtk(Ω, "foo", cellfields=["uh" => uh])
        Gridap.writevtk(Γ, "bar", cellfields=["uhath" => ûh, "ph"=>λh])
    end
    return uh, ûh, λh
end
