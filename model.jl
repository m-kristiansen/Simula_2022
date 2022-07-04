using Gridap
using GridapGmsh
"""
Solving stokes problem
  −2μ∇·ε(u) + ∇p = f
   ∇·u =0

with boundary conditions
   u            = g     on Γd
  (2με(u)-pI)·n = h     on Γn

Variational formulation
  2μ(ε(u), ε(v))_Ω − (∇·v, p)_Ω =(f, v)_Ω + (h, v)_Γn
  (∇·u, q)_Ω = 0

Boundary conditions is set on boundary box edges denoted Bottom, Left, Top, Right
"""

function solve_PDE(mesh, dirichlet_tags, neumann_tags, write::Bool=true)
#load model
    model = GmshDiscreteModel(mesh) #Load model

    # Define reference FE (Q2/P1(disc) pair)
    order = 2
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)

    # Define test FESpaces
    V = TestFESpace(model,reffeᵤ,dirichlet_tags=dirichlet_tags,conformity=:H1)
    Q = TestFESpace(model,reffeₚ,conformity=:L2,constraint=:zeromean)
    Y = MultiFieldFESpace([V,Q])

    # Define trial FESpaces
    U = TrialFESpace(V, VectorValue(0,1)) #Second argument is Dirichlet values
    P = TrialFESpace(Q)
    X = MultiFieldFESpace([U,P])

    # Define triangulation and integration measure
    degree = order
    Ωₕ = Triangulation(model)
    dΩ = Measure(Ωₕ,degree)

    #Define Neumann boundary
    Γ = BoundaryTriangulation(model, tags=neumann_tags)
    dΓ = Measure(Γ,degree)

    f(x) = x
    h(x) = VectorValue(1.0, 1.0)

    μ = 1e-5 #Dynamic viscosity

    # Define bilinear and linear form
    ε(u) = 0.5*(∇(u)+∇(u)')
    a((u,p),(v,q)) = ∫(2*μ*ε(u)⊙ε(v) - (∇⋅v)*p + q*(∇⋅u) )dΩ
    l((v,q)) = ∫( f⋅v )dΩ + ∫( h⋅v )dΓ

    # Build affine FE operator
    op = AffineFEOperator(a,l,X,Y)

    # Solve
    uh, ph = solve(op)

    # Export results to vtk
    if write
        Gridap.writevtk(Ωₕ,"results",order=2,cellfields=["uh"=>uh,"ph"=>ph])
    end

    return uh, ph
end


if PROGRAM_FILE == basename(@__FILE__)
    solve_PDE("graph.msh", ["Top"], ["Bottom", "Left", "Right"])
end
