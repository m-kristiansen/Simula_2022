using Gridap
using GridapGmsh
using LinearAlgebra: I

"""
Performing a convergence test, where we are going to use all the new concepts we have learned.
We will consider a manufactured solution that does not belong to the FE interpolation space.
In this test, we expect to see the optimal convergence order of the FE discretization.
"""


# Here, we define the manufactured functions
μ = 1E0 #Dynamic viscosity
u(x) = VectorValue(-cos(x[1])*sin(x[2]), cos(x[2])*sin(x[1]))
∇u(x) = TensorValue(sin(x[1])*sin(x[2]), -cos(x[1])*cos(x[2]), cos(x[1])*cos(x[2]), -sin(x[1])*sin(x[2]))
f(x) = -2*μ*VectorValue(cos(x[1])*sin(x[2]), -cos(x[2])*sin(x[1]))

#h(x) = VectorValue(sin(x[2])*sin(x[1]), 0)


import Gridap: ∇
∇(::typeof(u)) = ∇u #get grad u from a

function run(n)
"""
In order to perform the convergence test, we write in a function all the
code needed to perform a single computation and measure its error.
The input of this function is the number of cells in each direction and the interpolation order.
The output is the computed L^2 norm.
"""

    domain = (0,1,0,1)
    partition = (32,32)
    model = CartesianDiscreteModel(domain,partition)#compute mesh and bouding box
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"bottom",[1,2,5])
    add_tag_from_tags!(labels,"left",[7])
    add_tag_from_tags!(labels,"right",[8])
    add_tag_from_tags!(labels,"top",[3,6,4])
    add_tag_from_tags!(labels,"inside",[9])



    # Define reference FE (Q2/P1(disc) pair)
    order = 2
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)

    # Define test FESpaces
    V = TestFESpace(model,reffeᵤ,dirichlet_tags=["bottom", "left", "top"],conformity=:H1)
    Q = TestFESpace(model,reffeₚ,conformity=:L2, constraint=:zeromean)
    Y = MultiFieldFESpace([V,Q])

    # Define trial FESpaces
    U = TrialFESpace(V, u)
    P = TrialFESpace(Q)
    X = MultiFieldFESpace([U,P])

    # Define triangulation and integration measure
    degree = order
    Ωₕ = Triangulation(model)
    dΩ = Measure(Ωₕ,degree)

    #Define Neumann boundary
    Γ = BoundaryTriangulation(model, tags=["right"])
    dΓ = Measure(Γ,degree)
    n_Γ = get_normal_vector(Γ)

    # Define bilinear and linear form
    ε(u) = 0.5*(∇(u)+∇(u)')
    σ(u, p) = 2*μ*ε(u)-p*I
    a((u,p),(v,q)) = ∫(2*μ*ε(u)⊙ε(v) - (∇⋅v)*p + q*(∇⋅u) )dΩ
    l((u,p),(v,q)) = ∫( f⋅v )dΩ + ∫((σ(u, p)⋅n_Γ)⋅v )dΓ

    # Build affine FE operator
    op = AffineFEOperator(a,l,X,Y)

    # Solve
    uh, ph = solve(op)

    eu = u - uh
    println(sqrt(sum( ∫( eu⋅eu )dΩ )))

    p = 0 #pressure
    ep = p - ph

    el = sqrt(sum( ∫( eu⋅eu )dΩ )) + sqrt(sum( ∫( ep*ep )dΩ ))

    return el

end

function conv_test(ns)
"""
The following function does the convergence test.
It takes a vector of integers (representing the number of cells per direction in each computation) plus the interpolation order.
It returns the L^2 error norm for each computation as well as the corresponding cell size.
"""

    els = Float64[]
    hs = Float64[]

    for n in ns

        el = run(n)
        h = 1.0/n

        push!(els,el)
        push!(hs,h)

    end

    return (els, hs)

end

els, hs = conv_test([8,16,32,64,128])

using Plots

plot(hs,[els],
    xaxis=:log, yaxis=:log,
    label=["L2 k=2"],
    shape=:auto,
    xlabel="h",ylabel="error norm")


#get slope
function slope(hs,errors)
  x = log10.(hs)
  y = log10.(errors)
  linreg = hcat(fill!(similar(x), 1), x) \ y
  linreg[2]
end

# The slopes for the $L^2$ error norm is computed as
slope(hs,els)
