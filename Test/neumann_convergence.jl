using Gridap
using GridapGmsh
#using LinearAlgebra: I #may be of use later

"""
    -∇⋅σ(u, p) = f,         σ(u, p) = 2με(u)-pI,       ε(u) = 0.5(∇u+∇u^T)
       ∇⋅u     = 0

      σ⋅ν      = h  on Γn
      u        = g  on Γd

"""


# Here, we define the manufactured functions
μ = 1 #Dynamic viscosity
u(x) = VectorValue(-cos(x[1])*sin(x[2]), cos(x[2])*sin(x[1]))
∇u(x) = TensorValue(sin(x[1])*sin(x[2]), -cos(x[1])*cos(x[2]), cos(x[1])*cos(x[2]), -sin(x[1])*sin(x[2]))
f(x) = -2*μ*VectorValue(cos(x[1])*sin(x[2]), -cos(x[2])*sin(x[1]))

p(x) = 0 #pressure
σ(x) = 2*μ*TensorValue(sin(x[1])*sin(x[2]), 0, 0, -sin(x[1])*sin(x[2]))


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
    partition = (n,n)
    model = CartesianDiscreteModel(domain,partition)#compute mesh and bouding box
    model = simplexify(model)

    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"bottom",[1,5,2])
    add_tag_from_tags!(labels,"left",[1,7,3])
    add_tag_from_tags!(labels,"right",[2,4,8])
    add_tag_from_tags!(labels,"top",[3,6,4])
    add_tag_from_tags!(labels,"inside",[9])
    #Include corner tags for every side.

    """
    3 ------ 6 ------- 4
    |                  |
    |                  |
    7        9         8
    |                  |
    |                  |
    1 ------ 5 ------- 2
    """


    # Define reference FE (Q2/P1(disc) pair)
    order = 2
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    reffeₚ = ReferenceFE(lagrangian,Float64,order-1)

    # Define test FESpaces
    V = TestFESpace(model,reffeᵤ,dirichlet_tags=["bottom", "right"],conformity=:H1)
    Q = TestFESpace(model,reffeₚ,conformity=:H1)
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
    Γ = BoundaryTriangulation(model, tags=["left", "top"])
    dΓ = Measure(Γ,degree)
    nΓ = get_normal_vector(Γ)

    # Define bilinear and linear form
    ε(u) = 0.5*(∇(u)+∇(u)')
    a((u,p),(v,q)) = ∫(2*μ*ε(u)⊙ε(v) - (∇⋅v)*p + q*(∇⋅u) )dΩ
    l((v,q)) = ∫( f⋅v )dΩ+∫((σ⋅nΓ)⋅v )dΓ

    # Build affine FE operator
    op = AffineFEOperator(a,l,X,Y)

    # Solve
    uh, ph = solve(op)

    eu = u - uh
    println("u: ", sqrt(sum( ∫( eu⋅eu )dΩ )))

    ep = p - ph
    println("p: ", sqrt(sum( ∫( ep⋅ep )dΩ )))
    el = sqrt(sum( ∫( eu⋅eu )dΩ )) + sqrt(sum( ∫( ep*ep )dΩ ))

    #writevtk(Ωₕ,"results",order=2,cellfields=["uh"=>uh,"ph"=>ph, "eu"=>eu])

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

els, hs = conv_test([8, 16, 32, 64, 128])

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
