using Gridap
using GridapGmsh
"""
We want to solve the coupled 2D/1D problem
-Δu = f                on Ω
-Δû - λ = f̂            on Γ       λ = -∇u⋅ν
u-û = g                on Γ


For simplicity we use the traditional dirchlet boundary conditions with
u = 0 on ∂Ω / Γ
û = 0 on ∂Γ
"""
# Here, we define the manufactured functions
u(x) = sin(x[1])+cos(x[2])
f(x) = u(x)
∇u(x) = VectorValue(cos(x[1]), -sin(x[2]))
Δu(x) = u(x)

û(x) = sin(x[1])+cos(x[2])
∇û(x) = VectorValue(cos(x[1]), -sin(x[2]))
g(x) = u(x)-û(x)
#f̂ is defined after normal vector is computed

import Gridap: ∇
∇(::typeof(u)) = ∇u #get grad u from here


model = GmshDiscreteModel("curved.msh")
#model = simplexify(model)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"Γ0",[1])
add_tag_from_tags!(labels,"Γ1",[2])

order = 2 #Choose order of local interpolation, order 0 => piecewise linear, order 1 => hat function etc..

#Define local basis functions
reffe_u = ReferenceFE(lagrangian, Float64, order)
reffe_û = ReferenceFE(lagrangian, Float64, order)
reffe_μ = ReferenceFE(lagrangian,Float64, order-2)

#where should our solution live
Ω = Triangulation(model)
Γ = BoundaryTriangulation(model, tags = ["bottom"])

V = TestFESpace(Ω,reffe_u,dirichlet_tags = ["right", "left", "top"], conformity=:H1)
V̂ = TestFESpace(Γ,reffe_û,dirichlet_tags = ["Γ0","Γ1"], conformity=:H1)
M = TestFESpace(Γ,reffe_μ, conformity=:L2) #create test space for lagrange multiplier

Y = MultiFieldFESpace([V,V̂,M])

# Define global trial space
U = TrialFESpace(V,u) #include u to enforce dirichlet tradtionally
Û = TrialFESpace(V̂,û)
Λ = TrialFESpace(M) #trial space for lagrange multiplier
X = MultiFieldFESpace([U,Û,Λ])

#integration measure
degree = order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)
ν = get_normal_vector(Γ)

a((u,û,λ),(v,v̂,μ)) = ∫(∇(u)⋅∇(v))dΩ + ∫(λ*(v-v̂))dΓ + ∫(∇(û)⋅∇(v̂))dΓ + ∫((u-û)*μ)dΓ

l((v,v̂,μ)) = ∫(f*v)dΩ + ∫((u+∇u⋅ν)*v̂)dΓ + ∫(g*μ)dΓ #here we have used f̂ = -Δu-λ

# Build affine FE operator
op = AffineFEOperator(a,l,X,Y)

# Solve
ls = LUSolver()
solver = LinearFESolver(ls)
wh = solve(solver, op)
uh, ûh, λh = wh

eu = u - uh
println("error_u: ", sqrt(sum( ∫( eu*eu )dΩ )))

eû = û-ûh
println("error_û: ", sqrt(sum( ∫( eû*eû )dΓ )))

λ = ∇u⋅(-ν)
eλ = λ-λh

println("error_λ: ", sqrt(sum( ∫( eλ⋅eλ )dΓ )))

writevtk(Ω, "foo", cellfields=["uh" => uh])
writevtk(Γ, "bar", cellfields=["uhath" => ûh, "ph"=>λh])
writevtk(Γ,"normal",cellfields=["ν"=>ν])
