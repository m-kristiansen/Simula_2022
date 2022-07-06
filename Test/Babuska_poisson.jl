using Gridap
using GridapGmsh
"""
We want to solve the possion equation
-Δu + u = f on Ω
u = g       on Γd
∇u⋅ν = h    on Γn

with a veraity of mixed Neumann-Dirchlet boundary conditions enforced trough Lagrangian multipliers.
Hence we want to find w, λ such that we minimize
F(w) = ∫(∇w)^2 +w^2 -2∫fw

with the constraint
-2λ∫(w-g)d(∂Ω) = 0 i.e w=g on the boundary

Considering the pertubation w = u+εv and setting dF(u+εv)/dε = 0 we obtain the weak formulation
∫∇u⋅∇v + uv dΩ - ∫λv+μu dΓd = ∫fv dΩ - ∫μg dΓd + ∫hv dΓn
"""

# Here, we define the manufactured functions
# u(x) = x[1]^2 + x[2]^2
# f(x) = x[1]^2 + x[2]^2-4
# ∇u(x) = VectorValue(2*x[1],2*x[2])
# g(x) = u(x)
# λ(x) = 2

u(x) = x[1] + x[2]
f(x) = x[1] + x[2]
∇u(x) = VectorValue(1,1)
g(x) = u(x)
λ(x) = 2

import Gridap: ∇
∇(::typeof(u)) = ∇u #get grad u from here

domain = (0,1,0,1)
partition = (16,16)
model = CartesianDiscreteModel(domain,partition)
model = simplexify(model)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"bottom",[1,5,2])
add_tag_from_tags!(labels,"left",[1,7,3])
add_tag_from_tags!(labels,"right",[4,8,2])
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

order = 1 #Choose order of local interpolation, order 0 => piecewise linear, order 1 => hat function etc..

#Define local basis functions
reffe_u = ReferenceFE(lagrangian, Float64, 2)
reffe_μ = ReferenceFE(lagrangian,Float64, 0)
"""
By changing the order and conformity of M (test space for μ) the conflicting dirichlet conditions i.e.
The babuska and traditional, was solved.
"""

#where should our solution live
Ω = Triangulation(model)

#case 1
#Works with partition = (4,4)
#Γd = BoundaryTriangulation(model, tags=["left", "right", "top", "bottom"])
#Γn = BoundaryTriangulation(model, tags=["bottom", "top"]) #Neumann domain
#V = TestFESpace(Ω,reffe_u,conformity=:H1)

#case 2
#needs partition = (64,64)
Γd = BoundaryTriangulation(model, tags=["left", "right"])
#Γn = BoundaryTriangulation(model, tags=["bottom", "top"]) #Neumann domain
V = TestFESpace(Ω,reffe_u,dirichlet_tags = ["top", "bottom"], conformity=:H1)

#case 3
#works with partition = (8,8)
#Γd = BoundaryTriangulation(model, tags=["left"])
#Γn = BoundaryTriangulation(model, tags=["bottom", "top"]) #Neumann domain
#V = TestFESpace(Ω,reffe_u,dirichlet_tags = ["right"], conformity=:H1)


M = TestFESpace(Γd,reffe_μ, conformity=:L2) #create test space for lagrange multiplier
Y = MultiFieldFESpace([V,M])

# Define global trial space
U = TrialFESpace(V,u) #include u to enforce dirichlet tradtionally
Λ = TrialFESpace(M) #trial space for lagrange multiplier
X = MultiFieldFESpace([U,Λ])


#integration measure
degree = 2
dΩ = Measure(Ω,degree)
dΓd = Measure(Γd,degree)
dΓn = Measure(Γn,degree)
ν = get_normal_vector(Γn)

#Weak formulation
a((u,λ),(v,μ)) = ∫(∇(u)⋅∇(v) + u*v)dΩ - ∫(λ*v+μ*u)dΓd
l((v,μ)) = ∫(f*v)dΩ - ∫(μ*g)dΓd #+ ∫((∇u⋅ν)*v)dΓn

# Build affine FE operator
op = AffineFEOperator(a,l,X,Y)

# Solve
uh, λh = solve(op)

eu = u - uh
println("error_u: ", sqrt(sum( ∫( eu*eu )dΩ )))

#whole boundary
Γ = BoundaryTriangulation(model)
dΓ = Measure(Γ, order)


#define λ = (∇u⋅ν)
foo(x) = x[1] < 0.5 ? -1. : 1.

eλ = foo - λh
println("error_λ: ", sqrt(sum( ∫( eλ*eλ )dΓ )))

writevtk(Ω, "xxx", order=1, cellfields=["uh"=>uh, "eu" => eu])
writevtk(Γd, "yyy", order=1, cellfields=["λh"=>λh, "eλ" =>eλ])

using Plots
plot(λh.free_values.vector, label = "λ")
