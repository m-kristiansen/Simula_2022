using Gridap
using GridapGmsh

model = GmshDiscreteModel("t1.msh")

order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,reffe,dirichlet_tags=["Sides"])

U = TrialFESpace(V,4)
Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)
a(u,v) = ∫( ∇(v)⋅∇(u) )dΩ
l(v) = 0
op = AffineFEOperator(a,l,U,V)
uh = solve(op)
writevtk(Ω,"demo", cellfields=["uh"=>uh])
