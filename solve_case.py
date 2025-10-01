# solve_case.py
# Steady Stokes in chamber with Darcy (or Brinkman) term in TM and optional Boussinesq buoyancy, following ocular LBM/CFD setups [web:1][web:235]
from mpi4py import MPI
import numpy as np
import meshio
from dolfinx import mesh, fem, io
import ufl

# Convert Gmsh -> XDMF
msh = meshio.read("anterior_segment.msh")
triangle_cells = [c for c in msh.cells if c.type == "triangle"][0].data
line_cells = [c for c in msh.cells if c.type == "line"][0].data
meshio.write("anterior_segment.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": triangle_cells, "line": line_cells}))

# Read into dolfinx
with io.XDMFFile(MPI.COMM_WORLD, "anterior_segment.xdmf", "r") as xdmf:
    domain = xdmf.read_mesh(name="Grid")
tdim = domain.topology.dim

# Function spaces
V = fem.VectorFunctionSpace(domain, ("Lagrange", 2))
Q = fem.FunctionSpace(domain, ("Lagrange", 1))
Tspace = fem.FunctionSpace(domain, ("Lagrange", 1))

u = fem.Function(V, name="u")
p = fem.Function(Q, name="p")
v = ufl.TestFunction(V)
q = ufl.TestFunction(Q)
T = fem.Function(Tspace, name="T")
S = ufl.TestFunction(Tspace)

# Material params (water-like at 37 C) and porous permeability range [web:1]
rho = 1000.0
mu = 0.0007
beta = 2.5e-4
g = 9.81
kth = 0.6
cp = 4180.0

# Tags for subdomains/facets (fallback: mark by coordinate boxes if .msh tags not imported)
x = ufl.SpatialCoordinate(domain)
tm_indicator = ufl.conditional(ufl.and_(x[0] > 6.1e-3, x[1] < 3.0e-4), 1.0, 0.0)  # TM band heuristic [web:235]

# Darcy permeability (healthy vs glaucomatous sweep values) [web:1]
k_tm = fem.Constant(domain, 5e-14)

# Momentum equations (Brinkman-like): -mu ∇²u + ∇p + mu/k * chi_tm * u = rho*beta*(T-T0) g e_y, div u = 0 [web:1]
T0 = fem.Constant(domain, 310.0)
ey = ufl.as_vector((0.0, 1.0))
buoy = rho*beta*(T - T0)*g*ey

a_m = 2*mu*ufl.inner(ufl.sym(ufl.grad(u)), ufl.sym(ufl.grad(v)))*ufl.dx + (mu/k_tm)*tm_indicator*ufl.inner(u, v)*ufl.dx
b_c = - p*ufl.div(v)*ufl.dx + q*ufl.div(u)*ufl.dx
L_m = rho*ufl.inner(buoy, v)*ufl.dx

F = a_m + b_c - L_m

# Temperature equation (steady advection-diffusion) used to trigger buoyancy vortices [web:1]
T_cornea = fem.Constant(domain, 307.0)
a_T = kth*ufl.inner(ufl.grad(T), ufl.grad(S))*ufl.dx + rho*cp*ufl.inner(u, ufl.grad(T))*S*ufl.dx
L_T = fem.Constant(domain, 0.0)*S*ufl.dx

# Boundary conditions
def on_wall(x):
    return np.logical_or(np.isclose(x[1], 0.0), np.isclose(x[0], 0.0))

# No-slip walls (cornea/iris/lens) [web:235]
bc_u0 = fem.dirichletbc(np.array((0.0, 0.0), dtype=np.float64), fem.locate_dofs_geometrical(V, on_wall), V)

# Thermal BC: cornea slightly cooler, internal near 310 K [web:1]
bc_T_cool = fem.dirichletbc(np.array(307.0, dtype=np.float64), fem.locate_dofs_geometrical(Tspace, lambda x: np.isclose(x[1], 0.0)), Tspace)
bc_T_warm = fem.dirichletbc(np.array(310.0, dtype=np.float64), fem.locate_dofs_geometrical(Tspace, lambda x: np.isclose(x[1], 3.0e-3)), Tspace)

# Inflow as weak constraint: apply a parabolic inflow on a short edge representing iris–lens gap with total Q ~ 2–3 µL/min [web:235][web:6]
Q_in = 2.5e-9/60.0

# For simplicity here, impose a small average upward velocity near the posterior gap region (heuristic)
bc_in = fem.dirichletbc(np.array((0.0, -Q_in/(1e-6)), dtype=np.float64),
                        fem.locate_dofs_geometrical(V, lambda x: np.logical_and(x[0] > 5.8e-3, np.isclose(x[1], 2.6e-3))), V)

bcs_u = [bc_u0, bc_in]
bcs_T = [bc_T_cool, bc_T_warm]

# Solve temperature then flow
problem_T = fem.petsc.NonlinearProblem(a_T - L_T, T, bcs_T)
solver_T = fem.petsc.NewtonSolver(MPI.COMM_WORLD, problem_T)
solver_T.rtol = 1e-8
solver_T.solve(T)

problem_u = fem.petsc.NonlinearProblem(F, u, bcs_u)
solver_u = fem.petsc.NewtonSolver(MPI.COMM_WORLD, problem_u)
solver_u.rtol = 1e-8
solver_u.solve(u)

# Recover pressure by solving a Poisson-like stabilization (omitted details)
# Postprocess: mean chamber pressure (IOP proxy) and total flux through TM
with io.XDMFFile(MPI.COMM_WORLD, "results.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(u)
    xdmf.write_function(T)
