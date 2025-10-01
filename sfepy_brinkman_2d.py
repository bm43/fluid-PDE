# sfepy_brinkman_2d.py
# Conceptual 2D steady Stokes–Brinkman in SfePy with EVP outlet and porous angle [web:235][web:225]
from sfepy import data_dir
from sfepy.discrete.fem import Mesh, Domain, Field
from sfepy.discrete import Problem
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson

# Load gmsh mesh and define subdomains similarly, then set up Brinkman momentum and continuity equations
# Momentum: -mu ∇²u + ∇p + mu/K u = 0 in TM; -mu ∇²u + ∇p = 0 in chamber; div u = 0 everywhere [web:29]
# BCs: no-slip walls, inflow flux at iris–lens gap, pressure=EVP at TM outlet [web:225][web:235]
# See SfePy examples for 'navier_stokes' and modify for Brinkman; full listing omitted for brevity.
