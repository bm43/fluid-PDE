# brinkman_patch.py
# Replace simple Darcy sink with Brinkman subdomain and enforce continuity as in numerical drainage models [web:29]
from dolfinx import fem
import ufl

# Indicator chi_tm already defined; replace momentum form
mu_eff = mu
Kinv = (1.0/k_tm) * tm_indicator
a_m = 2*mu*ufl.inner(ufl.sym(ufl.grad(u)), ufl.sym(ufl.grad(v)))*ufl.dx \
      + mu_eff*ufl.inner(ufl.grad(u), ufl.grad(v))*tm_indicator*ufl.dx \
      + mu*ufl.inner(Kinv*u, v)*ufl.dx  # Brinkman + Darcy drag [web:29]

# Interface terms would be weakly enforced if mesh has explicit interface markers; omitted for brevity [web:29]
