# build_mesh.py
# Geometry inspired by AS-OCT-based anterior segment models with an annular TM band [web:235]
import gmsh, numpy as np

gmsh.initialize()
gmsh.model.add("anterior_segment")

# Parameters (meters)
R = 6.4e-3         # corneal/limbal radius ~ 12.8 mm diameter / 2 [web:235]
H = 3.0e-3         # anterior chamber depth ~ 3 mm [web:235]
tm_th = 2.0e-4     # TM thickness 0.2 mm idealized band [web:235]
gap = 5.0e-5       # iris–lens gap (inlet aperture) tens of microns [web:235]

# Construct a 2D axisymmetric meridional profile; revolve later if needed
# Cornea arc (simple circle sag)
c1 = gmsh.model.occ.addCircle(0, 0, 0, R)
# Iris chord and lens edge as simple arcs/lines
p_iris_a = gmsh.model.occ.addPoint(0.0, 0.4e-3, 0)
p_iris_b = gmsh.model.occ.addPoint(R-0.5e-3, 0.4e-3, 0)
l_iris = gmsh.model.occ.addLine(p_iris_a, p_iris_b)

# Chamber floor (lens/anterior face simplified)
p_lens_b = gmsh.model.occ.addPoint(R-0.5e-3, H-0.4e-3, 0)
l_side = gmsh.model.occ.addLine(p_iris_b, p_lens_b)
p_lens_a = gmsh.model.occ.addPoint(0.0, H-0.4e-3, 0)
l_lens = gmsh.model.occ.addLine(p_lens_b, p_lens_a)

# Close outer boundary to form a cavity polygon beneath cornea
p_c_left = gmsh.model.occ.addPoint(0, 0, 0)
p_c_right = gmsh.model.occ.addPoint(R, 0, 0)
l_base = gmsh.model.occ.addLine(p_c_left, p_c_right)
l_left = gmsh.model.occ.addLine(p_lens_a, p_c_left)
l_right = gmsh.model.occ.addLine(p_c_right, p_lens_b)

# For simplicity, approximate cornea by a quadratic spline from (0,0) to (R,0) with apex sag
p_cornea_apex = gmsh.model.occ.addPoint(R*0.5, -0.6e-3, 0)
spl_cornea = gmsh.model.occ.addSpline([p_c_left, p_cornea_apex, p_c_right])

# Make a closed curve loop of the chamber boundary
cloop = gmsh.model.occ.addCurveLoop([spl_cornea, l_right, l_lens, l_left])
face = gmsh.model.occ.addPlaneSurface([cloop])

# Mark a TM sub-region near the angle as a thin porous strip
# Here create a small rectangular inclusion near the right limbus
tm_x0 = R-0.3e-3
tm = gmsh.model.occ.addRectangle(tm_x0, 0.05e-3, 0, 0.25e-3, tm_th)

# Cut TM from chamber to create a separate surface for tagging
gmsh.model.occ.fragment([(2, face)], [(2, tm)])

gmsh.model.occ.synchronize()

# Physical groups for simulation
# 2D surfaces
entities = gmsh.model.getEntities(dim=2)
for dim, tag in entities:
    com = gmsh.model.occ.getCenterOfMass(dim, tag)
    if com[0] > tm_x0 and com[1] < tm_th + 0.05e-3:
        gmsh.model.addPhysicalGroup(2, [tag], 2)  # TM porous region
    else:
        gmsh.model.addPhysicalGroup(2, [tag], 1)  # fluid chamber

# Facets: inlet on a tiny iris–lens gap segment near (R-0.5 mm, ~0.4 mm)
# Simplify by tagging the short line l_side midpoint as inlet proxy
gmsh.model.addPhysicalGroup(1, [l_side], 11)   # inlet
gmsh.model.addPhysicalGroup(1, [spl_cornea, l_right, l_lens, l_left], 12)  # no-slip walls

gmsh.model.mesh.generate(2)
gmsh.write("anterior_segment.msh")
gmsh.finalize()
