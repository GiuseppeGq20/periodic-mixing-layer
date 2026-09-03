# ------------------------------------------------------------------------------
#
#  Gmsh Python  periodic mixing layer
#
#  Transfinite mesh
#
#              y
#              ^
#              |
#        +-----+-------------------+
#        |     |   gas             |
#        |     |                   |
#      H |     |   mixing layer    |
#        |     |-------------------|--+
#        |     |                   |  | H_INTERFACE
#        |     |   liquid          |  |
#        +-----+-------------------+--+-> x
#             0                    L
#
# ------------------------------------------------------------------------------

import numpy as np
import sys
from pathlib import Path

import gmsh

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mesh_gen.utils import *

# ---- GEOMETRIC DATA ---- #
H = 1
H_INTERFACE = 0.5 * H
L = 1
W = 0.1

h = [0, H_INTERFACE, H]

# ---- MESH SIZING PARAMETERS ---- #
N_POINTS = 20
N_LIQUID = N_POINTS
N_GAS = N_POINTS
N_X = N_POINTS


# ----           --- #
gmsh.initialize()

# 0D entities
p1 = vertical_points(h, 0)
p2 = vertical_points(h, L)
# 1D entities
l1 = vertical_lines(p1)
l2 = vertical_lines(p2)
lh = horizontal_lines(p1, p2)
# lh = np.array(lh)
print(lh)
curve_loop = quad_loops(l1, l2, lh[:-1], lh[1:])

# set transfinite curve
for l_h in lh:
    gmsh.model.geo.mesh.setTransfiniteCurve(l_h, N_X)
# vertical lines
spacing = [N_LIQUID, N_GAS]
for n, l_v in zip(spacing, l1):
    gmsh.model.geo.mesh.setTransfiniteCurve(l_v, n)
for n, l_v in zip(spacing, l2):
    gmsh.model.geo.mesh.setTransfiniteCurve(l_v, n)

# 2D entities

surfaces = []
for c_loop in curve_loop:
    j = gmsh.model.geo.addPlaneSurface([c_loop])
    gmsh.model.geo.mesh.setTransfiniteSurface(
        j, "Alternate")  # structured grid
    gmsh.model.geo.mesh.setRecombine(2, j)  # recombine quads
    surfaces.append(j)

# --- Extrusion --- #
# returns entities: [(dim, tag), ...]
# keeps order of the original definition
dim_tag_list = []
for surf in surfaces:
    dim_tag_list.append((2, surf))
ext = gmsh.model.geo.extrude(
    dim_tag_list,      # surface
    0, 0, W,           # direction
    numElements=[1],   # number of layers (IMPORTANT)
    recombine=True     # hex elements
)
gmsh.model.geo.synchronize()


# --- PHYSICAL GROUPS USING BOUNDING BOXES --- #

def surfaces_in_box(xmin, ymin, zmin, xmax, ymax, zmax, tol=1e-8):
    return [
        tag
        for dim, tag in gmsh.model.getEntitiesInBoundingBox(
            xmin - tol, ymin - tol, zmin - tol,
            xmax + tol, ymax + tol, zmax + tol,
            2
        )
    ]


# Boundary surfaces
front = surfaces_in_box(0, 0, 0, L, H, 0)
back = surfaces_in_box(0, 0, W, L, H, W)

inlet_liquid = surfaces_in_box(0, 0, 0, 0, H_INTERFACE, W)
inlet_gas = surfaces_in_box(0, H_INTERFACE, 0, 0, H, W)
outlet_liquid = surfaces_in_box(L, 0, 0, L, H_INTERFACE, W)
outlet_gas = surfaces_in_box(L, H_INTERFACE, 0, L, H, W)

bottom = surfaces_in_box(0, 0, 0, L, 0, W)
top = surfaces_in_box(0, H, 0, L, H, W)


def add_physical_surf(name, entities):
    if entities:
        gmsh.model.addPhysicalGroup(2, entities, name=name)


add_physical_surf("front", front)
add_physical_surf("back", back)
add_physical_surf("inlet_liquid", inlet_liquid)
add_physical_surf("inlet_gas", inlet_gas)
add_physical_surf("outlet_liquid", outlet_liquid)
add_physical_surf("outlet_gas", outlet_gas)
add_physical_surf("bottom", bottom)
add_physical_surf("top", top)

# Optional internal mixing-layer interface
interface = surfaces_in_box(
    0, H_INTERFACE, 0,
    L, H_INTERFACE, W
)
add_physical_surf("mixing_layer", interface)


# dummy 3D physical group needed by  openfoam's gmshToFoam conversion utility
volumes = gmsh.model.getEntities(3)
gmsh.model.addPhysicalGroup(3, [volume[1] for volume in volumes], -1, name="")

# ----           --- #
gmsh.model.mesh.generate()

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.write("mixing_layer.msh")
gmsh.finalize()

# write geometric data to file
geom_path = Path("geometry_data.txt")

with geom_path.open("w") as f:
    f.write(f"L {L};\n")
    f.write(f"H_INTERFACE {H_INTERFACE};\n")
    f.write(f"H {H};\n")
