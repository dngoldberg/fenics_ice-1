from fenics import *
import scipy.interpolate as interp
import numpy as np
import matplotlib.pyplot as plt
import test_domains
import fenics_util as fu
from IPython import embed


#Number of cells in grid
nx = 50;
ny = 400;

#Fenics mesh
Lx = 50e3
Ly = 400e3

mesh = RectangleMesh(Point(0,0), Point(Lx, Ly), nx, ny)
V = FunctionSpace(mesh, 'Lagrange',1)
v = Function(V)
n = V.dim()
d = mesh.geometry().dim()

dof_coordinates = V.tabulate_dof_coordinates()
dof_coordinates.resize((n, d))
dof_x = dof_coordinates[:, 0]
dof_y = dof_coordinates[:, 1]

#Sampling Mesh, identical to Fenics mesh
domain = test_domains.analytical2(Lx, Ly, nx=nx+1,ny=ny+1)
xcoord = domain.x
ycoord = domain.y
xycoord = (xcoord, ycoord)

embed()

#Data is not stored in an ordered manner on the fencis mesh.
#Using interpolation function to get correct grid ordering
bed_interp = interp.RegularGridInterpolator(xycoord, domain.bed)
surf_interp = interp.RegularGridInterpolator(xycoord, domain.surf)
bmelt_interp = interp.RegularGridInterpolator(xycoord, domain.bmelt)
bdrag_interp = interp.RegularGridInterpolator(xycoord, domain.bdrag)


#Coordinates of DOFS of fenics mesh in order data is stored
dof_xy = (dof_x, dof_y)
bed = bed_interp(dof_xy)
surf = surf_interp(dof_xy)
bmelt = bmelt_interp(dof_xy)
bdrag = bdrag_interp(dof_xy)

#Save mesh and data points at coordinates
dd = '../input/analytical2/'

File(''.join([dd,'analytical2_mesh.xml'])) << mesh

v.vector()[:] = bed.flatten()
File(''.join([dd,'analytical2_mesh_bed.xml'])) <<  v

v.vector()[:] = surf.flatten()
File(''.join([dd,'analytical2_mesh_surf.xml'])) <<  v

v.vector()[:] = bmelt.flatten()
File(''.join([dd,'analytical2_mesh_bmelt.xml'])) <<  v

v.vector()[:] = bdrag.flatten()
File(''.join([dd,'analytical2_mesh_bdrag.xml'])) <<  v