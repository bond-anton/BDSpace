import numpy as np
from mayavi import mlab

from Space.Coordinates import Cartesian
from Space.Curve.Parametric import Helix
from Space.Pathfinder import helix_between_two_points, arc_between_two_points
from Space import Visual

CS = Cartesian()
CS.rotate_axis_angle(np.ones(3), np.deg2rad(45))
print CS

fig = mlab.figure('CS demo', bgcolor=(0, 0, 0))
Visual.draw_CS_axes(fig, CS)
# Visual.draw_CS_box(fig, CS)

right_helix = Helix(name='Right Helix', CS=CS, r=2, h=0.5, start=0, stop=np.pi * 2, right=True)
left_helix = Helix(name='Left Helix', CS=CS, r=2, h=0.5, start=0, stop=np.pi * 2, right=False)

print right_helix, left_helix
print right_helix.length(symbolic=True), left_helix.length()

# t = np.linspace(0, np.pi*2, num=101, endpoint=True)
# points = right_helix.points(t)
# global_points = right_helix.CS.to_global(points)
# mlab.plot3d(global_points[:, 0], global_points[:, 1], global_points[:, 2], color=(1, 0, 0))

# points = left_helix.points(t)
# global_points = left_helix.CS.to_global(points)
# mlab.plot3d(global_points[:, 0], global_points[:, 1], global_points[:, 2], color=(0, 0, 1))


point1 = np.array([1, 1, 0])
point2 = np.array([2, 2, 0])
points = np.vstack((CS.to_global(point1), CS.to_global(point2)))
mlab.points3d(points[:, 0], points[:, 1], points[:, 2], scale_factor=0.1)
path = helix_between_two_points(CS, point1, point2, r=1, loops=7, right=True)
t = np.linspace(path.start, path.stop, num=101 * (path.stop - path.start) / (2 * np.pi), endpoint=True)
points = path.points(t)
global_points = path.CS.to_global(points)
mlab.plot3d(global_points[:, 0], global_points[:, 1], global_points[:, 2], color=(1, 0, 0))
Visual.draw_CS_axes(fig, path.CS)

path = arc_between_two_points(CS, point1, point2, r=1, right=True)

t = np.linspace(path.start, path.stop, num=101 * (path.stop - path.start) / (2 * np.pi), endpoint=True)
points = path.points(t)
global_points = path.CS.to_global(points)
mlab.plot3d(global_points[:, 0], global_points[:, 1], global_points[:, 2], color=(0, 1, 0))
Visual.draw_CS_axes(fig, path.CS)

mlab.show()
