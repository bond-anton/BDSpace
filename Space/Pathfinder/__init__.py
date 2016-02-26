import numpy as np
from Space.Coordinates import Cartesian
from Space.Coordinates.transforms import unit_vector
from Space.Curve.Parametric import Helix, Arc


def helix_between_two_points(CS, point1, point2, r=1, loops=1, right=True):
    distance = np.sqrt(np.dot(point2 - point1, point2 - point1))
    direction = unit_vector(point2 - point1)
    print direction
    origin = CS.to_parent(point1)
    helix_CS = Cartesian(basis=np.copy(CS.basis), origin=np.copy(origin), name='Helix CS')
    R = np.sqrt(direction[0]**2 + direction[1]**2)
    theta = np.arcsin(R)
    if R == 0:
        phi = 0
    else:
        phi = np.arccos(direction[0] / R)
    print R, np.rad2deg(theta), np.rad2deg(phi)
    helix_CS.rotate_axis_angle(helix_CS.basis[:, 2], -phi)
    helix_CS.rotate_axis_angle(helix_CS.basis[:, 1], -theta)
    h = distance / int(loops)
    path = Helix(name='Right Helix', CS=helix_CS, r=r, h=h, start=0, stop=np.pi * 2 * int(loops), right=right)
    return path


def arc_between_two_points(CS, point1, point2, r=1, right=True):
    direction = unit_vector(point2 - point1)
    origin = CS.to_parent(point1)
    arc_CS = Cartesian(basis=np.copy(CS.basis), origin=np.copy(origin), name='Arc CS')
    R = np.sqrt(direction[0]**2 + direction[1]**2)
    theta = np.arcsin(R)
    if R == 0:
        phi = 0
    else:
        phi = np.arccos(direction[0] / R)
    arc_CS.rotate_axis_angle(arc_CS.basis[:, 2], -phi)
    arc_CS.rotate_axis_angle(arc_CS.basis[:, 1], -theta)
    local_point2 = arc_CS.to_local(CS.to_parent(point2))
    origin2 = np.array([np.sqrt((2 * r)**2 - local_point2[2]**2), 0, 0])
    print 'Arc CS origin offset:', origin2
    distance = np.sqrt(np.dot(local_point2 - origin2, local_point2 - origin2))
    print distance
    direction = unit_vector(local_point2 - origin2)
    print direction
    origin = arc_CS.to_parent(origin2)
    arc_CS = Cartesian(basis=np.copy(arc_CS.basis), origin=np.copy(origin), name='Arc CS')
    R = np.sqrt(direction[0]**2 + direction[1]**2)
    theta = np.arcsin(R)
    if R == 0:
        phi = 0
    else:
        phi = np.arccos(direction[0] / R)
    print R, np.rad2deg(theta), np.rad2deg(phi)
    arc_CS.rotate_axis_angle(arc_CS.basis[:, 2], -phi)
    arc_CS.rotate_axis_angle(arc_CS.basis[:, 1], -theta)
    local_point1 = arc_CS.to_local(CS.to_parent(point1)) - np.array([0, 0, r])
    local_point2 = arc_CS.to_local(CS.to_parent(point2)) - np.array([0, 0, r])
    print local_point1
    print local_point2
    start = np.pi - np.arccos(local_point1[2] / r)
    stop = np.pi - np.arccos(local_point2[2] / r)
    print start, stop
    path = Arc(CS=arc_CS, r=r, start=start, stop=stop, right=True)
    return path
